# Local-FDR gating orchestration
#
# Entry points, sample preparation, gate-combination handling, and condition-
# level control flow. Density estimation, smoothing, filtering, threshold
# selection, and output assembly are kept in separate files.

# Calculate local FDR-based cutpoint for each level of bias
#' @keywords internal
.getCpUnsLoc <- function(
  exList,
  .data,
  chnlSettings,
  stage,
  pathProject
) {
  purrr::map(chnlSettings$biasUns, function(bias) {
    .debug("biasUns", bias) # nolint

    exListPrep <- .prepareDataWithBiasAndNoise(
      exList = exList,
      bias = bias,
      excMin = chnlSettings$excMin
    )

    # get gates for given level of bias across gate combination methods
    list("loc" = .getCpUnsLocGateCombn(
      exListOrig = exListPrep[["exListOrig"]],
      exListNoMin = exListPrep[["exListNoMin"]],
      exTblUnsBias = exListPrep[["exTblUnsBias"]],
      chnlSettings = chnlSettings,
      bias = bias,
      pathProject = pathProject,
      stage = stage
    ))
  }) |>
    purrr::flatten()
}

#' @keywords internal
.prepareDataWithBiasAndNoise <- function(
  exList,
  bias,
  excMin
) {
  # Keep the complete original data.
  exListOrig <- .prepareExListWithBiasAndNoise(
    exList = exList,
    ind = names(exList),
    excMin = FALSE,
    bias = 0
  ) |>
    .arrangeSamplesByExpr()

  # Remove minimum-expression values where requested. Row subsetting keeps the
  # ascending order of exListOrig.
  if (isTRUE(excMin)) {
    exListNoMin <- .prepareExListWithBiasAndNoise(
      exList = exListOrig,
      ind = names(exListOrig),
      excMin = TRUE,
      bias = 0
    )
  } else {
    exListNoMin <- exListOrig
  }

  # The biased unstimulated data differ from the already-prepared unstimulated
  # data only by a constant expression shift. Adding a constant preserves its
  # existing ordering.
  exTblUnsBias <- exListNoMin[[1L]]

  chnl <- attr(exTblUnsBias, "chnlCut")
  exTblUnsBias[[chnl]] <- .getCut(exTblUnsBias) + bias

  list(
    "exListOrig" = exListOrig,
    "exListNoMin" = exListNoMin,
    "exTblUnsBias" = exTblUnsBias
  )
}

# Get the unstim-based local fdr-method
# cutpoint for a given bias across gate combination methods
#' @keywords internal
.getCpUnsLocGateCombn <- function(
  exListOrig,
  exListNoMin,
  exTblUnsBias,
  chnlSettings,
  bias,
  pathProject,
  stage
) {
  .debug("getting gateCombn") # nolint
  cp <- list()

  # gate using prejoined stim data
  if ("prejoin" %in% chnlSettings$gateCombn) {
    .debug("prejoin") # nolint
    # join marker expression for stim samples and sort into ascending order
    cp[["prejoin"]] <- .getCpUnsLocSample(
      exListOrig = .prepareDataForPrejoinInd(exListOrig),
      exListNoMinStim = .prepareDataForPrejoinInd(exListNoMin)[-1],
      exTblUnsBias = exTblUnsBias,
      chnlSettings = chnlSettings,
      bias = bias,
      stage = stage,
      pathProject = pathProject,
      indStim = names(exListNoMin)[-1],
      exListOrigOutput = exListOrig,
      prejoin = TRUE
    )[["loc"]]
  }

  # gate each sample individually, then combine
  nonPrejoinCombn <- setdiff(chnlSettings$gateCombn, "prejoin")
  if (length(nonPrejoinCombn) > 0L) {
    .debug("non-prejoin") # nolint
    cpUnsListNonjoin <- .getCpUnsLocSample(
      exListOrig = exListOrig,
      exListNoMinStim = exListNoMin[-1],
      exTblUnsBias = exTblUnsBias,
      chnlSettings = chnlSettings,
      bias = bias,
      stage = stage,
      pathProject = pathProject
    )
    .debug("Combining cutpoints") # nolint
    cp <- c(cp, .getCpUnsLocCombineCpWithMeta(
      cp = cpUnsListNonjoin[["loc"]],
      gateCombn = nonPrejoinCombn
    ))
  }

  .debug("done getting gateCombn") # nolint
  list("cp" = cp)
}

#' @keywords internal
.prepareDataForPrejoinInd <- function(exList) {
  indUns <- names(exList)[1]
  indStim <- names(exList)[-1]
  exTblStim <- exList[seq.int(2, length(exList))] |>
    dplyr::bind_rows()
  exTblStim <- exTblStim[order(.getCut(exTblStim)), ]
  list(
    exList[[1]],
    exTblStim
  ) |>
    stats::setNames(c(indUns, .createCombinedIdentifier(indStim)))
}

# Sort each sample's rows into ascending order of the cut channel
#' @keywords internal
.arrangeSamplesByExpr <- function(exList) {
  purrr::map(exList, function(x) {
    cut <- .getCut(x)

    if (is.unsorted(cut)) {
      x[order(cut), , drop = FALSE]
    } else {
      x
    }
  })
}

#' @keywords internal
.getCpUnsLocCombineCpWithMeta <- function(cp, gateCombn) {
  meta <- .getCpUnsLocMetaFromCp(cp)
  locGenerated <- meta$locGenerated %in% TRUE
  stimGenerated <- locGenerated & !(meta$locSource %in% "unstim_summary")
  cpReal <- cp
  cpReal[!stimGenerated] <- NA_real_

  purrr::map(gateCombn, function(gateCombnCurr) {
    if (is.null(gateCombnCurr) || gateCombnCurr %in% c("no", "prejoin")) {
      return(.getCpUnsLocCpAttachMeta(cp, meta))
    }

    if (any(stimGenerated, na.rm = TRUE)) {
      cpCombined <- .combineCp(
        cp = cpReal,
        gateCombn = gateCombnCurr
      )[[1]]
      metaOut <- meta
      stimRow <- !(metaOut$locSource %in% "unstim_summary")
      metaOut$locGenerated[stimRow] <- TRUE
      stimDirectMatches <- (meta$locGeneratedDirect %in% TRUE) &
        (abs(
          suppressWarnings(as.numeric(cp)) -
            suppressWarnings(as.numeric(cpCombined))
        ) < 1e-7)
      metaOut$locGeneratedDirect[stimRow] <- stimDirectMatches[stimRow]
      metaOut$locSource[stimRow & !(metaOut$locGeneratedDirect %in% TRUE)] <-
        "combined"
      metaOut$locReason[stimRow & !(metaOut$locGeneratedDirect %in% TRUE)] <-
        "combined_from_generated_local_fdr_thresholds"
      metaOut$locGenerated[!stimRow] <- any(stimGenerated, na.rm = TRUE)
      metaOut$locGeneratedDirect[!stimRow] <- FALSE
      metaOut$locSource[!stimRow] <- "unstim_summary"
      metaOut$locReason[!stimRow] <-
        "summary_of_combined_generated_local_fdr_thresholds"
      return(.getCpUnsLocCpAttachMeta(cpCombined, metaOut))
    }

    cpCombined <- .combineCp(
      cp = cp,
      gateCombn = gateCombnCurr
    )[[1]]
    metaOut <- meta
    metaOut$locGenerated[] <- FALSE
    metaOut$locGeneratedDirect[] <- FALSE
    metaOut$locSource[] <- "not_calculated"
    metaOut$locReason[] <- "no_generated_local_fdr_threshold_to_combine"
    .getCpUnsLocCpAttachMeta(cpCombined, metaOut)
  }) |>
    stats::setNames(gateCombn)
}

#' @keywords internal
.getCpUnsLocMetaFromCp <- function(cp) {
  n <- length(cp)
  nm <- names(cp)
  if (is.null(nm)) {
    nm <- rep(NA_character_, n)
  }
  tibble::tibble(
    ind = as.character(nm),
    locGenerated = attr(cp, "locGenerated") %||% rep(FALSE, n),
    locGeneratedDirect = attr(cp, "locGeneratedDirect") %||% rep(FALSE, n),
    locSource = attr(cp, "locSource") %||% rep("not_calculated", n),
    locReason = attr(cp, "locReason") %||% rep(NA_character_, n)
  ) |>
    dplyr::mutate(
      locGenerated = .data$locGenerated %in% TRUE,
      locGeneratedDirect = .data$locGeneratedDirect %in% TRUE
    )
}

#' @keywords internal
.getCpUnsLocCpAttachMeta <- function(cp, meta) {
  if (nrow(meta) != length(cp)) {
    stop("Local-FDR metadata length does not match cutpoint vector length")
  }
  attr(cp, "locGenerated") <- meta$locGenerated %in% TRUE
  attr(cp, "locGeneratedDirect") <- meta$locGeneratedDirect %in% TRUE
  attr(cp, "locSource") <- as.character(meta$locSource)
  attr(cp, "locReason") <- as.character(meta$locReason)
  cp
}


# ------------------------------------------
# get cutpoints for a range of samples, and then individual samples
# ------------------------------------------

# Get cutpoint for a range of samples given the q-value and fdr
#' @keywords internal
.getCpUnsLocSample <- function(
  exListOrig,
  exListNoMinStim,
  exTblUnsBias,
  chnlSettings,
  bias,
  pathProject,
  stage,
  indStim = NULL,
  exListOrigOutput = NULL,
  prejoin = FALSE
) {
  .debug("getting loc gate at sample level") # nolint
  indStim <- indStim %||% names(exListNoMinStim)
  exListOrigOutput <- exListOrigOutput %||% exListOrig

  # get cutpoints for each sample
  cpUnsLocObjList <- purrr::map(
    seq_along(exListNoMinStim),
    function(i) {
      .debug("sample", i) # nolint

      exTblNoMinStim <- exListNoMinStim[[i]]
      ind <- .getInd(exTblNoMinStim)
      .debug("ind", ind) # nolint
      exTblUnsOrig <- exListOrig[[1]]
      exTblStimOrig <- exListOrig[[i + 1]]
      chnl <- chnlSettings$chnlCut %||%
        .getCpUnsLocGetChnl(exTblNoMinStim)
      stageChnl <- file.path(stage, chnl)
      .intSave(
        ind,
        stageChnl,
        pathProject,
        exTblNoMinStim,
        exTblUnsOrig,
        exTblStimOrig
      )

      # return early if there are too few cells
      tooFewCellsLgl <- .getCpUnsLocSampleCheckCellNumber(
        exTblStimNoMin = exTblNoMinStim,
        minCell = chnlSettings$minCell,
        exTblUnsBias = exTblUnsBias
      )
      if (tooFewCellsLgl) {
        objOut <- .getCpUnsLocTooFew(
          label = "too_few_cells_sample_fn",
          stage = stage,
          pathProject = pathProject,
          exTblNoMinStim = exTblNoMinStim,
          exTblUnsBias = exTblUnsBias,
          cpMin = chnlSettings$cpMin
        )
        return(objOut)
      }

      # remove any cytokine-positive cells from unstim using gates from
      # sample for which gates are required
      exTblUnsBiasRm <- .getCpUnsLocSampleUnsRmCytPos(
        exTblUnsOrig = exTblUnsOrig,
        chnlSettings = chnlSettings,
        exTblStimNoMin = exTblNoMinStim,
        bias = bias,
        exTblUnsBias = exTblUnsBias,
        stage = stage
      )
      # removal can leave too few unstim cells, so check again
      if (nrow(exTblUnsBiasRm) < chnlSettings$minCell) {
        objOut <- .getCpUnsLocTooFew(
          label = "too_few_cells_sample_fn",
          stage = stage,
          pathProject = pathProject,
          exTblNoMinStim = exTblNoMinStim,
          exTblUnsBias = exTblUnsBias,
          cpMin = chnlSettings$cpMin
        )
        return(objOut)
      }
      exTblUnsBias <- exTblUnsBiasRm
      .intSave(ind, stageChnl, pathProject, exTblUnsBias)

      .getCpUnsLocCondition(
        exTblUnsBias = exTblUnsBias,
        exTblStimNoMin = exTblNoMinStim,
        chnlSettings = chnlSettings,
        exTblStimOrig = exTblStimOrig,
        exTblUnsOrig = exTblUnsOrig,
        bias = bias,
        pathProject = pathProject,
        stage = stage
      )
    }
  ) |>
    stats::setNames(names(exListNoMinStim))

  chnl <- .getCpUnsLocGetChnl(exListNoMinStim[[1]])
  .getCpUnsLocOutput(
    cpUnsLocObjList = cpUnsLocObjList,
    indUns = names(exListOrig)[1],
    indStim = indStim,
    stage = stage,
    pathProject = pathProject,
    chnl = chnl,
    exListOrig = exListOrigOutput,
    prejoin = prejoin
  )
}

#' @keywords internal
.getCpUnsLocSampleCheckCellNumber <- function(
  exTblStimNoMin,
  minCell,
  exTblUnsBias
) {
  nrow(exTblStimNoMin) < minCell ||
    nrow(exTblUnsBias) < minCell
}

# Fallback result when there are too few cells to gate; label names the
# intermediate save marking where the check failed
#' @keywords internal
.getCpUnsLocTooFew <- function(
  label,
  stage,
  pathProject,
  exTblNoMinStim,
  exTblUnsBias,
  cpMin
) {
  chnl <- .getCpUnsLocGetChnl(exTblNoMinStim)
  stageChnl <- file.path(stage, chnl)
  .intSaveNm(
    label,
    NULL,
    .getInd(exTblNoMinStim),
    stageChnl,
    pathProject
  ) # nolint
  objOut <- .getCpUnsLocConditionCheckOut(
    cpMin = cpMin,
    exTblStimNoMin = exTblNoMinStim,
    exTblUnsBias = exTblUnsBias,
    stage = stage,
    msg = "Too few cells"
  )
  .intSaveNm(
    "cpInd",
    objOut$cp,
    .getInd(exTblNoMinStim),
    stageChnl,
    pathProject
  )
  objOut
}

#' @keywords internal
.getCpUnsLocSampleUnsRmCytPos <- function(
  exTblUnsOrig,
  chnlSettings,
  exTblStimNoMin,
  bias,
  exTblUnsBias,
  stage
) {
  if (stage == "init") {
    return(exTblUnsBias)
  }
  .debug("Removing cytokine-positive cells from unstim") # nolint

  # first filter
  gateTblGnInd <- chnlSettings$gateTbl |>
    dplyr::filter(
      ind == exTblStimNoMin$ind[1], # nolint
      .data$gateName == .env$chnlSettings$gateNameCurr # nolint
    )

  posIndVecButSinglePosCurr <-
    .getPosIndButSinglePosForOneCyt(
      ex = exTblUnsOrig,
      gateTbl = gateTblGnInd,
      chnlSingleExc = chnlSettings$chnlCut,
      chnl = NULL,
      gateTypeCytPos = ifelse(
        chnlSettings$calcCytPosGates,
        "cyt",
        "base"
      )
    )

  exTblUnsOrig <- exTblUnsOrig[
    !posIndVecButSinglePosCurr, ,
    drop = FALSE
  ]

  # re-apply bias, noise and exclude minimum after removing cytokine-positive cells
  .prepareExListWithBiasAndNoise(
    exList = stats::setNames(
      list(exTblUnsOrig),
      attr(exTblUnsOrig, "indUns")
    ),
    ind = attr(exTblUnsOrig, "indUns"),
    excMin = chnlSettings$excMin,
    bias = bias
  )[[1]]
}

#' @keywords internal
.getCpUnsLocGetChnl <- function(exTbl) {
  if (is.null(exTbl) || !is.data.frame(exTbl)) {
    return("unknown_chnl")
  }
  chnl <- attr(exTbl, "chnlCut")
  if (!is.null(chnl) && nzchar(chnl)) {
    return(chnl)
  }
  cutVec <- .getCut(exTbl)
  cnVec <- colnames(exTbl)
  if (is.null(cnVec) || length(cnVec) == 0L) {
    return("unknown_chnl")
  }
  matchInd <- vapply(
    cnVec,
    function(nm) identical(exTbl[[nm]], cutVec),
    logical(1)
  )
  if (any(matchInd)) {
    return(cnVec[matchInd][1])
  }
  cnVec[1]
}

#' @keywords internal
.getCpUnsLocCondition <- function(
  exTblUnsBias,
  exTblStimNoMin,
  chnlSettings,
  exTblStimOrig,
  exTblUnsOrig,
  plot = TRUE,
  probMin = 0.1,
  bias,
  pathProject,
  stage
) {
  .debug("getting loc gate for single sample") # nolint
  ind <- .getInd(exTblStimNoMin)
  chnl <- chnlSettings$chnlCut %||%
    .getCpUnsLocGetChnl(exTblStimNoMin)
  stageChnl <- file.path(stage, chnl)
  .debug("ind", ind) # nolint

  # return early if almost all stim expression is below cpMin (cell numbers
  # were already checked by .getCpUnsLocSample)
  cutStim <- .getCut(exTblStimNoMin)
  if (
    (stats::quantile(cutStim, 0.9) + 3 * stats::sd(cutStim)) <=
      chnlSettings$cpMin
  ) {
    objOut <- .getCpUnsLocTooFew(
      label = "too_few_cells_ind_fn",
      stage = stage,
      pathProject = pathProject,
      exTblNoMinStim = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      cpMin = chnlSettings$cpMin
    )
    return(objOut)
  }

  exTblStimThreshold <- exTblStimNoMin
  exTblUnsThreshold <- exTblUnsBias
  .intSave(
    .getInd(exTblStimNoMin),
    stageChnl,
    pathProject,
    exTblStimThreshold,
    exTblUnsThreshold
  )

  # get smoothed probabilities
  dataMod <- .getCpUnsLocGetProb(
    exTblStimNoMin = exTblStimNoMin,
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold,
    exTblUnsBias = exTblUnsBias,
    stage = stage,
    bias = bias,
    exTblUnsOrig = exTblUnsOrig,
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )
  .intSave(
    .getInd(exTblStimNoMin),
    stageChnl,
    pathProject,
    dataMod
  )

  # get threshold
  .getCpUnsLocGetCp(
    dataMod = dataMod,
    exTblStimNoMin = exTblStimNoMin,
    exTblStimOrig = exTblStimOrig,
    exTblUnsOrig = exTblUnsOrig,
    exTblUnsBias = exTblUnsBias,
    bias = bias,
    cpMin = chnlSettings$cpMin,
    stage = stage,
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )
}

#' @keywords internal
.getCpUnsLocConditionCheckOut <- function(
  cpMin,
  exTblStimNoMin,
  exTblUnsBias,
  stage,
  msg
) {
  .debug(msg) # nolint
  .getCpUnsLocConditionOut(
    cp = .getCpUnsLocConditionCpNonLoc(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias
    ),
    locGenerated = FALSE,
    locGeneratedDirect = FALSE,
    locSource = "not_calculated",
    locReason = msg
  )
}


#' @keywords internal
.getCpUnsLocConditionOut <- function(
  cp,
  locGenerated,
  locGeneratedDirect,
  locSource,
  locReason
) {
  list(
    cp = cp,
    locGenerated = isTRUE(locGenerated),
    locGeneratedDirect = isTRUE(locGeneratedDirect),
    locSource = locSource,
    locReason = locReason
  )
}


#' @keywords internal
.getCpUnsLocConditionCpNonLoc <- function(
  cpMin,
  exTblStimNoMin,
  exTblUnsBias
) {
  rangeVecStim <- range(.getCut(exTblStimNoMin))
  rangeVecUns <- range(.getCut(exTblUnsBias))
  rangeLen <- max(diff(rangeVecStim), diff(rangeVecUns))
  max(
    cpMin,
    rangeVecStim[[2]] + rangeLen / 5,
    rangeVecUns[[2]] + rangeLen / 3
  )
}
