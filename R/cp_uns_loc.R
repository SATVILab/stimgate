# Calculate local FDR-based cutpoint
#' @keywords internal
.getCpUnsLoc <- function(
  exList,
  .data,
  chnlSettings,
  stage,
  pathProject
) {
  # get cutpoints for each level of bias
  .getCpUnsLocBias(
    exList = exList,
    .data = .data,
    chnlSettings = chnlSettings,
    stage = stage,
    pathProject = pathProject
  )
}

# Get unstim-based local FDR cutpoint for each bias level
#' @keywords internal
.getCpUnsLocBias <- function(
  exList,
  .data,
  chnlSettings,
  stage,
  pathProject
) {
  # get ecdf of uns
  purrr::map(chnlSettings$biasUns, function(bias) {
    .debug("biasUns", bias) # nolint

    exListPrep <- .prepareDataWithBiasAndNoise(
      exList = exList,
      bias = bias,
      noiseSd = NULL,
      excMin = chnlSettings$excMin
    )

    # get gates for given level of bias across gate combination methods
    # --------------------------------------
    cpUnsGateCombnObj <- .getCpUnsLocGateCombn(
      exListOrig = exListPrep[["exListOrig"]],
      exListNoMin = exListPrep[["exListNoMin"]],
      exTblUnsBias = exListPrep[["exTblUnsBias"]],
      chnlSettings = chnlSettings,
      bias = bias,
      pathProject = pathProject,
      stage = stage
    )

    # extract and add bias label to gates
    # ---------------------------------------
    .getCpUnsLocGateLabel(cpUnsGateCombnObj, bias)
  }) |>
    purrr::flatten()
}

#' @keywords internal
.getCpUnsLocGateLabel <- function(cpUnsGateCombnObj, bias) {
  cpUnsGateCombnList <- cpUnsGateCombnObj[["cpUns"]]
  names(cpUnsGateCombnList) <-
    paste0(names(cpUnsGateCombnList))
  cpUnsGateCombnList
}

#' @keywords internal
.prepareDataWithBiasAndNoise <- function(
  exList,
  bias,
  excMin,
  noiseSd
) {
  # rename `cut` column` to `expr`
  # -------------------------------------
  exListOrig <- .prepareExListWithBiasAndNoise(
    # excMin should always be FALSE here, as we're trying
    # to keep all the original data
    exList = exList,
    ind = names(exList),
    excMin = FALSE,
    bias = 0,
    noiseSd = NULL
  ) |>
    .arrangeSamplesByDescExpr()

  # separate stim and uns samples, rename
  # `cut` column to `expr` and exclude min
  # values
  # -------------------------------------
  exListNoMin <- .prepareExListWithBiasAndNoise(
    exList = exList,
    ind = names(exList),
    excMin = excMin,
    bias = 0,
    noiseSd = NULL
  ) |>
    .arrangeSamplesByDescExpr()

  # adjust expression in unstim,
  # applying bias, excluding the min val
  # and /or adding noise
  # -------------------------------------
  exListBias <- .prepareExListWithBiasAndNoise(
    exList = exList,
    ind = names(exList)[1],
    excMin = excMin,
    bias = bias,
    noiseSd = NULL
  ) |>
    .arrangeSamplesByDescExpr()
  exTblUnsBias <- exListBias[[1]]

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

  # get cutpoints for prejoin gate combination method
  cpUnsListPrejoin <- .getCpUnsLocGateCombnPrejoin(
    exListNoMin = exListNoMin,
    exListOrig = exListOrig,
    exTblUnsBias = exTblUnsBias,
    chnlSettings = chnlSettings,
    bias = bias,
    pathProject = pathProject,
    stage = stage
  )

  # get cutpoint for non-prejoin grouping methods
  cpUnsListPrejoinNon <- .getCpUnsLocGateCombnPrejoinNon(
    nonPrejoinCombn = setdiff(chnlSettings$gateCombn, "prejoin"),
    exListNoMinStim = exListNoMin[-1],
    exListOrig = exListOrig,
    exTblUnsBias = exTblUnsBias,
    chnlSettings = chnlSettings,
    bias = bias,
    pathProject = pathProject,
    stage = stage
  )

  # merge above two lists
  cpUnsList <- .getCpUnsLocGateCombnMerge(
    cpUnsListPrejoin,
    cpUnsListPrejoinNon,
    stage
  )

  list(
    "cpUns" = list("loc" = cpUnsList),
    "pList" = list()
  )
}

#' @keywords internal
.getCpUnsLocGateCombnMerge <- function(
  cpUnsListPrejoin,
  cpUnsListPrejoinNon,
  stage
) {
  .debug("done getting gateCombn") # nolint

  combinedList <- cpUnsListPrejoin |>
    append(cpUnsListPrejoinNon)
  purrr::map(
    unique(names(combinedList)),
    function(x) {
      .debug("cutpoint name", paste0(x, collapse = "-")) # nolint
      cpUnsListPrejoin[[x]] |>
        append(cpUnsListPrejoinNon[[x]])
    }
  ) |>
    stats::setNames(unique(names(combinedList)))
}

# --------------------------------
# gate using prejoined .data
# --------------------------------
#' @keywords internal
.getCpUnsLocGateCombnPrejoin <- function(
  exListNoMin,
  exListOrig,
  exTblUnsBias,
  chnlSettings,
  bias,
  pathProject,
  stage
) {
  if (!"prejoin" %in% chnlSettings$gateCombn) {
    return(.getCpUnsLocGateCombnPrejoinNot())
  }
  .getCpUnsLocGateCombnPrejoinActual(
    exListNoMin = exListNoMin,
    exListOrig = exListOrig,
    exTblUnsBias = exTblUnsBias,
    chnlSettings = chnlSettings,
    bias = bias,
    pathProject = pathProject,
    stage = stage
  )
}

#' @keywords internal
.getCpUnsLocGateCombnPrejoinNot <- function() {
  list("cp" = list(), "pList" = list())
}

#' @keywords internal
.prepareDataForPrejoin <- function(exListOrig, exListNoMin) {
  exListNoMin <- .prepareDataForPrejoinInd(
    exList = exListNoMin
  )
  exListOrig <- .prepareDataForPrejoinInd(
    exList = exListOrig
  )

  list(
    "exListOrig" = exListOrig,
    "exListNoMin" = exListNoMin
  )
}

#' @keywords internal
.prepareDataForPrejoinInd <- function(exList) {
  indUns <- names(exList)[1]
  exTblStim <- exList[seq.int(2, length(exList))] |>
    dplyr::bind_rows()
  exTblStim <- exTblStim[order(.getCut(exTblStim)), ]
  list(
    exList[[1]],
    exTblStim
  ) |>
    stats::setNames(c(indUns, names(exList)[seq(2, length(exList))]))
}

#' @keywords internal
.getCpUnsLocGateCombnPrejoinActual <- function(
  exListNoMin,
  exListOrig,
  exTblUnsBias,
  chnlSettings,
  bias,
  pathProject,
  stage
) {
  .debug("prejoin") # nolint

  # get marker expression for stim samples,
  # join and then sort into descending order
  exListPrejoin <- .prepareDataForPrejoin(
    exListOrig = exListOrig,
    exListNoMin = exListNoMin
  )

  # get cutpoints for gate combn method for a range of fdr's
  .getCpUnsLocSample(
    exListOrig = exListPrejoin[["exListOrig"]],
    exListNoMinStim = exListPrejoin[["exListNoMin"]],
    exTblUnsBias = exTblUnsBias,
    chnlSettings = chnlSettings,
    bias = bias,
    stage = stage,
    pathProject = pathProject
  ) |>
    purrr::map(function(x) list("prejoin" = x))
}

# --------------------------------
# gate each sample individually
# --------------------------------
#' @keywords internal
.getCpUnsLocGateCombnPrejoinNon <- function(
  nonPrejoinCombn,
  exListNoMinStim,
  exListOrig,
  exTblUnsBias,
  chnlSettings,
  bias,
  pathProject,
  stage
) {
  if (length(nonPrejoinCombn) == 0L) {
    return(
      .getCpUnsLocGateCombnPrejoinNonNot()
    )
  }
  .getCpUnsLocGateCombnPrejoinNonActual(
    exListNoMinStim = exListNoMinStim,
    exListOrig = exListOrig,
    exTblUnsBias = exTblUnsBias,
    chnlSettings = chnlSettings,
    bias = bias,
    nonPrejoinCombn = nonPrejoinCombn,
    pathProject = pathProject,
    stage = stage
  )
}

#' @keywords internal
.getCpUnsLocGateCombnPrejoinNonNot <- function() {
  list("cp" = list(), "pList" = list())
}

#' @keywords internal
.getCpUnsLocGateCombnPrejoinNonActual <- function(
  exListNoMinStim,
  exListOrig,
  exTblUnsBias,
  chnlSettings,
  bias,
  pathProject,
  stage,
  nonPrejoinCombn
) {
  .debug("non-prejoin") # nolint
  cpUnsListNonjoin <- .getCpUnsLocSample(
    exListOrig = exListOrig,
    exListNoMinStim = exListNoMinStim,
    exTblUnsBias = exTblUnsBias,
    chnlSettings = chnlSettings,
    bias = bias,
    stage = stage,
    pathProject = pathProject
  )

  cpUnsListNonjoin <- .getCpUnsLocGateCombnPrejoinNonActualCombn(
    stage,
    cpUnsListNonjoin,
    nonPrejoinCombn
  )
  list("cp" = cpUnsListNonjoin, "pList" = list())
}

#' @keywords internal
.arrangeSamplesByDescExpr <- function(exList) {
  # arrange in descending order of expression
  exList |>
    purrr::map(function(x) x[order(.getCut(x)), ]) # nolint
}

#' @keywords internal
.getCpUnsLocGateCombnPrejoinNonActualCombn <- function(
  stage,
  cpUnsListNonjoin, # nolint
  nonPrejoinCombnVec
) {
  .debug("Combining cutpoints") # nolint
  .getCpUnsLocCombineCpWithMeta(
    cp = cpUnsListNonjoin[["loc"]],
    gateCombn = nonPrejoinCombnVec
  )
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
      metaOut$locGeneratedDirect[stimRow] <- meta$locGeneratedDirect[
        stimRow
      ] %in%
        TRUE
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
  stage
) {
  .debug("getting loc gate at sample level") # nolint
  force(pathProject)

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
        objOut <- .getCpUnsLocSampleTooFew(
          stage = stage,
          pathProject = pathProject,
          exTblNoMinStim = exTblNoMinStim,
          exTblUnsBias = exTblUnsBias,
          cpMin = chnlSettings$cpMin
        )
        return(objOut)
      }

      # remove any cytokine-positive cells from unstim using gates from
      # sample for which single-positive gates are required
      exTblUnsBias <- .getCpUnsLocSampleUnsRmCytPos(
        exTblUnsOrig = exTblUnsOrig,
        chnlSettings = chnlSettings,
        exTblStimNoMin = exTblNoMinStim,
        bias = bias,
        exTblUnsBias = exTblUnsBias,
        stage = stage
      )
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

  .pathProject <- pathProject
  chnl <- .getCpUnsLocGetChnl(exListNoMinStim[[1]])
  .getCpUnsLocOutput(
    cpUnsLocObjList = cpUnsLocObjList,
    indUns = names(exListOrig)[1],
    indStim = names(exListNoMinStim),
    stage = stage,
    pathProject = .pathProject,
    chnl = chnl,
    exListOrig = exListOrig
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

#' @keywords internal
.getCpUnsLocSampleTooFew <- function(
  stage,
  pathProject,
  exTblNoMinStim,
  exTblUnsBias,
  cpMin
) {
  chnl <- .getCpUnsLocGetChnl(exTblNoMinStim)
  stageChnl <- file.path(stage, chnl)
  .intSaveNm(
    "too_few_cells_sample_fn",
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
      ),
      gateTypeSinglePos = "base"
    )

  exTblUnsOrig <- exTblUnsOrig[
    !posIndVecButSinglePosCurr,
    ,
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
    bias = bias,
    noiseSd = NULL
  )[[1]]
}

#' @keywords internal
.getCpUnsLocPListEmpty <- function() {
  lapply(seq_len(3), function(x) list()) |>
    stats::setNames(c("pLocDens", "pLocProb", "pLocCtb"))
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

  # estimate densities for stim and unstim over stim range
  if (
    .getCpUnsLocCheckEarly(
      exTblStimNoMin,
      chnlSettings$minCell,
      chnlSettings$cpMin
    )
  ) {
    objOut <- .getCpUnsLocConditionTooFew(
      stage = stage,
      pathProject = pathProject,
      exTblNoMinStim = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      cpMin = chnlSettings$cpMin
    )
    return(objOut)
  }

  # stop expr being higher than maxX to prevent really far away values creating modes
  exTblStimThreshold <- .getCpUnsLocSetMaxExpr(
    exTblStimNoMin,
    .getCpUnsLocConditionMaxDensX(exTblStimNoMin)
  )
  exTblUnsThreshold <- .getCpUnsLocSetMaxExpr(
    exTblUnsBias,
    .getCpUnsLocConditionMaxDensX(exTblStimNoMin)
  )
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

.getCpUnsLocConditionTooFew <- function(
  stage,
  pathProject,
  exTblNoMinStim,
  exTblUnsBias,
  cpMin
) {
  chnl <- .getCpUnsLocGetChnl(exTblNoMinStim)
  stageChnl <- file.path(stage, chnl)
  .intSaveNm(
    "too_few_cells_ind_fn",
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
  return(objOut)
}

# initial checks
# ---------------------
#' @keywords internal
.getCpUnsLocConditionCheckNCell <- function(exTblStimNoMin, minCell) {
  nrow(exTblStimNoMin) < minCell
}

#' @keywords internal
.getCpUnsLocConditionMaxDensX <- function(exTblStimNoMin) {
  max(.getCut(exTblStimNoMin)) -
    0.05 * (diff(range(.getCut(exTblStimNoMin))))
}

#' @keywords internal
.getCpUnsLocConditionCheckMaxX <- function(exTblStimNoMin, cpMin) {
  .getCpUnsLocConditionMaxDensX(exTblStimNoMin) <= cpMin
}

#' @keywords internal
.getCpUnsLocCheckEarly <- function(exTblStimNoMin, minCell, cpMin) {
  .getCpUnsLocConditionCheckNCell(exTblStimNoMin, minCell) ||
    .getCpUnsLocConditionCheckMaxX(exTblStimNoMin, cpMin)
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
  locReason,
  pList = .getCpUnsLocPListEmpty()
) {
  list(
    cp = cp,
    pList = pList,
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

#' @keywords internal
.getCpUnsLocSetMaxExpr <- function(.data, maxX) {
  .data[, attr(.data, "chnlCut")] <- pmin(.getCut(.data), maxX)
  .data
}

# get probabilities
#' @keywords internal
.getCpUnsLocGetProb <- function(
  exTblStimNoMin,
  exTblStimThreshold,
  exTblUnsThreshold,
  exTblUnsBias,
  bias,
  exTblUnsOrig,
  stage,
  pathProject,
  chnlSettings
) {
  ind <- .getInd(exTblStimNoMin)
  chnl <- .getCpUnsLocGetChnl(exTblStimNoMin)
  stageChnl <- file.path(stage, chnl)

  # get raw densities
  densTblRaw <- .getCpUnsLocGetDensRaw(
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold,
    stage = stage,
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )
  .intSave(ind, stageChnl, pathProject, densTblRaw)

  # get probabilities
  probTblList <- .getCpUnsLocGetProbTbl(
    densTblRaw = densTblRaw,
    stage = stage,
    cpMin = chnlSettings$cpMin + bias,
    exVecStimThreshold = .getCut(exTblStimThreshold),
    exVecUnsThreshold = .getCut(exTblUnsThreshold)
  )
  .intSave(ind, stageChnl, pathProject, probTblList)

  # get .data to smooth over
  dataMod <- .getCpUnsLocGetDataMod(
    exTblStimThreshold = exTblStimThreshold,
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsThreshold = exTblUnsThreshold,
    exTblUnsBias = exTblUnsBias,
    probTblList = probTblList,
    cpMin = chnlSettings$cpMin + bias,
    stage = stage
  )
  .intSave(ind, stageChnl, pathProject, dataMod)

  # smooth
  .getCpUnsLocGetProbSmooth(
    dataMod = dataMod,
    stage = stage,
    pathProject = pathProject,
    chnl = chnl,
    chnlSettings = chnlSettings
  )
}

# get densTblRaw
# -----------------------
#' @keywords internal
.getCpUnsLocGetDensRaw <- function(
  exTblStimThreshold,
  exTblUnsThreshold,
  stage,
  pathProject,
  chnlSettings
) {
  .debug("Calculating densities") # nolint

  densList <- .getCpUnsLocGetDensRawDensities(
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold,
    stage = stage,
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )

  # put raw densities into table
  .getCpUnsLocGetDensRawTabulate(
    stimX = densList$stim$x,
    stimY = densList$stim$y,
    unsX = densList$uns$x,
    unsY = densList$uns$y
  )
}

#' @keywords internal
.getCpUnsLocGetDensRawDensities <- function(
  exTblStimThreshold,
  exTblUnsThreshold,
  stage,
  pathProject,
  chnlSettings
) {
  useAdaptive <- .getCpUnsLocUseAdaptiveBw(chnlSettings)

  if (isTRUE(useAdaptive)) {
    densAdaptive <- .getCpUnsLocGetDensRawDensitiesAdaptive(
      exTblStimThreshold = exTblStimThreshold,
      exTblUnsThreshold = exTblUnsThreshold,
      chnlSettings = chnlSettings
    )

    if (!is.null(densAdaptive)) {
      chnl <- .getCpUnsLocGetChnl(exTblStimThreshold)
      stageChnl <- file.path(stage, chnl)
      .intSaveNm(
        "bwCpUnsLocAdaptive",
        densAdaptive$bw,
        .getInd(exTblStimThreshold),
        stageChnl,
        pathProject
      )
      return(densAdaptive)
    }
  }

  bw <- .getCpUnsLocGetDensRawDensitiesBw(
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold,
    bw = chnlSettings$bw,
    bwMin = chnlSettings$bwMin,
    bwMax = chnlSettings$bwMax,
    bwFallback = chnlSettings$bwFallback,
    bwMtd = chnlSettings$bwMtd,
    bwAdj = chnlSettings$bwAdj,
    bwNcellMin = chnlSettings$bwNcellMin,
    bwNcellMax = chnlSettings$bwNcellMax,
    chnlSettings = chnlSettings
  )
  chnl <- .getCpUnsLocGetChnl(exTblStimThreshold)
  stageChnl <- file.path(stage, chnl)
  .intSaveNm(
    "bwCpUnsLoc",
    bw,
    .getInd(exTblStimThreshold),
    stageChnl,
    pathProject
  )
  densStim <- .getCpUnsLocGetDensRawDensitiesStim(
    exTblStimThreshold = exTblStimThreshold,
    bw = bw
  )
  densUns <- .getCpUnsLocGetDensRawDensitiesUns(
    exTblUnsThreshold = exTblUnsThreshold,
    densStim = densStim,
    bw = bw
  )
  list(stim = densStim, uns = densUns, bw = bw)
}


#' @keywords internal
.getCpUnsLocUseAdaptiveBw <- function(chnlSettings) {
  (isTRUE(chnlSettings$bwAdaptive %||% FALSE) ||
    .getCpUnsLocHasManualAdaptiveBw(chnlSettings)) &&
    is.null(chnlSettings$bw)
}

#' @keywords internal
.getCpUnsLocHasManualAdaptiveBw <- function(chnlSettings) {
  is.finite(suppressWarnings(as.numeric(chnlSettings$bwAdaptiveCore)[1])) ||
    is.finite(suppressWarnings(as.numeric(chnlSettings$bwAdaptiveExtra)[1])) ||
    is.finite(suppressWarnings(as.numeric(chnlSettings$bwAdaptiveCrossover)[1]))
}

#' @keywords internal
.getCpUnsLocManualAdaptiveBwOnGrid <- function(chnlSettings, grid) {
  bwCore <- suppressWarnings(as.numeric(chnlSettings$bwAdaptiveCore)[1])
  bwExtra <- suppressWarnings(as.numeric(chnlSettings$bwAdaptiveExtra)[1])
  crossover <- suppressWarnings(as.numeric(chnlSettings$bwAdaptiveCrossover)[1])
  transitionWidth <- suppressWarnings(
    as.numeric(chnlSettings$bwAdaptiveTransitionWidth %||% 0)[1]
  )

  if (
    !is.finite(bwCore) ||
      bwCore <= 0 ||
      !is.finite(bwExtra) ||
      bwExtra <= 0 ||
      !is.finite(crossover)
  ) {
    return(NULL)
  }

  bw <- .bwNormBwFromCrossover(
    bin = grid,
    bwCore = bwCore,
    bwExtra = bwExtra,
    crossover = crossover,
    transitionWidth = transitionWidth
  )

  bw <- .getCpUnsLocRepairBwGrid(bw)

  list(
    bin = grid,
    bw = bw,
    bwCore = bwCore,
    bwExtra = bwExtra,
    bwAdaptiveCoreManual = bwCore,
    bwAdaptiveExtraManual = bwExtra,
    bwAdaptiveCrossover = crossover,
    bwAdaptiveTransitionWidth = transitionWidth,
    adaptive = TRUE,
    manual = TRUE
  )
}

#' @keywords internal
.getCpUnsLocGetDensRawDensitiesAdaptive <- function(
  exTblStimThreshold,
  exTblUnsThreshold,
  chnlSettings
) {
  xStim <- .getCut(exTblStimThreshold)
  xUns <- .getCut(exTblUnsThreshold)
  xStim <- suppressWarnings(as.numeric(xStim))
  xUns <- suppressWarnings(as.numeric(xUns))
  xStim <- xStim[is.finite(xStim)]
  xUns <- xUns[is.finite(xUns)]

  if (
    length(xStim) < 20L ||
      length(unique(xStim)) < 5L ||
      length(xUns) < 20L ||
      length(unique(xUns)) < 5L
  ) {
    return(NULL)
  }

  densityN <- .bwAsSafeSampleN(
    chnlSettings$bwAdaptiveDensityN %||%
      chnlSettings$normDensityN %||%
      512L,
    default = 512L,
    lower = 64L
  )

  grid <- .getCpUnsLocAdaptiveGrid(
    x = c(xStim, xUns),
    n = densityN,
    padFrac = chnlSettings$bwAdaptivePadFrac %||% 0.15
  )

  if (length(grid) < 10L) {
    return(NULL)
  }

  manualBwObj <- .getCpUnsLocManualAdaptiveBwOnGrid(
    chnlSettings = chnlSettings,
    grid = grid
  )

  if (!is.null(manualBwObj)) {
    bwStimObj <- manualBwObj
    bwUnsObj <- manualBwObj
    bwStimGrid <- manualBwObj$bw
    bwUnsGrid <- manualBwObj$bw
  } else {
    bwStimObj <- .getCpUnsLocGetDensRawDensitiesBwAdaptiveOne(
      x = xStim,
      chnlSettings = chnlSettings
    )
    bwUnsObj <- .getCpUnsLocGetDensRawDensitiesBwAdaptiveOne(
      x = xUns,
      chnlSettings = chnlSettings
    )

    bwStimGrid <- .getCpUnsLocAdaptiveBwOnGrid(
      bwObj = bwStimObj,
      grid = grid,
      fallback = chnlSettings$bwFallback
    )
    bwUnsGrid <- .getCpUnsLocAdaptiveBwOnGrid(
      bwObj = bwUnsObj,
      grid = grid,
      fallback = chnlSettings$bwFallback
    )
  }

  if (
    length(bwStimGrid) != length(grid) ||
      length(bwUnsGrid) != length(grid) ||
      all(!is.finite(bwStimGrid)) ||
      all(!is.finite(bwUnsGrid))
  ) {
    return(NULL)
  }

  bwStimGrid <- .getCpUnsLocRepairBwGrid(bwStimGrid)
  bwUnsGrid <- .getCpUnsLocRepairBwGrid(bwUnsGrid)

  # Preliminary densities use each sample's own adaptive bandwidth curve. These
  # are only weighting curves for constructing the shared bandwidth curve.
  densStimPre <- .getCpUnsLocDensityAdaptiveGrid(
    x = xStim,
    grid = grid,
    bwGrid = bwStimGrid,
    normalise = TRUE,
    probGMin = attr(exTblStimThreshold, "probGMin")
  )

  densUnsPre <- .getCpUnsLocDensityAdaptiveGrid(
    x = xUns,
    grid = grid,
    bwGrid = bwUnsGrid,
    normalise = TRUE,
    probGMin = attr(exTblUnsThreshold, "probGMin")
  )

  if (is.null(densStimPre) || is.null(densUnsPre)) {
    return(NULL)
  }

  denom <- densStimPre$y + densUnsPre$y
  bwShared <- ifelse(
    is.finite(denom) & denom > 0,
    (densStimPre$y * bwStimGrid + densUnsPre$y * bwUnsGrid) / denom,
    rowMeans(cbind(bwStimGrid, bwUnsGrid), na.rm = TRUE)
  )
  bwShared <- .getCpUnsLocRepairBwGrid(bwShared)

  # Final densities use the same grid and the same location-specific bandwidth
  # vector, then are normalised separately over the full padded grid.
  densStim <- .getCpUnsLocDensityAdaptiveGrid(
    x = xStim,
    grid = grid,
    bwGrid = bwShared,
    normalise = TRUE,
    probGMin = attr(exTblStimThreshold, "probGMin")
  )

  densUns <- .getCpUnsLocDensityAdaptiveGrid(
    x = xUns,
    grid = grid,
    bwGrid = bwShared,
    normalise = TRUE,
    probGMin = attr(exTblUnsThreshold, "probGMin")
  )

  if (is.null(densStim) || is.null(densUns)) {
    return(NULL)
  }

  list(
    stim = densStim,
    uns = densUns,
    bw = list(
      adaptive = TRUE,
      grid = grid,
      stim = bwStimObj,
      uns = bwUnsObj,
      stimGrid = bwStimGrid,
      unsGrid = bwUnsGrid,
      sharedGrid = bwShared,
      densStimWeight = densStimPre$y,
      densUnsWeight = densUnsPre$y
    )
  )
}

#' @keywords internal
.getCpUnsLocGetDensRawDensitiesBwAdaptiveOne <- function(
  x,
  chnlSettings
) {
  .bwCalcOne(
    x = x,
    bwMtd = chnlSettings$bwMtd %||% "hpi1Norm",
    bwAdj = chnlSettings$bwAdj %||% 1,
    bwNcellMin = chnlSettings$bwNcellMin,
    bwNcellMax = chnlSettings$bwNcellMax,
    normPeakFrac = chnlSettings$normPeakFrac %||% 0.1,
    normPeakMinRel = chnlSettings$normPeakMinRel %||% 0.75,
    normExtraFrac = chnlSettings$normExtraFrac %||% 0.2,
    normExtraMax = chnlSettings$normExtraMax %||% Inf,
    normExtraJitterFrac = chnlSettings$normExtraJitterFrac %||% 0.25,
    normLambda = chnlSettings$normLambda %||% seq(-2, 2, length.out = 81),
    normDensityN = chnlSettings$normDensityN %||% 512L,
    normExcessBwMtd = chnlSettings$normExcessBwMtd %||% "hpi3",
    normExcessNcell = chnlSettings$normExcessNcell %||% 10000L,
    normAdaptiveNcell = chnlSettings$normAdaptiveNcell %||%
      chnlSettings$bwAdaptiveNcell %||%
      2500L,
    bwAdaptiveCore = chnlSettings$bwAdaptiveCore,
    bwAdaptiveExtra = chnlSettings$bwAdaptiveExtra,
    bwAdaptiveCrossover = chnlSettings$bwAdaptiveCrossover,
    bwAdaptiveTransitionWidth = chnlSettings$bwAdaptiveTransitionWidth %||% 0,
    normMtd = chnlSettings$normMtd %||% "moments",
    adaptive = TRUE
  )
}

#' @keywords internal
.getCpUnsLocAdaptiveGrid <- function(
  x,
  n = 512L,
  padFrac = 0.15
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(numeric(0L))
  }

  n <- .bwAsSafeSampleN(n, default = 512L, lower = 16L)
  rangeVec <- range(x, na.rm = TRUE)
  rangeWidth <- diff(rangeVec)

  if (!is.finite(rangeWidth) || rangeWidth <= 0) {
    rangeWidth <- .bwRobustSd(x)
  }
  if (!is.finite(rangeWidth) || rangeWidth <= 0) {
    rangeWidth <- 1
  }

  padFrac <- suppressWarnings(as.numeric(padFrac)[1])
  if (!is.finite(padFrac) || padFrac < 0) {
    padFrac <- 0.15
  }

  pad <- padFrac * rangeWidth
  seq(
    from = rangeVec[[1]] - pad,
    to = rangeVec[[2]] + pad,
    length.out = n
  )
}

#' @keywords internal
.getCpUnsLocAdaptiveBwOnGrid <- function(
  bwObj,
  grid,
  fallback = NULL
) {
  grid <- suppressWarnings(as.numeric(grid))
  grid <- grid[is.finite(grid)]

  if (
    is.list(bwObj) &&
      all(c("bin", "bw") %in% names(bwObj)) &&
      length(bwObj$bin) >= 2L &&
      length(bwObj$bin) == length(bwObj$bw)
  ) {
    bwGrid <- stats::approx(
      x = suppressWarnings(as.numeric(bwObj$bin)),
      y = suppressWarnings(as.numeric(bwObj$bw)),
      xout = grid,
      rule = 2
    )$y
    return(.getCpUnsLocRepairBwGrid(bwGrid))
  }

  bwScalar <- suppressWarnings(as.numeric(bwObj)[1])
  if (!is.finite(bwScalar) || bwScalar <= 0) {
    bwScalar <- suppressWarnings(as.numeric(fallback)[1])
  }
  if (!is.finite(bwScalar) || bwScalar <= 0) {
    bwScalar <- NA_real_
  }

  rep(bwScalar, length(grid))
}

#' @keywords internal
.getCpUnsLocRepairBwGrid <- function(bwGrid) {
  bwGrid <- suppressWarnings(as.numeric(bwGrid))

  if (length(bwGrid) == 0L) {
    return(bwGrid)
  }

  good <- is.finite(bwGrid) & bwGrid > 0
  if (!any(good)) {
    return(rep(NA_real_, length(bwGrid)))
  }

  bwGrid[!good] <- stats::median(bwGrid[good], na.rm = TRUE)
  pmax(bwGrid, .Machine$double.eps)
}

#' @keywords internal
.getCpUnsLocDensityAdaptiveGrid <- function(
  x,
  grid,
  bwGrid,
  normalise = TRUE,
  probGMin = NULL
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  grid <- suppressWarnings(as.numeric(grid))
  bwGrid <- suppressWarnings(as.numeric(bwGrid))

  ok <- is.finite(grid) & is.finite(bwGrid) & bwGrid > 0
  grid <- grid[ok]
  bwGrid <- bwGrid[ok]

  if (
    length(x) < 2L ||
      length(unique(x)) < 2L ||
      length(grid) < 2L ||
      length(grid) != length(bwGrid)
  ) {
    return(NULL)
  }

  y <- vapply(
    seq_along(grid),
    function(i) {
      mean(stats::dnorm(
        x = grid[[i]],
        mean = x,
        sd = bwGrid[[i]]
      ))
    },
    numeric(1)
  )

  y <- pmax(y, 0)

  y <- .getCpUnsLocDensityNormalizeAndScale(
    grid = grid,
    y = y,
    normalise = normalise,
    probGMin = probGMin
  )

  out <- list(
    x = grid,
    y = y,
    bw = bwGrid,
    n = length(x),
    call = match.call(),
    data.name = deparse(substitute(x)),
    has.na = FALSE
  )
  class(out) <- "density"
  out
}

#' @keywords internal
.getCpUnsLocTrapz <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2L || length(x) != length(y)) {
    return(NA_real_)
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

#' @keywords internal
.getCpUnsLocDensityNormalizeAndScale <- function(
  grid,
  y,
  normalise = TRUE,
  probGMin = NULL
) {
  grid <- suppressWarnings(as.numeric(grid))
  y <- suppressWarnings(as.numeric(y))

  if (length(grid) != length(y)) {
    return(y)
  }

  y <- pmax(y, 0)

  if (isTRUE(normalise)) {
    area <- .getCpUnsLocTrapz(grid, y)
    if (is.finite(area) && area > 0) {
      y <- y / area
    }
  }

  probGMin <- suppressWarnings(as.numeric(probGMin)[1])
  if (is.finite(probGMin)) {
    probGMin <- max(0, min(1, probGMin))
    y <- y * probGMin
  }

  y
}


#' @keywords internal
.getCpUnsLocGetDensRawDensitiesBw <- function(
  exTblStimThreshold,
  exTblUnsThreshold,
  bw,
  bwMin,
  bwMax,
  bwFallback,
  bwMtd,
  bwAdj,
  bwNcellMin,
  bwNcellMax,
  chnlSettings = NULL
) {
  if (!is.null(bw)) {
    return(bw)
  }
  bwStim <- .getCpUnsLocGetDensRawDensitiesBwInit(
    .data = .getCut(exTblStimThreshold),
    bwMin = bwMin,
    bwMax = bwMax,
    bwFallback = bwFallback,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax,
    chnlSettings = chnlSettings
  )
  bwUns <- .getCpUnsLocGetDensRawDensitiesBwInit(
    .data = .getCut(exTblUnsThreshold),
    bwMin = bwMin,
    bwMax = bwMax,
    bwFallback = bwFallback,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax,
    chnlSettings = chnlSettings
  )
  min(bwUns, bwStim)
}


#' @keywords internal
.getCpUnsLocGetDensRawDensitiesStim <- function(
  exTblStimThreshold,
  bw
) {
  densObj <- density(.getCut(exTblStimThreshold), bw = bw)
  if (is.null(attr(exTblStimThreshold, "probGMin"))) {
    return(densObj)
  }
  densObj$y <- densObj$y * attr(exTblStimThreshold, "probGMin")
  densObj
}

#' @keywords internal
.getCpUnsLocGetDensRawDensitiesUns <- function(
  exTblUnsThreshold,
  densStim,
  bw
) {
  densObj <- density(
    .getCut(exTblUnsThreshold),
    from = min(densStim$x),
    to = max(densStim$x),
    bw = bw
  )
  if (is.null(attr(exTblUnsThreshold, "probGMin"))) {
    return(densObj)
  }
  densObj$y <- densObj$y * attr(exTblUnsThreshold, "probGMin")
  densObj
}

#' @keywords internal
.getCpUnsLocGetDensRawTabulate <- function(
  stimX,
  stimY,
  unsX,
  unsY
) {
  densTblRawStim <- tibble::tibble(xStim = stimX, yStim = stimY)
  densTblRawWide <- .getCpUnsLocGetDensRawTabulateUnsInterp(
    .data = densTblRawStim,
    unsX = unsX,
    unsY = unsY
  )
  .getCpUnsLocGetDensRawTabulateFormat(densTblRawWide)
}

#' @keywords internal
.getCpUnsLocGetDensRawTabulateUnsInterp <- function(
  .data,
  unsX,
  unsY
) {
  yUns <- stats::approx(
    x = suppressWarnings(as.numeric(unsX)),
    y = suppressWarnings(as.numeric(unsY)),
    xout = suppressWarnings(as.numeric(.data$xStim)),
    rule = 2
  )$y

  .data |>
    dplyr::mutate(yUns = .env$yUns)
}


#' @keywords internal
.getCpUnsLocGetDensRawTabulateFormat <- function(.data) {
  .data |>
    tidyr::pivot_longer(yStim:yUns, names_to = "stim", values_to = "dens") |> # nolint
    dplyr::mutate(stim = ifelse(stim == "yStim", "yes", "no")) # nolint
}

# get probabilities
# -------------------
#' @keywords internal
.getCpUnsLocGetProbTbl <- function(
  densTblRaw,
  stage,
  cpMin,
  exVecStimThreshold,
  exVecUnsThreshold
) {
  .debug("Normalising probabilities") # nolint

  probTbl <- .getCpUnsLocGetProbTblInit(densTblRaw, cpMin)

  maxXFilter <- if (length(exVecStimThreshold) < 20L) {
    Inf
  } else {
    quantile(exVecStimThreshold, (length(exVecStimThreshold) - 20L) / length(exVecStimThreshold))
  }

  probTblPos <- .getCpUnsLocProbTblFilter(
    exVecStimThreshold = exVecStimThreshold,
    exVecUnsThreshold = exVecUnsThreshold,
    probTbl = probTbl,
    maxXFilter = maxXFilter,
    stage = stage
  )
  list(all = probTbl, pos = probTblPos)
}

#' @keywords internal
.getCpUnsLocGetProbTblInit <- function(densTblRaw, cpMin) {
  densTblRaw |>
    tidyr::pivot_wider(
      id_cols = xStim,
      names_from = stim,
      values_from = dens # nolint
    ) |>
    dplyr::mutate(
      probStim = 1 - no / yes, # nolint
      probStim = ifelse(yes == 0 & no == 0, 0, probStim), # nolint
      probStimNorm = pmin(1, probStim), # nolint
      probStimNorm = pmax(0, probStimNorm) # nolint
    ) |>
    dplyr::filter(xStim > cpMin)
}

#' @keywords internal
.getCpUnsLocGetDensRawDensitiesBwInit <- function(
  .data,
  bwMin,
  bwMax,
  bwFallback,
  bwMtd,
  bwAdj,
  bwNcellMin,
  bwNcellMax,
  chnlSettings = NULL
) {
  chnlSettings <- chnlSettings %||% list()
  .data <- suppressWarnings(as.numeric(.data))
  .data <- .data[is.finite(.data)]

  if (length(.data) < 2L || length(unique(.data)) < 2L) {
    return(bwFallback)
  }

  bwCalc <- .bwCalcOne(
    x = .data,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax,
    normPeakFrac = chnlSettings$normPeakFrac %||% 0.1,
    normPeakMinRel = chnlSettings$normPeakMinRel %||% 0.75,
    normExtraFrac = chnlSettings$normExtraFrac %||% 0.2,
    normExtraMax = chnlSettings$normExtraMax %||% Inf,
    normExtraJitterFrac = chnlSettings$normExtraJitterFrac %||% 0.25,
    normLambda = chnlSettings$normLambda %||% seq(-2, 2, length.out = 81),
    normDensityN = chnlSettings$normDensityN %||% 512L,
    normExcessBwMtd = chnlSettings$normExcessBwMtd %||% "hpi3",
    normExcessNcell = chnlSettings$normExcessNcell %||% 10000L,
    normAdaptiveNcell = chnlSettings$normAdaptiveNcell %||%
      chnlSettings$bwAdaptiveNcell %||%
      2500L,
    bwAdaptiveCore = chnlSettings$bwAdaptiveCore,
    bwAdaptiveExtra = chnlSettings$bwAdaptiveExtra,
    bwAdaptiveCrossover = chnlSettings$bwAdaptiveCrossover,
    bwAdaptiveTransitionWidth = chnlSettings$bwAdaptiveTransitionWidth %||% 0,
    normMtd = chnlSettings$normMtd %||% "moments",
    adaptive = FALSE
  )

  if (!is.finite(bwCalc) || bwCalc <= 0) {
    return(bwFallback)
  }

  max(bwMin, min(as.numeric(bwCalc)[1], bwMax))
}


#' @keywords internal
.getCpUnsLocProbTblFilter <- function(
  exVecStimThreshold,
  exVecUnsThreshold,
  probTbl,
  maxXFilter,
  stage
) {
  .debug("Filtering before smoothing") # nolint

  densityExcMinStim <- density(exVecStimThreshold)
  densTblStim <- tibble::tibble(
    x = densityExcMinStim$x,
    y = densityExcMinStim$y
  )
  peakStim <- densTblStim |>
    dplyr::filter(y == max(y)) |> # nolint
    dplyr::pull("x") # nolint
  densityExcMinUns <- density(exVecUnsThreshold)
  densTblUns <- tibble::tibble(
    x = densityExcMinUns$x,
    y = densityExcMinUns$y
  )
  peakUns <- densTblUns |>
    dplyr::filter(y == max(y)) |> # nolint
    dplyr::pull("x") # nolint
  peak <- max(peakStim, peakUns)

  windowWidth <- 0.15 * diff(quantile(probTbl$xStim, c(0.05, 0.5)))

  probTbl <- probTbl |>
    dplyr::filter(xStim > peak + windowWidth) # nolint

  if (nrow(probTbl) < 20L) {
    return(probTbl)
  }

  probTbl <- probTbl |> dplyr::arrange(xStim)

  probFilterIdx <- which(probTbl$xStim < maxXFilter)
  probIncIdx <- setdiff(seq_len(nrow(probTbl)), probFilterIdx)
  probVec <- probTbl$probStimNorm[probFilterIdx] |>
    stats::setNames(probTbl$xStim)
  probVecGe0025FilterIdx <- which(cumsum(probVec) >= 0.025)
  probVecGe0075 <- probVec >= 0.075

  probTbl <- probTbl |>
    dplyr::mutate(
      minorResponseInd = probStimNorm >= 0.025,
      moderateResponseInd = probStimNorm >= 0.075,
      nRemaining = dplyr::n() - seq_len(dplyr::n()) + 1
    )
  probTbl |>
    dplyr::filter(
      cumsum(minorResponseInd) > 0 # nolint
    ) |>
    dplyr::mutate(
      probLargerCount = purrr::map_int(xStim, function(x) {
        sum(probTbl$moderateResponseInd[probTbl$xStim >= x])
      }),
      probLargerProp = probLargerCount / nRemaining # nolint
    ) |>
    dplyr::filter(
      probLargerProp > 0.25 # nolint
    ) |>
    dplyr::select(
      -c(
        probLargerProp,
        minorResponseInd,
        moderateResponseInd,
        nRemaining,
        probLargerCount
      )
    )
}

#' @keywords internal
.getCpUnsLocGetMinProbX <- function(probTblPos) {
  min(probTblPos$xStim)
}

#' @keywords internal
.getCpUnsLocCheckResponse <- function(probTblPos, exTblStimOrig) {
  nrow(probTblPos) == 0 ||
    max(.getCut(exTblStimOrig)) < .getCpUnsLocGetMinProbX(probTblPos)
}

#' @keywords internal
.getCpUnsLocGetDataMod <- function(
  exTblStimThreshold,
  exTblStimNoMin,
  exTblUnsThreshold,
  exTblUnsBias,
  probTblList,
  cpMin,
  stage
) {
  if (.getCpUnsLocCheckResponse(probTblList$pos, exTblStimNoMin)) {
    return(.getCpUnsLocConditionCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "No responding cells" # nolint
    ))
  }
  margin <- .getCpUnsLocGetDataModMargin(
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsNoMin = exTblUnsThreshold
  )
  binVec <- density(
    c(.getCut(exTblStimThreshold), .getCut(exTblUnsThreshold))
  )$x

  dataMod <- exTblStimThreshold
  dataMod <- dataMod[
    .getCut(dataMod) >=
      (min(.getCpUnsLocGetMinProbX(probTblList$pos) - margin)),
  ]
  if (nrow(dataMod) == 0L) {
    return(.getCpUnsLocConditionCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "No responding cells" # nolint
    ))
  }
  probVec <- try(
    approx(
      x = probTblList$pos$xStim,
      y = probTblList$pos$probStimNorm,
      xout = dataMod[[1]],
      method = "linear",
      f = 0.5,
      rule = 2
    )$y,
    silent = TRUE
  )
  if (inherits(probVec, "try-error")) {
    return(.getCpUnsLocConditionCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "No responding cells" # nolint
    ))
  }
  attr(dataMod, "binVec") <- binVec

  dataMod <- dataMod |>
    dplyr::mutate(probSmooth = probVec)

  .thinDataMod(dataMod, maxCellsPerBin = 20)
}

#' @keywords internal
.getCpUnsLocGetDataModMargin <- function(
  exTblStimNoMin,
  exTblUnsNoMin
) {
  abs(max(
    diff(.getCut(exTblStimNoMin)),
    diff(.getCut(exTblUnsNoMin))
  )) *
    0.05
}

# smooth
# ---------------------
#' @keywords internal
.getCpUnsLocGetProbSmooth <- function(
  dataMod,
  stage,
  pathProject,
  chnl,
  chnlSettings = list()
) {
  stageChnl <- file.path(stage, chnl)

  if (!.getCpUnsLocGetProbSmoothCheckNCell(dataMod)) {
    .intSaveNm(
      "not_enough_cells_to_smooth",
      NULL,
      .getInd(dataMod),
      stageChnl,
      pathProject
    )
    dataModOut <- .getCpUnsLocGetProbSmoothCheckNCellOut(dataMod)
    .intSaveNm(
      "probSmoothOut",
      dataModOut,
      .getInd(dataMod),
      stageChnl,
      pathProject
    )
    return(dataModOut)
  }

  smoothObj <- .getCpUnsLocGetProbSmoothActual(
    dataMod = dataMod,
    stage = stage,
    chnlSettings = chnlSettings
  )
  predVec <- .getCpUnsLocGetProbSmoothObjPred(smoothObj)
  dataModOut <- dataMod |> dplyr::mutate(pred = predVec)
  dataModOut <- .getCpUnsLocGetProbSmoothAttachDeriv(
    dataMod = dataModOut,
    smoothObj = smoothObj
  )
  .intSaveNm(
    "probSmoothOut",
    dataModOut,
    .getInd(dataMod),
    stageChnl,
    pathProject
  )
  dataModOut
}

#' @keywords internal
.getCpUnsLocGetProbSmoothCheckNCell <- function(dataMod) {
  is.data.frame(dataMod) && nrow(dataMod) > 10
}

#' @keywords internal
.getCpUnsLocGetProbSmoothCheckNCellOut <- function(dataMod) {
  if (is.data.frame(dataMod)) {
    dataMod |> dplyr::mutate(pred = probSmooth - 1e-4) # nolint
  } else {
    dataMod
  }
}

#' @keywords internal
.getCpUnsLocGetProbSmoothObjPred <- function(smoothObj) {
  if (is.list(smoothObj) && "pred" %in% names(smoothObj)) {
    return(smoothObj$pred)
  }
  smoothObj
}

#' @keywords internal
.getCpUnsLocGetProbSmoothAttachDeriv <- function(dataMod, smoothObj) {
  if (!is.list(smoothObj)) {
    return(dataMod)
  }
  if (!is.null(smoothObj$derivTbl)) {
    attr(dataMod, "locProbDerivTbl") <- smoothObj$derivTbl
  }
  if (!is.null(smoothObj$method)) {
    attr(dataMod, "locProbSmoothMethod") <- smoothObj$method
  }
  dataMod
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActual <- function(
  dataMod,
  stage,
  chnlSettings = list()
) {
  fit1 <- .getCpUnsLocGetProbSmoothActualFirst(dataMod, stage)
  .getCpUnsLocGetProbSmoothActualFirstResponse(
    fit = fit1,
    dataMod = dataMod,
    stage = stage,
    chnlSettings = chnlSettings
  )
}

.getCpUnsLocGetProbSmoothActualFirst <- function(dataMod, stage) {
  .debug("Smoothing I")

  idxMod <- attr(dataMod, "idxMod") %||% seq_len(nrow(dataMod))
  chnl <- .getCpUnsLocGetChnl(dataMod)

  dataMod <- dataMod[idxMod, ] |>
    dplyr::mutate(
      probSmooth = pmin(probSmooth, 0.999),
      probSmooth = pmax(probSmooth, 0.001)
    )

  n <- nrow(dataMod)

  if (n <= 4L) {
    return(try(stop(), silent = TRUE)) # nolint
  }

  k <- min(n - 1L, 20L)

  fml <- as.formula(paste0(
    "probSmooth ~ s(`",
    chnl,
    "`, bs = 'mpi', k = ",
    k,
    ", m = c(2, 1))"
  ))

  try(
    scam::scam(
      fml,
      family = "quasibinomial",
      data = dataMod,
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        devtol.fit = 0.5,
        steptol.fit = 1e-1,
        maxHalf = 5,
        bfgs = list(steptol.bfgs = 1e-1),
        maxit = 1e1
      )
    ),
    silent = TRUE
  )
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualFirstResponse <- function(
  fit,
  dataMod,
  stage,
  chnlSettings = list()
) {
  if (.getCpUnsLocGetProbSmoothActualCheck(fit, dataMod)) {
    .debug("Smoothed") # nolint
    return(
      .getCpUnsLocGetProbSmoothActualResponseSuccess(
        fit = fit,
        dataMod = dataMod,
        chnlSettings = chnlSettings,
        method = "scam_mpi"
      )
    )
  }
  .getCpUnsLocGetProbSmoothActualFirstResponseFailure(
    stage = stage,
    dataMod = dataMod,
    chnlSettings = chnlSettings
  )
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualCheck <- function(fit, dataMod) {
  if (inherits(fit, "try-error") || is.null(fit)) {
    return(FALSE)
  }
  outList <- .getCpUnsLocGetProbSmoothActualResponseSuccess(
    fit = fit,
    dataMod = dataMod
  )
  !(all(outList$pred > 0.99) || outList$meanAbsError > 0.3) # nolint
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualResponseSuccess <- function(
  fit,
  dataMod,
  chnlSettings = list(),
  method = NA_character_
) {
  # predict on full dataset
  predVec <- predict(fit, newdata = dataMod, type = "response")
  meanAbsError <- mean(abs(predVec - dataMod$probSmooth))
  derivTbl <- .getCpUnsLocGetProbSmoothDerivativeTbl(
    fit = fit,
    dataMod = dataMod,
    chnlSettings = chnlSettings
  )
  list(
    "pred" = predVec,
    "meanAbsError" = meanAbsError,
    "derivTbl" = derivTbl,
    "method" = method
  )
}

#' @keywords internal
.getCpUnsLocGetProbSmoothDerivativeTbl <- function(
  fit,
  dataMod,
  chnlSettings = list()
) {
  chnl <- .getCpUnsLocGetChnl(dataMod)
  x <- suppressWarnings(as.numeric(.getCut(dataMod)))
  x <- x[is.finite(x)]
  if (length(x) < 4L || diff(range(x)) <= 0) {
    return(NULL)
  }

  nGrid <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locFlatDerivGridN",
    512L
  )
  nGrid <- suppressWarnings(as.integer(nGrid[1]))
  if (!is.finite(nGrid) || nGrid < 25L) {
    nGrid <- 512L
  }

  epsFrac <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locFlatDerivEpsFrac",
    1e-5
  )
  epsFrac <- suppressWarnings(as.numeric(epsFrac[1]))
  if (!is.finite(epsFrac) || epsFrac <= 0) {
    epsFrac <- 1e-5
  }

  xRange <- range(x, na.rm = TRUE)
  xWidth <- diff(xRange)
  xGrid <- seq(xRange[1], xRange[2], length.out = nGrid)
  eps <- max(
    xWidth * epsFrac,
    sqrt(.Machine$double.eps) * max(abs(xRange), 1, na.rm = TRUE)
  )
  if (!is.finite(eps) || eps <= 0) {
    return(NULL)
  }

  xLeft <- pmax(xRange[1], xGrid - eps)
  xRight <- pmin(xRange[2], xGrid + eps)
  denom <- xRight - xLeft
  if (any(!is.finite(denom) | denom <= 0)) {
    return(NULL)
  }

  newData <- dataMod[rep(1L, length(xGrid)), , drop = FALSE]
  newData[[chnl]] <- xGrid
  predGrid <- try(
    predict(fit, newdata = newData, type = "response"),
    silent = TRUE
  )
  newData[[chnl]] <- xLeft
  predLeft <- try(
    predict(fit, newdata = newData, type = "response"),
    silent = TRUE
  )
  newData[[chnl]] <- xRight
  predRight <- try(
    predict(fit, newdata = newData, type = "response"),
    silent = TRUE
  )
  if (
    inherits(predGrid, "try-error") ||
      inherits(predLeft, "try-error") ||
      inherits(predRight, "try-error")
  ) {
    return(NULL)
  }

  deriv <- (as.numeric(predRight) - as.numeric(predLeft)) / denom
  deriv[!is.finite(deriv)] <- 0
  deriv <- pmax(0, deriv)
  tibble::tibble(
    x = xGrid,
    pred = as.numeric(predGrid),
    deriv = deriv
  )
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualFirstResponseFailure <- function(
  stage,
  dataMod,
  chnlSettings = list()
) {
  fit2 <- .getCpUnsLocGetProbSmoothActualSecond(dataMod, stage)
  if (.getCpUnsLocGetProbSmoothActualCheck(fit2, dataMod)) {
    .debug("Smoothed") # nolint
    return(
      .getCpUnsLocGetProbSmoothActualResponseSuccess(
        fit = fit2,
        dataMod = dataMod,
        chnlSettings = chnlSettings,
        method = "scam_micv"
      )
    )
  }
  .getCpUnsLocGetProbSmoothActualThird(dataMod, stage)
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualSecond <- function(dataMod, stage) {
  .debug("Smoothing II") # nolint

  idxMod <- attr(dataMod, "idxMod") %||% seq_len(nrow(dataMod))
  chnl <- .getCpUnsLocGetChnl(dataMod)

  dataMod <- dataMod[idxMod, ] |>
    dplyr::mutate(
      probSmooth = pmin(probSmooth, 0.999),
      probSmooth = pmax(probSmooth, 0.001)
    )

  n <- nrow(dataMod)

  if (n <= 4L) {
    return(NULL)
  }

  k <- min(n - 1L, 20L)

  fml <- as.formula(paste0(
    "probSmooth ~ s(`",
    chnl,
    "`, bs = 'micv', k = ",
    k,
    ", m = c(2, 1))"
  ))

  try(
    suppressWarnings(scam::scam(
      fml,
      family = "binomial",
      data = dataMod,
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        devtol.fit = 0.01
      )
    )),
    silent = TRUE
  )
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualThird <- function(dataMod, stage) {
  .debug("Failed to smooth") # nolint
  list(
    "pred" = dataMod$probSmooth - 0.0001,
    "meanAbsError" = NA_real_,
    "derivTbl" = NULL,
    "method" = "probSmooth_fallback"
  )
}


#' @keywords internal
.getCpUnsLocThresholdOrigin <- function(
  locGenerated,
  locGeneratedDirect,
  locSource,
  locReason = NA_character_
) {
  locGenerated <- locGenerated %in% TRUE
  locGeneratedDirect <- locGeneratedDirect %in% TRUE
  locSource <- as.character(locSource %||% NA_character_)
  locReason <- as.character(locReason %||% NA_character_)

  if (locGenerated && locGeneratedDirect && locSource %in% "direct") {
    return("condition_detected_response")
  }
  if (locGenerated && locSource %in% "combined") {
    return("sample_imputed_from_other_stim_conditions")
  }
  if (locGenerated && locSource %in% "prejoin") {
    return("prejoin_generated_from_joined_stim_conditions")
  }
  if (locGenerated && locSource %in% "unstim_summary") {
    return("unstim_summary_from_generated_thresholds")
  }
  if (!locGenerated) {
    return("not_generated_fallback")
  }
  if (!is.na(locReason) && nzchar(locReason)) {
    return(paste0("generated_", locReason))
  }
  "generated_unknown_source"
}

#' @keywords internal
.getCpUnsLocPropBsAtCp <- function(cp, exTblStim, exTblUns) {
  cp <- suppressWarnings(as.numeric(cp))[1]
  if (
    !is.finite(cp) ||
      is.null(exTblStim) ||
      !is.data.frame(exTblStim) ||
      is.null(exTblUns) ||
      !is.data.frame(exTblUns) ||
      nrow(exTblStim) == 0L ||
      nrow(exTblUns) == 0L
  ) {
    return(tibble::tibble(
      nCellStim = if (is.data.frame(exTblStim)) {
        nrow(exTblStim)
      } else {
        NA_integer_
      },
      nCellUns = if (is.data.frame(exTblUns)) nrow(exTblUns) else NA_integer_,
      propStim = NA_real_,
      propUns = NA_real_,
      propBs = NA_real_
    ))
  }
  propStim <- sum(.getCut(exTblStim) >= cp, na.rm = TRUE) / nrow(exTblStim)
  propUns <- sum(.getCut(exTblUns) >= cp, na.rm = TRUE) / nrow(exTblUns)
  tibble::tibble(
    nCellStim = nrow(exTblStim),
    nCellUns = nrow(exTblUns),
    propStim = propStim,
    propUns = propUns,
    propBs = propStim - propUns
  )
}

#' @keywords internal
.getCpUnsLocSelectedThresholdRow <- function(dataThreshold, cp) {
  if (!is.data.frame(dataThreshold) || nrow(dataThreshold) == 0L) {
    return(NULL)
  }
  cp <- suppressWarnings(as.numeric(cp))[1]
  if (!is.finite(cp)) {
    return(NULL)
  }
  cutVec <- suppressWarnings(as.numeric(.getCut(dataThreshold)))
  idx <- which(is.finite(cutVec))
  if (length(idx) == 0L) {
    return(NULL)
  }
  idx <- idx[which.min(abs(cutVec[idx] - cp))]
  dataThreshold[idx, , drop = FALSE]
}

#' @keywords internal
.getCpUnsLocConditionDetailRow <- function(
  cpObj,
  dataThreshold,
  exTblStimOrig,
  exTblUnsOrig,
  exTblStimNoMin,
  bias,
  stage,
  chnl
) {
  cp <- suppressWarnings(as.numeric(cpObj$cp))[1]
  selectedRow <- .getCpUnsLocSelectedThresholdRow(dataThreshold, cp)
  freqTbl <- .getCpUnsLocPropBsAtCp(
    cp = cp,
    exTblStim = exTblStimOrig,
    exTblUns = exTblUnsOrig
  )

  propBsEst <- NA_real_
  propBsDiff <- NA_real_
  if (!is.null(selectedRow) && "propBsDiff" %in% names(selectedRow)) {
    propBsDiff <- suppressWarnings(as.numeric(selectedRow$propBsDiff[1]))
    if (is.finite(propBsDiff) && "propBs" %in% names(selectedRow)) {
      propBsEst <- suppressWarnings(as.numeric(selectedRow$propBs[1])) -
        propBsDiff
    }
  }

  tibble::tibble(
    detailLevel = "condition",
    stage = stage,
    chnl = chnl,
    ind = as.character(.getInd(exTblStimNoMin)),
    threshold = cp,
    thresholdOrigin = .getCpUnsLocThresholdOrigin(
      locGenerated = cpObj$locGenerated,
      locGeneratedDirect = cpObj$locGeneratedDirect,
      locSource = cpObj$locSource,
      locReason = cpObj$locReason
    ),
    locGenerated = cpObj$locGenerated %in% TRUE,
    locGeneratedDirect = cpObj$locGeneratedDirect %in% TRUE,
    locSource = as.character(cpObj$locSource %||% NA_character_),
    locReason = as.character(cpObj$locReason %||% NA_character_),
    bias = suppressWarnings(as.numeric(bias))[1],
    propBsEst = propBsEst,
    propBsDiff = propBsDiff
  ) |>
    dplyr::bind_cols(freqTbl)
}

#' @keywords internal
.getCpUnsLocExByInd <- function(exList, ind, pos = NULL) {
  ind <- as.character(ind)
  if (!is.null(names(exList)) && ind %in% names(exList)) {
    return(exList[[ind]])
  }
  if (!is.null(pos) && is.finite(pos) && pos >= 1L && pos <= length(exList)) {
    return(exList[[pos]])
  }
  NULL
}

#' @keywords internal
.getCpUnsLocSampleDetailTbl <- function(
  cpVec,
  exListOrig,
  indUns,
  indStim,
  stage,
  chnl
) {
  meta <- .getCpUnsLocMetaFromCp(cpVec)
  exTblUns <- .getCpUnsLocExByInd(exListOrig, indUns, pos = 1L)

  purrr::map_df(seq_along(indStim), function(i) {
    indCurr <- as.character(indStim[[i]])
    cp <- suppressWarnings(as.numeric(cpVec[indCurr]))[1]
    metaRow <- meta[meta$ind == indCurr, , drop = FALSE]
    if (nrow(metaRow) == 0L) {
      metaRow <- tibble::tibble(
        ind = indCurr,
        locGenerated = FALSE,
        locGeneratedDirect = FALSE,
        locSource = "not_calculated",
        locReason = NA_character_
      )
    }
    exTblStim <- .getCpUnsLocExByInd(exListOrig, indCurr, pos = i + 1L)
    freqTbl <- .getCpUnsLocPropBsAtCp(
      cp = cp,
      exTblStim = exTblStim,
      exTblUns = exTblUns
    )
    tibble::tibble(
      detailLevel = "sample",
      stage = stage,
      chnl = chnl,
      ind = indCurr,
      threshold = cp,
      thresholdOrigin = .getCpUnsLocThresholdOrigin(
        locGenerated = metaRow$locGenerated[1],
        locGeneratedDirect = metaRow$locGeneratedDirect[1],
        locSource = metaRow$locSource[1],
        locReason = metaRow$locReason[1]
      ),
      locGenerated = metaRow$locGenerated[1] %in% TRUE,
      locGeneratedDirect = metaRow$locGeneratedDirect[1] %in% TRUE,
      locSource = as.character(metaRow$locSource[1] %||% NA_character_),
      locReason = as.character(metaRow$locReason[1] %||% NA_character_),
      propBsEst = NA_real_,
      propBsDiff = NA_real_
    ) |>
      dplyr::bind_cols(freqTbl)
  })
}

# get cp
#' @keywords internal
.getCpUnsLocGetCp <- function(
  dataMod,
  exTblStimOrig,
  exTblStimNoMin,
  exTblUnsOrig,
  exTblUnsBias,
  bias,
  cpMin,
  stage,
  pathProject,
  chnlSettings = list()
) {
  ind <- .getInd(exTblStimNoMin)
  chnl <- .getCpUnsLocGetChnl(exTblStimNoMin)
  stageChnl <- file.path(stage, chnl)
  if (!is.data.frame(dataMod)) {
    .intSaveNm("noDataModDf", NULL, ind, stageChnl, pathProject)
    .intSaveNm("cpInd", dataMod, ind, stageChnl, pathProject)
    objOut <- .getCpUnsLocGetCpEnsureMeta(
      obj = dataMod,
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      reason = "data_mod_not_available"
    )
    locDetailCondition <- .getCpUnsLocConditionDetailRow(
      cpObj = objOut,
      dataThreshold = NULL,
      exTblStimOrig = exTblStimOrig,
      exTblUnsOrig = exTblUnsOrig,
      exTblStimNoMin = exTblStimNoMin,
      bias = bias,
      stage = stage,
      chnl = chnl
    )
    .intSaveNm(
      "locDetailCondition",
      locDetailCondition,
      ind,
      stageChnl,
      pathProject
    )
    return(objOut)
  }

  trimObj <- .getCpUnsLocGetCpTrimBeforeThreshold(
    dataMod = dataMod,
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsBias = exTblUnsBias,
    cpMin = cpMin,
    stage = stage,
    chnlSettings = chnlSettings
  )
  .intSaveNm("dataModTrimInfo", trimObj$info, ind, stageChnl, pathProject)
  .intSaveNm("dataModTrim", trimObj$dataMod, ind, stageChnl, pathProject)

  if (!is.null(trimObj$cp)) {
    cpInd <- trimObj$cp
    .intSave(ind, stageChnl, pathProject, cpInd)
    .debug("Completed loc gate for single sample") # nolint
    objOut <- .getCpUnsLocConditionOut(
      cp = cpInd,
      locGenerated = FALSE,
      locGeneratedDirect = FALSE,
      locSource = "not_calculated",
      locReason = trimObj$info$reason %||% "trim_returned_non_local_cutpoint"
    )
    locDetailCondition <- .getCpUnsLocConditionDetailRow(
      cpObj = objOut,
      dataThreshold = NULL,
      exTblStimOrig = exTblStimOrig,
      exTblUnsOrig = exTblUnsOrig,
      exTblStimNoMin = exTblStimNoMin,
      bias = bias,
      stage = stage,
      chnl = chnl
    )
    .intSaveNm(
      "locDetailCondition",
      locDetailCondition,
      ind,
      stageChnl,
      pathProject
    )
    return(objOut)
  }

  dataMod <- trimObj$dataMod
  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    cpInd <- .getCpUnsLocConditionCpNonLoc(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias
    )
    .intSave(ind, stageChnl, pathProject, cpInd)
    .debug("Completed loc gate for single sample") # nolint
    objOut <- .getCpUnsLocConditionOut(
      cp = cpInd,
      locGenerated = FALSE,
      locGeneratedDirect = FALSE,
      locSource = "not_calculated",
      locReason = "empty_data_mod_after_trimming"
    )
    locDetailCondition <- .getCpUnsLocConditionDetailRow(
      cpObj = objOut,
      dataThreshold = NULL,
      exTblStimOrig = exTblStimOrig,
      exTblUnsOrig = exTblUnsOrig,
      exTblStimNoMin = exTblStimNoMin,
      bias = bias,
      stage = stage,
      chnl = chnl
    )
    .intSaveNm(
      "locDetailCondition",
      locDetailCondition,
      ind,
      stageChnl,
      pathProject
    )
    return(objOut)
  }

  dataThreshold <- .getCpUnsLocGetCpDataThreshold(
    dataMod = dataMod,
    exTblStimOrig = exTblStimOrig,
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsBias = exTblUnsBias,
    exTblUnsOrig = exTblUnsOrig,
    bias = bias,
    stage = stage,
    pathProject = pathProject
  )
  .intSave(ind, stageChnl, pathProject, dataThreshold)
  cpObj <- .getCpUnsLocGetCpActual(
    dataThreshold = dataThreshold,
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsBias = exTblUnsBias,
    cpMin = cpMin,
    stage = stage
  )
  locDetailCondition <- .getCpUnsLocConditionDetailRow(
    cpObj = cpObj,
    dataThreshold = dataThreshold,
    exTblStimOrig = exTblStimOrig,
    exTblUnsOrig = exTblUnsOrig,
    exTblStimNoMin = exTblStimNoMin,
    bias = bias,
    stage = stage,
    chnl = chnl
  )
  .intSaveNm(
    "locDetailCondition",
    locDetailCondition,
    ind,
    stageChnl,
    pathProject
  )
  .intSave(ind, stageChnl, pathProject, cpObj$cp)
  .debug("Completed loc gate for single sample") # nolint
  cpObj
}

#' @keywords internal
.getCpUnsLocGetCpEnsureMeta <- function(
  obj,
  cpMin,
  exTblStimNoMin,
  exTblUnsBias,
  stage,
  reason
) {
  if (is.list(obj) && "cp" %in% names(obj)) {
    obj$locGenerated <- obj$locGenerated %||% FALSE
    obj$locGeneratedDirect <- obj$locGeneratedDirect %||% FALSE
    obj$locSource <- obj$locSource %||% "not_calculated"
    obj$locReason <- obj$locReason %||% reason
    obj$pList <- obj$pList %||% .getCpUnsLocPListEmpty()
    return(obj)
  }
  .getCpUnsLocConditionOut(
    cp = .getCpUnsLocConditionCpNonLoc(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias
    ),
    locGenerated = FALSE,
    locGeneratedDirect = FALSE,
    locSource = "not_calculated",
    locReason = reason
  )
}


#' @keywords internal
.getCpUnsLocGetCpTrimBeforeThreshold <- function(
  dataMod,
  exTblStimNoMin,
  exTblUnsBias,
  cpMin,
  stage,
  chnlSettings
) {
  info <- list(
    applied = FALSE,
    reason = "not_trimmed"
  )

  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    info$reason <- "no_data_mod"
    return(list("dataMod" = dataMod, "cp" = NULL, "info" = info))
  }

  dataMod <- dataMod[order(.getCut(dataMod)), , drop = FALSE]
  probCol <- .getCpUnsLocGetCpTrimProbCol(dataMod, chnlSettings)
  probVec <- .getCpUnsLocGetCpTrimProbVec(dataMod, probCol)
  maxProb <- suppressWarnings(max(probVec, na.rm = TRUE))
  minPeakProb <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locMinPeakProb",
    0.25
  )

  info$probCol <- probCol
  info$maxProb <- maxProb
  info$minPeakProb <- minPeakProb

  if (!is.finite(maxProb) || maxProb < minPeakProb) {
    info$applied <- TRUE
    info$reason <- "max_response_probability_below_minimum"
    return(list(
      "dataMod" = dataMod[0, , drop = FALSE],
      "cp" = .getCpUnsLocConditionCpNonLoc(
        cpMin = cpMin,
        exTblStimNoMin = exTblStimNoMin,
        exTblUnsBias = exTblUnsBias
      ),
      "info" = info
    ))
  }

  antiObj <- .getCpUnsLocGetCpTrimAntimode(
    dataMod = dataMod,
    exTblStimNoMin = exTblStimNoMin,
    chnlSettings = chnlSettings,
    probCol = probCol
  )
  dataMod <- antiObj$dataMod
  info$antimode <- antiObj$info

  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    info$applied <- TRUE
    info$reason <- "all_cells_removed_by_antimode_trim"
    return(list(
      "dataMod" = dataMod,
      "cp" = .getCpUnsLocConditionCpNonLoc(
        cpMin = cpMin,
        exTblStimNoMin = exTblStimNoMin,
        exTblUnsBias = exTblUnsBias
      ),
      "info" = info
    ))
  }

  flatObj <- .getCpUnsLocGetCpTrimFlatLeft(
    dataMod = dataMod,
    exTblStimNoMin = exTblStimNoMin,
    chnlSettings = chnlSettings,
    probCol = probCol
  )
  dataMod <- flatObj$dataMod
  info$flatLeft <- flatObj$info

  info$applied <- isTRUE(antiObj$info$applied) ||
    isTRUE(flatObj$info$applied)
  if (isTRUE(info$applied)) {
    info$reason <- "trimmed_before_threshold"
  }

  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    info$reason <- "all_cells_removed_by_left_flat_trim"
    return(list(
      "dataMod" = dataMod,
      "cp" = .getCpUnsLocConditionCpNonLoc(
        cpMin = cpMin,
        exTblStimNoMin = exTblStimNoMin,
        exTblUnsBias = exTblUnsBias
      ),
      "info" = info
    ))
  }

  list("dataMod" = dataMod, "cp" = NULL, "info" = info)
}

#' @keywords internal
.getCpUnsLocGetCpTrimSetting <- function(
  chnlSettings,
  nm,
  default
) {
  if (!is.null(chnlSettings) && !is.null(chnlSettings[[nm]])) {
    return(chnlSettings[[nm]])
  }
  default
}

#' @keywords internal
.getCpUnsLocGetCpTrimAsFiniteNumeric <- function(x, positive = FALSE) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (isTRUE(positive)) {
    x <- x[x > 0]
  }
  if (length(x) == 0L) {
    return(NA_real_)
  }
  min(x)
}

#' @keywords internal
.getCpUnsLocGetCpTrimProbCol <- function(dataMod, chnlSettings) {
  probCol <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locProbCol",
    "pred"
  )
  if (probCol %in% names(dataMod)) {
    return(probCol)
  }
  if ("pred" %in% names(dataMod)) {
    return("pred")
  }
  "probSmooth"
}

#' @keywords internal
.getCpUnsLocGetCpTrimProbVec <- function(dataMod, probCol) {
  probVec <- dataMod[[probCol]]
  probVec <- pmin(1, pmax(0, probVec))
  probVec[!is.finite(probVec)] <- NA_real_
  probVec
}

#' @keywords internal
.getCpUnsLocGetCpTrimAntimode <- function(
  dataMod,
  exTblStimNoMin,
  chnlSettings,
  probCol
) {
  info <- list(
    applied = FALSE,
    reason = "not_multimodal"
  )

  exprVec <- .getCut(exTblStimNoMin)
  exprVec <- exprVec[is.finite(exprVec)]
  if (length(exprVec) < 5L || length(unique(exprVec)) < 3L) {
    info$reason <- "too_few_expression_values"
    return(list("dataMod" = dataMod, "info" = info))
  }

  dipP <- .getCpUnsLocGetCpTrimDipP(exprVec)
  dipAlpha <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locDipAlpha",
    0.2
  )
  info$dipP <- dipP
  info$dipAlpha <- dipAlpha

  if (!is.finite(dipP) || dipP >= dipAlpha) {
    info$reason <- "dip_test_not_significant"
    return(list("dataMod" = dataMod, "info" = info))
  }

  densObj <- .getCpUnsLocGetCpTrimDensity(
    exprVec = exprVec,
    chnlSettings = chnlSettings
  )
  if (is.null(densObj)) {
    info$reason <- "density_failed"
    return(list("dataMod" = dataMod, "info" = info))
  }

  minObj <- .getCpUnsLocGetCpTrimAntimodes(
    densObj = densObj,
    chnlSettings = chnlSettings
  )
  info$bw <- attr(densObj, "bw")
  info$peakHeight <- minObj$peakHeight
  info$antimodeHeightFrac <- minObj$heightFrac
  info$antimodeX <- minObj$antimodeX

  if (length(minObj$antimodeX) == 0L) {
    info$reason <- "no_deep_antimodes"
    return(list("dataMod" = dataMod, "info" = info))
  }

  regionObj <- .getCpUnsLocGetCpTrimAntimodeRegions(
    dataMod = dataMod,
    antimodeX = minObj$antimodeX,
    probCol = probCol,
    chnlSettings = chnlSettings
  )
  info <- c(info, regionObj$info)
  if (!isTRUE(regionObj$info$applied)) {
    return(list("dataMod" = dataMod, "info" = info))
  }

  list("dataMod" = regionObj$dataMod, "info" = info)
}

#' @keywords internal
.getCpUnsLocGetCpTrimDipP <- function(exprVec) {
  if (!requireNamespace("diptest", quietly = TRUE)) {
    return(NA_real_)
  }
  dipObj <- try(
    suppressWarnings(diptest::dip.test(exprVec)),
    silent = TRUE
  )
  if (inherits(dipObj, "try-error")) {
    return(NA_real_)
  }
  dipObj$p.value
}

#' @keywords internal
.getCpUnsLocGetCpTrimDensity <- function(exprVec, chnlSettings) {
  bw <- .getCpUnsLocGetCpTrimDensityBw(
    exprVec = exprVec,
    chnlSettings = chnlSettings
  )
  densObj <- try(
    suppressWarnings(stats::density(exprVec, bw = bw, n = 512)),
    silent = TRUE
  )
  if (inherits(densObj, "try-error")) {
    return(NULL)
  }
  attr(densObj, "bw") <- densObj$bw
  densObj
}

#' @keywords internal
.getCpUnsLocGetCpTrimDensityBw <- function(exprVec, chnlSettings) {
  hpiBw <- try(
    suppressWarnings(ks::hpi(x = exprVec)),
    silent = TRUE
  )
  hpiBw <- if (inherits(hpiBw, "try-error")) NA_real_ else hpiBw
  hpiBw <- .getCpUnsLocGetCpTrimAsFiniteNumeric(hpiBw, positive = TRUE)

  preCalcBw <- .getCpUnsLocGetCpTrimAsFiniteNumeric(
    .getCpUnsLocGetCpTrimSetting(
      chnlSettings,
      "bwCluster",
      NA_real_
    ),
    positive = TRUE
  )
  if (!is.finite(preCalcBw)) {
    preCalcBw <- .getCpUnsLocGetCpTrimAsFiniteNumeric(
      .getCpUnsLocGetCpTrimSetting(
        chnlSettings,
        "bw",
        NA_real_
      ),
      positive = TRUE
    )
  }
  if (!is.finite(preCalcBw)) {
    preCalcBw <- .getCpUnsLocGetCpTrimAsFiniteNumeric(
      .getCpUnsLocGetCpTrimSetting(
        chnlSettings,
        "bwMax",
        NA_real_
      ),
      positive = TRUE
    )
  }

  bwVec <- c(hpiBw, preCalcBw / 4)
  bwVec <- bwVec[is.finite(bwVec) & bwVec > 0]
  if (length(bwVec) == 0L) {
    return("nrd0")
  }
  min(bwVec)
}

#' @keywords internal
.getCpUnsLocGetCpTrimAntimodes <- function(densObj, chnlSettings) {
  y <- densObj$y
  x <- densObj$x
  heightFrac <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locAntimodeHeightFrac",
    1 / 6
  )
  peakHeight <- max(y, na.rm = TRUE)
  if (!is.finite(peakHeight) || length(y) < 3L) {
    return(list(
      "antimodeX" = numeric(0),
      "peakHeight" = peakHeight,
      "heightFrac" = heightFrac
    ))
  }

  idx <- seq.int(2L, length(y) - 1L)
  minIdx <- idx[y[idx] <= y[idx - 1L] & y[idx] < y[idx + 1L]]
  minIdx <- minIdx[y[minIdx] < peakHeight * heightFrac]

  list(
    "antimodeX" = sort(unique(x[minIdx])),
    "peakHeight" = peakHeight,
    "heightFrac" = heightFrac
  )
}

#' @keywords internal
.getCpUnsLocGetCpTrimAntimodeRegions <- function(
  dataMod,
  antimodeX,
  probCol,
  chnlSettings
) {
  info <- list(
    applied = FALSE,
    reason = "no_low_left_antimode_region"
  )
  x <- .getCut(dataMod)
  probVec <- .getCpUnsLocGetCpTrimProbVec(dataMod, probCol)
  region <- findInterval(
    x,
    c(-Inf, antimodeX, Inf),
    rightmost.closed = TRUE
  )
  regionId <- sort(unique(region))
  regionMean <- purrr::map_dbl(regionId, function(id) {
    mean(probVec[region == id], na.rm = TRUE)
  })
  regionN <- purrr::map_int(regionId, function(id) sum(region == id))
  regionXMin <- purrr::map_dbl(regionId, function(id) min(x[region == id]))
  regionXMax <- purrr::map_dbl(regionId, function(id) max(x[region == id]))

  if (length(regionId) == 0L || all(!is.finite(regionMean))) {
    info$reason <- "no_regions_with_probability"
    return(list("dataMod" = dataMod, "info" = info))
  }

  peakRegion <- regionId[which.max(regionMean)]
  peakRegionProb <- max(regionMean, na.rm = TRUE)
  lowRel <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locAntimodeLowRel",
    0.25
  )
  lowAbs <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locAntimodeLowAbs",
    0.15
  )
  dropLgl <- regionId < peakRegion &
    is.finite(regionMean) &
    (regionMean < peakRegionProb * lowRel | regionMean < lowAbs)
  dropLgl[is.na(dropLgl)] <- FALSE
  dropRegion <- regionId[dropLgl]

  regionSummary <- tibble::tibble(
    region = regionId,
    meanProb = regionMean,
    n = regionN,
    xMin = regionXMin,
    xMax = regionXMax,
    drop = regionId %in% dropRegion
  )

  info$regionSummary <- regionSummary
  info$peakRegion <- peakRegion
  info$peakRegionProb <- peakRegionProb
  info$lowRel <- lowRel
  info$lowAbs <- lowAbs
  info$dropRegion <- dropRegion

  if (length(dropRegion) == 0L) {
    return(list("dataMod" = dataMod, "info" = info))
  }

  keep <- !region %in% dropRegion
  info$applied <- TRUE
  info$reason <- "dropped_low_probability_left_antimode_regions"
  list("dataMod" = dataMod[keep, , drop = FALSE], "info" = info)
}


#' @keywords internal
.getCpUnsLocGetCpTrimFlatLeftDerivData <- function(
  probData,
  dataMod,
  probCol
) {
  derivData <- .getCpUnsLocGetCpTrimFlatLeftDerivFitted(
    probData = probData,
    dataMod = dataMod,
    probCol = probCol
  )
  if (!is.null(derivData)) {
    return(derivData)
  }
  .getCpUnsLocGetCpTrimFlatLeftDerivObserved(probData)
}

#' @keywords internal
.getCpUnsLocGetCpTrimFlatLeftDerivFitted <- function(
  probData,
  dataMod,
  probCol
) {
  if (!identical(probCol, "pred")) {
    return(NULL)
  }
  derivTbl <- attr(dataMod, "locProbDerivTbl")
  if (
    is.null(derivTbl) ||
      !all(c("x", "deriv") %in% names(derivTbl)) ||
      nrow(derivTbl) < 2L
  ) {
    return(NULL)
  }
  derivTbl <- derivTbl |>
    dplyr::filter(is.finite(.data$x), is.finite(.data$deriv)) |>
    dplyr::arrange(.data$x)
  if (nrow(derivTbl) < 2L || diff(range(derivTbl$x)) <= 0) {
    return(NULL)
  }

  derivVec <- try(
    stats::approx(
      x = derivTbl$x,
      y = derivTbl$deriv,
      xout = probData$x,
      rule = 2
    )$y,
    silent = TRUE
  )
  if (inherits(derivVec, "try-error")) {
    return(NULL)
  }
  tibble::tibble(
    x = probData$x,
    deriv = pmax(0, as.numeric(derivVec)),
    source = "fitted_pred_derivative"
  )
}

#' @keywords internal
.getCpUnsLocGetCpTrimFlatLeftDerivObserved <- function(probData) {
  if (nrow(probData) < 4L || diff(range(probData$x)) <= 0) {
    return(NULL)
  }
  deriv <- diff(probData$prob) / diff(probData$x)
  deriv[!is.finite(deriv)] <- 0
  deriv <- pmax(0, deriv)
  tibble::tibble(
    x = probData$x[-nrow(probData)],
    deriv = deriv,
    source = "observed_probability_slope"
  )
}


#' @keywords internal
.getCpUnsLocGetCpTrimMarginalLeft <- function(
  dataMod,
  chnlSettings,
  probCol,
  startX
) {
  info <- list(
    applied = FALSE,
    reason = "marginal_left_trim_not_run",
    startX = startX
  )

  if (!is.data.frame(dataMod) || nrow(dataMod) < 4L) {
    info$reason <- "too_few_cells_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  dataMod <- dataMod[order(.getCut(dataMod)), , drop = FALSE]
  x <- suppressWarnings(as.numeric(.getCut(dataMod)))
  probVec <- .getCpUnsLocGetCpTrimProbVec(dataMod, probCol)
  dm <- tibble::tibble(x = x, prob = probVec) |>
    dplyr::filter(is.finite(.data$x), is.finite(.data$prob)) |>
    dplyr::mutate(prob = pmin(1, pmax(0, .data$prob))) |>
    dplyr::arrange(.data$x)

  if (
    nrow(dm) < 4L ||
      !is.finite(startX) ||
      diff(range(dm$x, na.rm = TRUE)) <= 0
  ) {
    info$reason <- "insufficient_data_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  xRight <- dm$x[dm$x >= startX]
  if (length(xRight) == 0L) {
    info$reason <- "no_cells_right_of_marginal_start"
    return(list("dataMod" = dataMod, "info" = info))
  }

  refQuantile <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locMarginalRefQuantile",
    0.75
  )
  refQuantile <- suppressWarnings(as.numeric(refQuantile[1]))
  if (!is.finite(refQuantile) || refQuantile <= 0 || refQuantile > 1) {
    refQuantile <- 0.75
  }

  xRefHi <- suppressWarnings(stats::quantile(
    xRight,
    probs = refQuantile,
    na.rm = TRUE,
    names = FALSE
  ))
  if (!is.finite(xRefHi) || xRefHi <= startX) {
    xRefHi <- max(xRight, na.rm = TRUE)
  }
  if (!is.finite(xRefHi) || xRefHi <= startX) {
    info$reason <- "right_reference_interval_too_short"
    return(list("dataMod" = dataMod, "info" = info))
  }

  breaks <- .getCpUnsLocGetCpTrimMarginalBreaks(
    dataMod = dataMod,
    xMin = min(dm$x, na.rm = TRUE),
    xMax = max(dm$x, na.rm = TRUE),
    startX = startX,
    xRefHi = xRefHi
  )
  if (length(breaks) < 2L) {
    info$reason <- "insufficient_bins_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  binTbl <- tibble::tibble(
    left = breaks[-length(breaks)],
    right = breaks[-1]
  ) |>
    dplyr::filter(.data$right > .data$left)

  rightRefBins <- binTbl |>
    dplyr::filter(
      .data$right > .env$startX,
      .data$left < .env$xRefHi
    )
  refNBin <- nrow(rightRefBins)
  if (refNBin == 0L) {
    info$reason <- "no_right_reference_bins_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  rightMaskForBinDensity <- dm$x >= startX & dm$x <= xRefHi
  rightMaskForPurity <- dm$x >= startX
  refNCellBinRegion <- sum(rightMaskForBinDensity, na.rm = TRUE)
  refNCell <- sum(rightMaskForPurity, na.rm = TRUE)
  refExpectedResp <- sum(dm$prob[rightMaskForPurity], na.rm = TRUE)
  if (refNCell == 0L || !is.finite(refExpectedResp)) {
    info$reason <- "no_right_reference_cells_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  refCellsPerBin <- refNCellBinRegion / refNBin
  refPurity <- refExpectedResp / refNCell
  if (
    !is.finite(refCellsPerBin) ||
      refCellsPerBin < 0 ||
      !is.finite(refPurity) ||
      refPurity <= 0
  ) {
    info$reason <- "invalid_right_reference_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  cellBinRatio <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locMarginalCellBinRatio",
    2
  )
  cellBinRatio <- suppressWarnings(as.numeric(cellBinRatio[1]))
  if (!is.finite(cellBinRatio) || cellBinRatio <= 0) {
    cellBinRatio <- 2
  }

  purityRel <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locMarginalPurityRel",
    0.5
  )
  purityRel <- suppressWarnings(as.numeric(purityRel[1]))
  if (!is.finite(purityRel) || purityRel < 0 || purityRel > 1) {
    purityRel <- 0.5
  }

  maxLeftBinCells <- cellBinRatio * refCellsPerBin
  minLeftBinPurity <- purityRel * refPurity

  leftBins <- binTbl |>
    dplyr::filter(.data$right <= .env$startX) |>
    dplyr::arrange(dplyr::desc(.data$right))
  if (nrow(leftBins) == 0L) {
    info$reason <- "no_left_bins_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  currentCut <- startX
  scanList <- vector("list", nrow(leftBins))
  stopReason <- "accepted_all_left_bins"

  for (i in seq_len(nrow(leftBins))) {
    left <- leftBins$left[[i]]
    right <- leftBins$right[[i]]
    mask <- dm$x >= left & dm$x < right
    nCell <- sum(mask, na.rm = TRUE)
    expectedResp <- sum(dm$prob[mask], na.rm = TRUE)
    purity <- if (nCell > 0L) expectedResp / nCell else NA_real_

    acceptBin <- nCell == 0L ||
      (is.finite(purity) &&
        nCell <= maxLeftBinCells &&
        purity >= minLeftBinPurity)

    scanList[[i]] <- tibble::tibble(
      left = left,
      right = right,
      nCell = nCell,
      expectedResp = expectedResp,
      purity = purity,
      refNBin = refNBin,
      refCellsPerBin = refCellsPerBin,
      refPurity = refPurity,
      maxLeftBinCells = maxLeftBinCells,
      minLeftBinPurity = minLeftBinPurity,
      relCellsPerBin = ifelse(refCellsPerBin > 0, nCell / refCellsPerBin, Inf),
      relPurity = purity / refPurity,
      acceptBin = acceptBin
    )

    if (!isTRUE(acceptBin)) {
      stopReason <- "stopped_at_low_purity_or_high_cell_bin"
      break
    }
    currentCut <- left
  }

  scanTbl <- dplyr::bind_rows(scanList)
  finalStartX <- currentCut
  keep <- is.finite(x) & x >= finalStartX

  info$applied <- sum(!keep, na.rm = TRUE) > 0L
  info$reason <- if (isTRUE(info$applied)) {
    "dropped_left_region_by_marginal_purity_cell_trim"
  } else {
    "marginal_purity_cell_trim_kept_all_cells"
  }
  info$stopReason <- stopReason
  info$finalStartX <- finalStartX
  info$xRefHi <- xRefHi
  info$refQuantile <- refQuantile
  info$refNBin <- refNBin
  info$refNCell <- refNCell
  info$refNCellBinRegion <- refNCellBinRegion
  info$refExpectedResp <- refExpectedResp
  info$refCellsPerBin <- refCellsPerBin
  info$refPurity <- refPurity
  info$cellBinRatio <- cellBinRatio
  info$purityRel <- purityRel
  info$maxLeftBinCells <- maxLeftBinCells
  info$minLeftBinPurity <- minLeftBinPurity
  info$scanTbl <- scanTbl

  list("dataMod" = dataMod[keep, , drop = FALSE], "info" = info)
}

#' @keywords internal
.getCpUnsLocGetCpTrimMarginalBreaks <- function(
  dataMod,
  xMin,
  xMax,
  startX,
  xRefHi
) {
  binVec <- attr(dataMod, "binVec")
  if (is.null(binVec)) {
    binVec <- stats::pretty(.getCut(dataMod), n = 80)
  }
  binVec <- suppressWarnings(as.numeric(binVec))
  binVec <- binVec[is.finite(binVec)]
  breaks <- sort(unique(c(xMin, xMax, startX, xRefHi, binVec)))
  breaks <- breaks[
    is.finite(breaks) &
      breaks >= xMin &
      breaks <= xMax
  ]
  sort(unique(breaks))
}

#' @keywords internal
.getCpUnsLocGetCpTrimFlatLeft <- function(
  dataMod,
  exTblStimNoMin,
  chnlSettings,
  probCol
) {
  info <- list(
    applied = FALSE,
    reason = "no_derivative_based_left_trim"
  )

  if (!is.data.frame(dataMod) || nrow(dataMod) < 4L) {
    info$reason <- "too_few_cells_for_derivative_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  dataMod <- dataMod[order(.getCut(dataMod)), , drop = FALSE]
  x <- suppressWarnings(as.numeric(.getCut(dataMod)))
  probVec <- .getCpUnsLocGetCpTrimProbVec(dataMod, probCol)
  probData <- tibble::tibble(x = x, prob = probVec) |>
    dplyr::filter(is.finite(x), is.finite(prob)) |>
    dplyr::group_by(x) |>
    dplyr::summarise(prob = mean(prob), .groups = "drop") |>
    dplyr::arrange(x)

  if (nrow(probData) < 4L || diff(range(probData$x)) <= 0) {
    info$reason <- "insufficient_unique_x_for_derivative_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  derivData <- .getCpUnsLocGetCpTrimFlatLeftDerivData(
    probData = probData,
    dataMod = dataMod,
    probCol = probCol
  )
  if (is.null(derivData) || nrow(derivData) == 0L) {
    info$reason <- "no_probability_derivative_available"
    return(list("dataMod" = dataMod, "info" = info))
  }

  deriv <- derivData$deriv
  deriv[!is.finite(deriv)] <- 0
  deriv <- pmax(0, deriv)
  peakDeriv <- max(deriv, na.rm = TRUE)
  if (!is.finite(peakDeriv) || peakDeriv <= 0) {
    info$reason <- "no_positive_probability_derivative"
    return(list("dataMod" = dataMod, "info" = info))
  }

  hardDerivFrac <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locFlatHardDerivFrac",
    1 / 4
  )
  hardDerivFrac <- suppressWarnings(as.numeric(hardDerivFrac[1]))
  if (!is.finite(hardDerivFrac) || hardDerivFrac < 0 || hardDerivFrac > 1) {
    hardDerivFrac <- 1 / 4
  }

  hardDerivThreshold <- peakDeriv * hardDerivFrac
  hardIncIdx <- which(deriv >= hardDerivThreshold)
  if (length(hardIncIdx) == 0L) {
    info$reason <- "no_derivative_above_hard_flat_threshold"
    return(list("dataMod" = dataMod, "info" = info))
  }
  hardStartX <- derivData$x[min(hardIncIdx)]

  derivFrac <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locFlatDerivFrac",
    1 / 2
  )
  derivFrac <- suppressWarnings(as.numeric(derivFrac[1]))
  if (!is.finite(derivFrac) || derivFrac < 0 || derivFrac > 1) {
    derivFrac <- 1 / 2
  }

  derivThreshold <- peakDeriv * derivFrac
  incIdx <- which(deriv >= derivThreshold)
  if (length(incIdx) == 0L) {
    info$reason <- "no_derivative_above_marginal_flat_threshold"
    return(list("dataMod" = dataMod, "info" = info))
  }

  startX <- derivData$x[min(incIdx)]
  startX <- max(startX, hardStartX, na.rm = TRUE)

  hardKeep <- is.finite(x) & x >= hardStartX
  if (sum(hardKeep, na.rm = TRUE) == 0L) {
    info$reason <- "all_cells_removed_by_hard_derivative_trim"
    return(list("dataMod" = dataMod[0, , drop = FALSE], "info" = info))
  }

  dataModHard <- dataMod[hardKeep, , drop = FALSE]
  marginalObj <- .getCpUnsLocGetCpTrimMarginalLeft(
    dataMod = dataModHard,
    chnlSettings = chnlSettings,
    probCol = probCol,
    startX = startX
  )

  hardDropped <- sum(!hardKeep, na.rm = TRUE) > 0L
  marginalApplied <- isTRUE(marginalObj$info$applied)

  info$applied <- hardDropped || marginalApplied
  info$reason <- if (isTRUE(info$applied)) {
    "dropped_left_region_by_derivative_or_marginal_trim"
  } else {
    "derivative_and_marginal_trim_kept_all_cells"
  }
  info$hardStartX <- hardStartX
  info$startX <- startX
  info$marginalStartX <- startX
  info$hardDerivFrac <- hardDerivFrac
  info$derivFrac <- derivFrac
  info$hardDerivThreshold <- hardDerivThreshold
  info$derivThreshold <- derivThreshold
  info$peakDeriv <- peakDeriv
  info$derivSource <- unique(derivData$source)[1]
  info$nDerivPoint <- nrow(derivData)
  info$nHardDropped <- sum(!hardKeep, na.rm = TRUE)
  info$marginal <- marginalObj$info

  list("dataMod" = marginalObj$dataMod, "info" = info)
}

#' @keywords internal
.getCpUnsLocGetCpDataThreshold <- function(
  dataMod,
  exTblStimOrig,
  exTblStimNoMin,
  exTblUnsBias,
  exTblUnsOrig,
  bias,
  pathProject,
  stage
) {
  # final filtering step
  dataCount <- .getCpUnsLocGetCpDataThresholdCount(dataMod)
  probBsEst <- .getCpUnsLocGetCpDataThresholdPropBsEst(
    dataCount = dataCount,
    exTblStimOrig = exTblStimOrig
  )
  .intSaveNm(
    "probBsEstConditionRaw",
    probBsEst,
    .getInd(exTblStimNoMin),
    file.path(stage, .getCpUnsLocGetChnl(exTblStimNoMin)),
    pathProject
  )
  .getCpUnsLocGetCpDataThresholdActual(
    dataCount = dataCount,
    propBsEst = probBsEst,
    exTblStimOrig = exTblStimOrig,
    exTblUnsBias = exTblUnsBias,
    exTblUnsOrig = exTblUnsOrig,
    bias = bias
  )
}

#' @keywords internal
.getCpUnsLocGetCpDataThresholdCount <- function(dataMod) {
  if (nrow(dataMod) == 1L) {
    minVal <- min(.getCut(dataMod)) - 1
  } else {
    minVal <- min(.getCut(dataMod))
  }
  dataMod <- dataMod[.getCut(dataMod) > minVal, ]
  dataMod <- dataMod[order(.getCut(dataMod)), ]
  dataMod |>
    dplyr::mutate(nRow = seq_len(dplyr::n())) |>
    dplyr::filter(cumsum(pred > probSmooth) != nRow) |> # nolint
    dplyr::select(-nRow)
}

#' @keywords internal
.getCpUnsLocGetCpDataThresholdPropBsEst <- function(
  dataCount,
  exTblStimOrig
) {
  sum(dataCount$pred) / nrow(exTblStimOrig)
}

#' @keywords internal
.getCpUnsLocGetCpDataThresholdActual <- function(
  dataCount,
  propBsEst,
  exTblStimOrig,
  exTblUnsBias,
  exTblUnsOrig,
  bias
) {
  dataCount <- dataCount[order(.getCut(dataCount)), ]
  propStimVec <- purrr::map_dbl(.getCut(dataCount), function(x) {
    sum(.getCut(exTblStimOrig) >= x) / nrow(exTblStimOrig)
  })
  propUnsVec <- purrr::map_dbl(.getCut(dataCount), function(x) {
    sum(.getCut(exTblUnsOrig) >= x) / nrow(exTblUnsOrig)
  })
  dataCount |>
    dplyr::mutate(
      propStim = propStimVec,
      propUns = propUnsVec
    ) |>
    dplyr::mutate(
      propBs = propStim - propUns, # nolint
      propBsDiff = propBs - propBsEst # nolint
    )
}

#' @keywords internal
.getCpUnsLocGetCpActual <- function(
  dataThreshold,
  exTblStimNoMin,
  exTblUnsBias,
  cpMin,
  stage
) {
  if (nrow(dataThreshold) == 0L) {
    return(.getCpUnsLocConditionCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "Too few responding cells"
    ))
  }
  dataThreshold <- dataThreshold |>
    dplyr::filter(abs(propBsDiff) == min(abs(propBsDiff))) |> # nolint
    dplyr::slice(1) |>
    .getCut()

  if (!is.finite(dataThreshold)) {
    return(.getCpUnsLocConditionCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "No finite local-FDR threshold"
    ))
  }

  .getCpUnsLocConditionOut(
    cp = dataThreshold,
    locGenerated = TRUE,
    locGeneratedDirect = TRUE,
    locSource = "direct",
    locReason = "local_fdr_threshold_selected"
  )
}


#' @keywords internal
.createCombinedIdentifier <- function(indStim) {
  if (is.null(indStim) || length(indStim) == 0) {
    "empty_batch"
  } else {
    paste(indStim, collapse = "_")
  }
}

#' @keywords internal
.getCpUnsLocOutput <- function(
  cpUnsLocObjList,
  indUns,
  indStim,
  stage,
  pathProject,
  chnl,
  exListOrig
) {
  stageChnl <- file.path(stage, chnl)
  cpVec <- .getCpUnsLocSampleCpRep(
    stage = stage,
    cpUnsLocObjList = cpUnsLocObjList,
    indUns = indUns,
    indStim = indStim,
    pathProject = pathProject,
    chnl = chnl
  )
  indCombined <- .createCombinedIdentifier(indStim)
  .intSave(indCombined, stageChnl, pathProject, cpVec)
  .intSaveNm(
    "cpUnsLocMeta",
    .getCpUnsLocMetaFromCp(cpVec),
    indCombined,
    stageChnl,
    pathProject
  )
  locDetailSample <- .getCpUnsLocSampleDetailTbl(
    cpVec = cpVec,
    exListOrig = exListOrig,
    indUns = indUns,
    indStim = indStim,
    stage = stage,
    chnl = chnl
  )
  .intSaveNm(
    "locDetailSample",
    locDetailSample,
    indCombined,
    stageChnl,
    pathProject
  )
  .debug("done getting loc gate at sample level") # nolint
  list(
    "loc" = cpVec,
    "pList" = list()
  )
}


#' @keywords internal
.getCpUnsLocSampleCpRep <- function(
  stage,
  cpUnsLocObjList,
  indUns,
  indStim,
  pathProject,
  chnl
) {
  .debug("Possibly re-using calculated cutpoints") # nolint
  indCombined <- .createCombinedIdentifier(indStim)
  stageChnl <- file.path(stage, chnl)

  cpVec <- purrr::map_dbl(cpUnsLocObjList, ~ .x[["cp"]])
  meta <- .getCpUnsLocSampleMeta(
    cpUnsLocObjList = cpUnsLocObjList,
    ind = names(cpUnsLocObjList)
  )
  .intSaveNm(
    "cpVecBeforeRep",
    cpVec,
    indCombined,
    stageChnl,
    pathProject
  )
  .intSaveNm(
    "cpMetaBeforeRep",
    meta,
    indCombined,
    stageChnl,
    pathProject
  )

  if (length(cpVec) != length(indStim)) {
    .intSaveNm(
      "prejoinedCpUsed",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpVec <- stats::setNames(
      rep(cpVec, length(indStim)),
      indStim
    )
    meta <- .getCpUnsLocSampleMetaRep(
      meta = meta,
      indStim = indStim,
      source = "prejoin"
    )
  } else {
    .intSaveNm(
      "individualCpUsed",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpVec <- stats::setNames(cpVec, indStim)
    meta$ind <- as.character(indStim)
  }
  .intSaveNm(
    "cpVecAfterRep",
    cpVec,
    indCombined,
    stageChnl,
    pathProject
  )
  .intSaveNm(
    "cpMetaAfterRep",
    meta,
    indCombined,
    stageChnl,
    pathProject
  )

  locGenerated <- meta$locGenerated %in% TRUE
  if (any(locGenerated, na.rm = TRUE)) {
    .intSaveNm(
      "addingUnsCp",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpUns <- mean(cpVec[locGenerated], na.rm = TRUE)
    metaUns <- tibble::tibble(
      ind = as.character(indUns),
      locGenerated = TRUE,
      locGeneratedDirect = FALSE,
      locSource = "unstim_summary",
      locReason = "mean_of_generated_local_fdr_thresholds"
    )
    cpVec <- c(cpVec, stats::setNames(cpUns, indUns))
    meta <- dplyr::bind_rows(meta, metaUns)
  } else {
    .intSaveNm(
      "addNaForUnsCp",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpVec <- c(cpVec, stats::setNames(NA_real_, indUns))
    meta <- dplyr::bind_rows(
      meta,
      tibble::tibble(
        ind = as.character(indUns),
        locGenerated = FALSE,
        locGeneratedDirect = FALSE,
        locSource = "unstim_summary",
        locReason = "no_generated_local_fdr_thresholds"
      )
    )
  }
  .intSaveNm(
    "cpVecAfterUnsAdded",
    cpVec,
    indCombined,
    stageChnl,
    pathProject
  )
  .intSaveNm(
    "cpMetaAfterUnsAdded",
    meta,
    indCombined,
    stageChnl,
    pathProject
  )

  .getCpUnsLocCpAttachMeta(cpVec, meta)
}

#' @keywords internal
.getCpUnsLocSampleMeta <- function(cpUnsLocObjList, ind) {
  tibble::tibble(
    ind = as.character(ind),
    locGenerated = purrr::map_lgl(
      cpUnsLocObjList,
      ~ .x[["locGenerated"]] %||% FALSE
    ),
    locGeneratedDirect = purrr::map_lgl(
      cpUnsLocObjList,
      ~ .x[["locGeneratedDirect"]] %||% FALSE
    ),
    locSource = purrr::map_chr(
      cpUnsLocObjList,
      ~ .x[["locSource"]] %||% "not_calculated"
    ),
    locReason = purrr::map_chr(
      cpUnsLocObjList,
      ~ .x[["locReason"]] %||% NA_character_
    )
  )
}

#' @keywords internal
.getCpUnsLocSampleMetaRep <- function(meta, indStim, source) {
  if (nrow(meta) == 0L) {
    return(tibble::tibble(
      ind = as.character(indStim),
      locGenerated = FALSE,
      locGeneratedDirect = FALSE,
      locSource = "not_calculated",
      locReason = "no_local_fdr_objects"
    ))
  }
  metaOne <- meta[1, , drop = FALSE]
  metaOut <- metaOne[rep(1, length(indStim)), , drop = FALSE]
  metaOut$ind <- as.character(indStim)
  metaOut$locSource <- ifelse(
    metaOut$locGenerated %in% TRUE,
    source,
    metaOut$locSource
  )
  metaOut
}


#' @keywords internal
.thinDataMod <- function(dataMod, maxCellsPerBin = 20) {
  binVec <- attr(dataMod, "binVec")
  if (is.null(binVec)) {
    return(dataMod)
  }

  breaks <- c(-Inf, binVec[-length(binVec)] + diff(binVec) / 2, Inf)

  x_vals <- .getCut(dataMod)
  bin_indices <- findInterval(x_vals, breaks)

  # Sample down if a bin has too many cells
  idxMod <- unlist(lapply(split(seq_along(x_vals), bin_indices), function(idx) {
    if (length(idx) > maxCellsPerBin) {
      sample(idx, size = maxCellsPerBin, replace = FALSE)
    } else {
      idx
    }
  }))

  attr(dataMod, "idxMod") <- sort(unname(idxMod))

  dataMod
}
