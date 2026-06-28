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
    chnl = chnl
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
  .getCpUnsLocGetProbSmooth(dataMod, stage, pathProject, chnl)
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
  bw <- .getCpUnsLocGetDensRawDensitiesBw(
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold,
    bw = chnlSettings$bw,
    bwMin = chnlSettings$bwMin,
    bwMax = chnlSettings$bwMax,
    bwMtd = chnlSettings$bwMtd,
    bwAdj = chnlSettings$bwAdj,
    bwNcellMin = chnlSettings$bwNcellMin,
    bwNcellMax = chnlSettings$bwNcellMax
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
.getCpUnsLocGetDensRawDensitiesBw <- function(
  exTblStimThreshold,
  exTblUnsThreshold,
  bw,
  bwMin,
  bwMax,
  bwMtd,
  bwAdj,
  bwNcellMin,
  bwNcellMax
) {
  if (!is.null(bw)) {
    return(bw)
  }
  bwStim <- .getCpUnsLocGetDensRawDensitiesBwInit(
    .data = .getCut(exTblStimThreshold),
    bwMin = bwMin,
    bwMax = bwMax,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax
  )
  bwUns <- .getCpUnsLocGetDensRawDensitiesBwInit(
    .data = .getCut(exTblUnsThreshold),
    bwMin = bwMin,
    bwMax = bwMax,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax
  )
  min(bwUns, bwStim)
}

#' @keywords internal
.getCpUnsLocGetDensRawDensitiesBwInit <- function(
  .data,
  bwMin,
  bwMax,
  bwMtd,
  bwAdj,
  bwNcellMin,
  bwNcellMax
) {
  if (!is.null(bwNcellMin) && length(.data) < bwNcellMin) {
    iqrX <- diff(quantile(.data, c(0.75, 0.25), na.rm = TRUE))
    sdX <- abs(iqrX) / 1.5
    .data <- sample(.data, replace = TRUE, size = bwNcellMin) +
      rnorm(bwNcellMin, mean = 0, sd = sdX / 10)
  }
  if (!is.null(bwNcellMax) && length(.data) > bwNcellMax) {
    .data <- sample(.data, size = bwNcellMax, replace = FALSE)
  }
  bwCalc <- switch(
    bwMtd,
    "ndr0" = try(
      suppressWarnings(stats::density(.data, "nrd0")),
      silent = TRUE
    ),
    "sj" = try(suppressWarnings(stats::density(.data, "SJ")), silent = TRUE),
    try(
      suppressWarnings(
        ks::hpi(.data, deriv.order = as.numeric(gsub("hpi", "", bwMtd)))
      ),
      silent = TRUE
    )
  )
  if (inherits(bwCalc, "try-error")) {
    bwMin
  } else {
    max(bwMin, min(bwCalc * bwAdj, bwMax))
  }
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
  .data |>
    dplyr::mutate(
      yUns = purrr::map_dbl(
        xStim, # nolint
        function(marker) .interp(val = marker, x = unsX, y = unsY) # nolint
      )
    )
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

  probTblPos <- .getCpUnsLocProbTblFilter(
    exVecStimThreshold = exVecStimThreshold,
    exVecUnsThreshold = exVecUnsThreshold,
    probTbl = probTbl,
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
.getCpUnsLocProbTblFilter <- function(
  exVecStimThreshold,
  exVecUnsThreshold,
  probTbl,
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
  chnl
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

  predVec <- .getCpUnsLocGetProbSmoothActual(dataMod, stage)
  dataModOut <- dataMod |> dplyr::mutate(pred = predVec)
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
.getCpUnsLocGetProbSmoothActual <- function(dataMod, stage) {
  fit1 <- .getCpUnsLocGetProbSmoothActualFirst(dataMod, stage)
  .getCpUnsLocGetProbSmoothActualFirstResponse(
    fit = fit1,
    dataMod = dataMod,
    stage = stage
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
  stage
) {
  if (.getCpUnsLocGetProbSmoothActualCheck(fit, dataMod)) {
    .debug("Smoothed") # nolint
    return(
      .getCpUnsLocGetProbSmoothActualResponseSuccess(
        fit = fit,
        dataMod = dataMod
      )$pred
    )
  }
  .getCpUnsLocGetProbSmoothActualFirstResponseFailure(
    stage = stage,
    dataMod = dataMod
  )
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualCheck <- function(fit, dataMod) {
  if (inherits(fit, "try-error")) {
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
  dataMod
) {
  # predict on full dataset
  predVec <- predict(fit, newdata = dataMod, type = "response")
  meanAbsError <- mean(abs(predVec - dataMod$probSmooth))
  list("pred" = predVec, "meanAbsError" = meanAbsError)
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualFirstResponseFailure <- function(
  stage,
  dataMod
) {
  fit2 <- .getCpUnsLocGetProbSmoothActualSecond(dataMod, stage)
  if (.getCpUnsLocGetProbSmoothActualCheck(fit2, dataMod)) {
    .debug("Smoothed") # nolint
    return(
      .getCpUnsLocGetProbSmoothActualResponseSuccess(
        fit = fit2,
        dataMod = dataMod
      )$pred
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
  dataMod$probSmooth - 0.0001
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
    return(.getCpUnsLocGetCpEnsureMeta(
      obj = dataMod,
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      reason = "data_mod_not_available"
    ))
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
    return(.getCpUnsLocConditionOut(
      cp = cpInd,
      locGenerated = FALSE,
      locGeneratedDirect = FALSE,
      locSource = "not_calculated",
      locReason = trimObj$info$reason %||% "trim_returned_non_local_cutpoint"
    ))
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
    return(.getCpUnsLocConditionOut(
      cp = cpInd,
      locGenerated = FALSE,
      locGeneratedDirect = FALSE,
      locSource = "not_calculated",
      locReason = "empty_data_mod_after_trimming"
    ))
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
.getCpUnsLocGetCpTrimFlatLeft <- function(
  dataMod,
  exTblStimNoMin,
  chnlSettings,
  probCol
) {
  info <- list(
    applied = FALSE,
    reason = "no_long_low_flat_left_region"
  )

  if (!is.data.frame(dataMod) || nrow(dataMod) < 4L) {
    info$reason <- "too_few_cells_for_derivative_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  dataMod <- dataMod[order(.getCut(dataMod)), , drop = FALSE]
  x <- .getCut(dataMod)
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

  deriv <- diff(probData$prob) / diff(probData$x)
  deriv[!is.finite(deriv)] <- 0
  deriv <- pmax(0, deriv)
  peakDeriv <- max(deriv, na.rm = TRUE)
  if (!is.finite(peakDeriv) || peakDeriv <= 0) {
    info$reason <- "no_positive_probability_derivative"
    return(list("dataMod" = dataMod, "info" = info))
  }

  derivFrac <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locFlatDerivFrac",
    0.1
  )
  derivThreshold <- peakDeriv * derivFrac
  incIdx <- which(deriv >= derivThreshold)
  if (length(incIdx) == 0L) {
    info$reason <- "no_derivative_above_flat_threshold"
    return(list("dataMod" = dataMod, "info" = info))
  }

  startX <- probData$x[min(incIdx)]
  lowMask <- is.finite(x) & x < startX
  highMask <- is.finite(x) & x >= startX
  if (sum(lowMask) == 0L || sum(highMask) == 0L) {
    info$reason <- "left_flat_region_too_short"
    return(list("dataMod" = dataMod, "info" = info))
  }

  lowMean <- mean(probVec[lowMask], na.rm = TRUE)
  peakProb <- max(probVec, na.rm = TRUE)
  lowRel <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locLeftLowRel",
    0.25
  )
  lowAbs <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locLeftLowAbs",
    0.15
  )
  lowEnough <- is.finite(lowMean) &&
    (lowMean < peakProb * lowRel || lowMean < lowAbs)

  exprVec <- .getCut(exTblStimNoMin)
  exprVec <- exprVec[is.finite(exprVec)]
  nLow <- sum(exprVec < startX)
  nHigh <- sum(exprVec >= startX)
  cellFrac <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locLeftCellFrac",
    0.5
  )
  manyCells <- nHigh > 0L && nLow > cellFrac * nHigh

  minProbX <- probData$x[which.min(probData$prob)]
  maxProbX <- probData$x[which.max(probData$prob)]
  transitionLength <- abs(maxProbX - minProbX)
  lowLength <- startX - min(x, na.rm = TRUE)
  lengthFrac <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locLeftLengthFrac",
    0.5
  )
  longRegion <- is.finite(transitionLength) &&
    transitionLength > 0 &&
    is.finite(lowLength) &&
    lowLength > lengthFrac * transitionLength

  info$startX <- startX
  info$lowMean <- lowMean
  info$peakProb <- peakProb
  info$lowRel <- lowRel
  info$lowAbs <- lowAbs
  info$nLow <- nLow
  info$nHigh <- nHigh
  info$cellFrac <- cellFrac
  info$manyCells <- manyCells
  info$lowLength <- lowLength
  info$transitionLength <- transitionLength
  info$lengthFrac <- lengthFrac
  info$longRegion <- longRegion
  info$peakDeriv <- peakDeriv
  info$derivThreshold <- derivThreshold

  if (!lowEnough || !(manyCells || longRegion)) {
    return(list("dataMod" = dataMod, "info" = info))
  }

  info$applied <- TRUE
  info$reason <- "dropped_long_low_probability_left_region"
  list("dataMod" = dataMod[highMask, , drop = FALSE], "info" = info)
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
    "probBsEst",
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
  chnl
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
