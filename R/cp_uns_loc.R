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
.getCpUnsLocConditionCheckMaxX <- function(exTblStimNoMin, cpMin) {
  (quantile(.getCut(exTblStimNoMin), 0.9) + 3 * sd(.getCut(exTblStimNoMin))) <=
    cpMin
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

  # Keep the bandwidth object used here so the later antimode density can use
  # exactly half of the same fixed or adaptive bandwidth.
  densTblRaw <- .getCpUnsLocGetDensRawTabulate(
    stimX = densList$stim$x,
    stimY = densList$stim$y,
    unsX = densList$uns$x,
    unsY = densList$uns$y
  )
  attr(densTblRaw, "locDensityBw") <- densList$bw
  densTblRaw
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

  probTblPos <- .getCpUnsLocProbTblFilter(
    densTbl = densTblRaw,
    probTbl = probTbl,
    exVecStim = exVecStimThreshold,
    exVecUns = exVecUnsThreshold,
    stage = stage
  )

  list(
    all = probTbl,
    pos = probTblPos,
    densityBw = attr(densTblRaw, "locDensityBw")
  )
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
  densTbl,
  probTbl,
  exVecStim,
  exVecUns,
  stage
) {
  .debug("Filtering before smoothing") # nolint
  densTblStim <- densTbl |>
    dplyr::filter(stim == "yes")
  densTblUns <- densTbl |>
    dplyr::filter(stim == "no")
  peakStimIdx <- .getPeakMainLeftIdx(densTblStim$dens)
  peakStimX <- densTblStim$xStim[peakStimIdx]
  peakUnsIdx <- .getPeakMainLeftIdx(densTblUns$dens)
  peakUnsX <- densTblUns$xStim[peakUnsIdx]
  peakX <- max(peakStimX, peakUnsX)

  windowWidthStim <- 1 /
    3 *
    abs(diff(quantile(exVecStim[exVecStim < peakStimX], c(0.05, 1))))
  windowWidthUns <- 1 /
    3 *
    abs(diff(quantile(exVecUns[exVecUns < peakUnsX], c(0.05, 1))))
  windowWidth <- max(windowWidthStim, windowWidthUns, na.rm = TRUE)

  probTbl <- probTbl |>
    dplyr::filter(xStim > peakX + windowWidth) # nolint

  if (nrow(probTbl) <= 5L) {
    return(probTbl)
  }

  # Find the threshold index
  valid_idx <- probTbl |>
    dplyr::arrange(xStim) |>
    dplyr::mutate(
      ge_0025 = probStimNorm >= 0.025,
      ge_0075 = probStimNorm >= 0.075,
      n_remaining = rev(seq_len(dplyr::n())),

      # Calculate tail proportions
      prop_0025 = rev(cumsum(rev(ge_0025))) / n_remaining,
      prop_0075 = rev(cumsum(rev(ge_0075))) / n_remaining,

      # Do both conditions match?
      both_met = prop_0025 >= 0.90 & prop_0075 >= 0.25
    ) |>
    # Get the row index of the first time this becomes true
    dplyr::pull(both_met) |>
    which() |>
    which.min()

  # Slice the table from that first valid point all the way to the end
  if (length(valid_idx) > 0) {
    probTbl |>
      dplyr::arrange(xStim) |>
      dplyr::slice(valid_idx:dplyr::n())
  } else {
    probTbl |> dplyr::slice(0) # Empty if nothing matches
  }
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
      msg = "No responding cells"
    ))
  }

  minProbXPos <- min(probTblList$pos$xStim, na.rm = TRUE)

  margin <- .getCpUnsLocGetDataModMargin(
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsNoMin = exTblUnsThreshold
  )

  minMod <- minProbXPos - margin

  # Keep the lower values in dataMod to anchor/clamp the spline at the periphery
  dataMod <- exTblStimThreshold[
    .getCut(exTblStimThreshold) >= minMod,
    ,
    drop = FALSE
  ]

  # Interpolate using the 'all' table so baseline cells get real probabilities
  probVec <- try(
    approx(
      x = probTblList$all$xStim,
      y = probTblList$all$probStimNorm,
      xout = .getCut(dataMod),
      method = "linear",
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
      msg = "No responding cells"
    ))
  }

  dataMod <- dataMod |>
    dplyr::mutate(probSmooth = probVec)

  # get the bins that we'll need for thinning
  binVec <- .getCpUnsLocGetDataModBinVec(
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold
  )
  attr(dataMod, "binVec") <- binVec

  # Attach this for filtering when it comes to estimating
  # the final response proportion
  attr(dataMod, "minProbXPos") <- minProbXPos

  # Retain the exact bandwidth object used for the stimulated and unstimulated
  # density estimates. The later antimode density uses half this bandwidth.
  attr(dataMod, "locDensityBw") <- probTblList$densityBw

  .thinDataMod(dataMod, maxCellsPerBin = 20)
}

#' @keywords internal
.getCpUnsLocGetDataModMargin <- function(
  exTblStimNoMin,
  exTblUnsNoMin
) {
  spanStim <- diff(quantile(
    .getCut(exTblStimNoMin),
    probs = c(0.05, 0.95),
    na.rm = TRUE
  ))
  spanUns <- diff(quantile(
    .getCut(exTblUnsNoMin),
    probs = c(0.05, 0.95),
    na.rm = TRUE
  ))

  max(spanStim, spanUns) * 0.05
}

.getCpUnsLocGetDataModBinVec <- function(
  exTblStimThreshold,
  exTblUnsThreshold
) {
  stimVals <- .getCut(exTblStimThreshold)
  unsVals <- .getCut(exTblUnsThreshold)

  rng <- range(stimVals, unsVals, na.rm = TRUE)

  seq.int(from = rng[1], to = rng[2], length.out = 512L)
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
  densityBw <- attr(dataMod, "locDensityBw")

  if (!.getCpUnsLocGetProbSmoothCheckNCell(dataMod)) {
    .intSaveNm(
      "not_enough_cells_to_smooth",
      NULL,
      .getInd(dataMod),
      stageChnl,
      pathProject
    )
    dataModOut <- .getCpUnsLocGetProbSmoothCheckNCellOut(dataMod)
    attr(dataModOut, "locDensityBw") <- densityBw
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
  attr(dataModOut, "locDensityBw") <- densityBw
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
        devtol.fit = 0.5,
        steptol.fit = 1e-1,
        maxHalf = 5,
        bfgs = list(steptol.bfgs = 1e-1),
        maxit = 1e1
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

  # The lower-margin values were retained to anchor the monotone smoother.
  # Keep them available while identifying the antimode, global, and marginal
  # filtering thresholds. They are removed only when the final response
  # proportion is calculated.
  info$minProbXPos <- attr(dataMod, "minProbXPos")
  info$clampingCellsRetainedForThresholdSelection <- TRUE

  dataMod <- .getCpUnsLocGetCpTrimSubset(
    dataMod = dataMod,
    keep = order(.getCut(dataMod))
  )
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

  # First remove values below the right-most eligible antimode, when one is
  # found. The antimode upper limit uses (alpha, omega, psi) =
  # (2 / 3, 0.15, 0.75) by default.
  antiObj <- .getCpUnsLocGetCpTrimAntimode(
    dataMod = dataMod,
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

  # Apply the global response-probability filter and then the marginal quality
  # filter sequentially. Their default derivative parameters are
  # (0.05, 0.15, 0.05) and (0.5, 0.15, 0.75), respectively.
  flatObj <- .getCpUnsLocGetCpTrimLowProbLeft(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol
  )
  dataMod <- flatObj$dataMod
  info$lowProbLeft <- flatObj$info
  info$flatLeft <- flatObj$info

  info$applied <- isTRUE(antiObj$info$applied) ||
    isTRUE(flatObj$info$applied)
  if (isTRUE(info$applied)) {
    info$reason <- "trimmed_before_threshold"
  }

  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    info$reason <- "all_cells_removed_by_low_probability_left_trim"
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
.getCpUnsLocGetCpTrimRiseFrac <- function(
  chnlSettings,
  settingName,
  default
) {
  riseFrac <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    settingName,
    default
  )
  riseFrac <- suppressWarnings(as.numeric(riseFrac)[1])
  if (!is.finite(riseFrac) || riseFrac <= 0 || riseFrac > 1) {
    return(default)
  }
  riseFrac
}

#' Return the first available setting from a precedence-ordered name vector
#' @keywords internal
.getCpUnsLocGetCpTrimFirstSetting <- function(
  chnlSettings,
  settingNames,
  default
) {
  for (nm in settingNames) {
    value <- .getCpUnsLocGetCpTrimSetting(
      chnlSettings = chnlSettings,
      nm = nm,
      default = NULL
    )
    if (!is.null(value)) {
      return(value)
    }
  }
  default
}

#' Validate a probability or relative-height parameter
#' @keywords internal
.getCpUnsLocGetCpTrimUnitParameter <- function(
  value,
  default,
  allowZero = FALSE
) {
  value <- suppressWarnings(as.numeric(value)[1])
  validLower <- if (isTRUE(allowZero)) value >= 0 else value > 0
  if (!is.finite(value) || !validLower || value > 1) {
    return(default)
  }
  value
}

#' Obtain stage-specific derivative parameters
#'
#' The preferred setting names use the symbols from the appendix: `Alpha`,
#' `Omega`, and `Psi`. Earlier setting names remain supported for backwards
#' compatibility.
#' @keywords internal
.getCpUnsLocGetCpTrimStageDerivParams <- function(
  chnlSettings,
  stage
) {
  stage <- match.arg(stage, c("antimode", "global", "marginal"))
  stageTitle <- switch(
    stage,
    antimode = "Antimode",
    global = "Global",
    marginal = "Marginal"
  )
  defaults <- switch(
    stage,
    antimode = c(alpha = 2 / 3, omega = 0.15, psi = 0.75),
    global = c(alpha = 0.05, omega = 0.15, psi = 0.05),
    marginal = c(alpha = 0.5, omega = 0.15, psi = 0.75)
  )

  alpha <- .getCpUnsLocGetCpTrimFirstSetting(
    chnlSettings = chnlSettings,
    settingNames = c(
      paste0("loc", stageTitle, "DerivAlpha"),
      "locDerivAlpha",
      paste0("loc", stageTitle, "DerivPeakMinRel"),
      "locDerivPeakMinRel"
    ),
    default = defaults[["alpha"]]
  )
  omega <- .getCpUnsLocGetCpTrimFirstSetting(
    chnlSettings = chnlSettings,
    settingNames = c(
      paste0("loc", stageTitle, "DerivOmega"),
      "locDerivOmega",
      paste0("loc", stageTitle, "DerivPeakProbMin"),
      "locDerivPeakProbMin"
    ),
    default = defaults[["omega"]]
  )
  psi <- .getCpUnsLocGetCpTrimFirstSetting(
    chnlSettings = chnlSettings,
    settingNames = c(
      paste0("loc", stageTitle, "DerivPsi"),
      "locDerivPsi",
      paste0("loc", stageTitle, "DerivRiseFrac")
    ),
    default = defaults[["psi"]]
  )

  list(
    alpha = .getCpUnsLocGetCpTrimUnitParameter(
      alpha,
      defaults[["alpha"]]
    ),
    omega = .getCpUnsLocGetCpTrimUnitParameter(
      omega,
      defaults[["omega"]],
      allowZero = TRUE
    ),
    psi = .getCpUnsLocGetCpTrimUnitParameter(
      psi,
      defaults[["psi"]]
    )
  )
}


#' Get the minimum response probability required at a rise threshold
#' @keywords internal
.getCpUnsLocGetCpTrimRiseProbMin <- function(
  chnlSettings,
  settingName,
  default = 0
) {
  default <- suppressWarnings(as.numeric(default)[1])
  if (!is.finite(default) || default < 0 || default > 1) {
    default <- 0
  }

  commonDefault <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locDerivRiseProbMin",
    default
  )
  commonDefault <- suppressWarnings(as.numeric(commonDefault)[1])
  if (
    !is.finite(commonDefault) ||
      commonDefault < 0 ||
      commonDefault > 1
  ) {
    commonDefault <- default
  }

  probMin <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    settingName,
    commonDefault
  )
  probMin <- suppressWarnings(as.numeric(probMin)[1])
  if (!is.finite(probMin) || probMin < 0 || probMin > 1) {
    return(commonDefault)
  }
  probMin
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
.getCpUnsLocGetCpTrimSubset <- function(dataMod, keep) {
  attrNames <- c(
    "chnlCut",
    "ind",
    "indUns",
    "binVec",
    "minProbXPos",
    "locProbDerivTbl",
    "locProbSmoothMethod",
    "locDensityBw"
  )
  attrList <- lapply(attrNames, function(nm) attr(dataMod, nm))
  names(attrList) <- attrNames

  out <- dataMod[keep, , drop = FALSE]
  for (nm in attrNames) {
    if (!is.null(attrList[[nm]])) {
      attr(out, nm) <- attrList[[nm]]
    }
  }
  out
}

#' Find the left-most valid peak of a non-negative derivative
#'
#' This helper is deliberately generic. It operates only on numeric vectors and
#' has no knowledge of StimGate data frames, marker columns, or model objects.
#' @keywords internal
.getCpUnsLocDerivativePeak <- function(
  x,
  prob,
  deriv,
  alpha = 0.75
) {
  info <- list(reason = "no_valid_derivative_peak")

  x <- suppressWarnings(as.numeric(x))
  prob <- suppressWarnings(as.numeric(prob))
  deriv <- suppressWarnings(as.numeric(deriv))

  if (length(x) != length(prob) || length(x) != length(deriv)) {
    info$reason <- "derivative_peak_input_lengths_differ"
    return(list("index" = NA_integer_, "data" = NULL, "info" = info))
  }

  keep <- is.finite(x) & is.finite(prob) & is.finite(deriv)
  x <- x[keep]
  prob <- prob[keep]
  deriv <- pmax(0, deriv[keep])

  if (length(x) < 3L) {
    info$reason <- "too_few_finite_derivative_points"
    return(list("index" = NA_integer_, "data" = NULL, "info" = info))
  }

  ord <- order(x)
  peakData <- data.frame(
    x = x[ord],
    prob = pmin(1, pmax(0, prob[ord])),
    deriv = deriv[ord]
  )

  alpha <- suppressWarnings(as.numeric(alpha)[1])
  if (!is.finite(alpha) || alpha <= 0 || alpha > 1) {
    alpha <- 0.75
  }

  globalMaxDeriv <- max(peakData$deriv, na.rm = TRUE)
  if (!is.finite(globalMaxDeriv) || globalMaxDeriv <= 0) {
    info$reason <- "no_positive_probability_derivative"
    info$globalMaxDeriv <- globalMaxDeriv
    return(list("index" = NA_integer_, "data" = peakData, "info" = info))
  }

  # Treat a flat-topped local maximum as one peak and represent it by the
  # left-most point of the plateau. Endpoints are not local peaks. The existing
  # global-maximum fallback is retained for cases with no internal peak.
  derivRuns <- rle(peakData$deriv)
  runEnd <- cumsum(derivRuns$lengths)
  runStart <- runEnd - derivRuns$lengths + 1L
  runValue <- derivRuns$values

  if (length(runValue) >= 3L) {
    leftValue <- c(Inf, runValue[-length(runValue)])
    rightValue <- c(runValue[-1L], Inf)
    peakRun <- runStart > 1L &
      runEnd < nrow(peakData) &
      runValue > leftValue &
      runValue > rightValue
    peakIdx <- runStart[peakRun]
  } else {
    peakIdx <- integer(0L)
  }

  usedGlobalFallback <- length(peakIdx) == 0L
  if (isTRUE(usedGlobalFallback)) {
    peakIdx <- which.max(peakData$deriv)
  }

  maxPeakDeriv <- max(peakData$deriv[peakIdx], na.rm = TRUE)
  peakIdxEligible <- peakIdx[
    peakData$deriv[peakIdx] >= alpha * maxPeakDeriv
  ]

  info$alpha <- alpha
  info$peakMinRel <- alpha
  info$globalMaxDeriv <- globalMaxDeriv
  info$maxPeakDeriv <- maxPeakDeriv
  info$usedGlobalMaximumFallback <- usedGlobalFallback
  info$peakSummary <- data.frame(
    idx = peakIdx,
    x = peakData$x[peakIdx],
    prob = peakData$prob[peakIdx],
    deriv = peakData$deriv[peakIdx],
    relToMaxPeak = peakData$deriv[peakIdx] / maxPeakDeriv,
    relToGlobal = peakData$deriv[peakIdx] / globalMaxDeriv,
    eligible = peakIdx %in% peakIdxEligible
  )

  if (length(peakIdxEligible) == 0L) {
    info$reason <- "no_derivative_peak_met_relative_height"
    return(list("index" = NA_integer_, "data" = peakData, "info" = info))
  }

  peakIdxSelected <- min(peakIdxEligible)
  info$reason <- "identified_leftmost_valid_derivative_peak"
  info$peakIdx <- peakIdxSelected
  info$peakX <- peakData$x[peakIdxSelected]
  info$peakProb <- peakData$prob[peakIdxSelected]
  info$peakDeriv <- peakData$deriv[peakIdxSelected]

  list(
    "index" = peakIdxSelected,
    "data" = peakData,
    "info" = info
  )
}

#' Find where a derivative first reaches a fraction of its selected peak
#'
#' This is also generic and calls `.getCpUnsLocDerivativePeak()` to select the
#' reference peak before locating the rise threshold.
#' @keywords internal
.getCpUnsLocDerivativeRiseThreshold <- function(
  x,
  prob,
  deriv,
  alpha = 0.75,
  omega = 0.15,
  psi,
  thresholdProbMin = 0
) {
  peakObj <- .getCpUnsLocDerivativePeak(
    x = x,
    prob = prob,
    deriv = deriv,
    alpha = alpha
  )
  info <- peakObj$info

  if (is.null(peakObj$data) || !is.finite(peakObj$index)) {
    return(list(
      "riseThresholdX" = NA_real_,
      "risingFastX" = NA_real_,
      "info" = info
    ))
  }

  omega <- suppressWarnings(as.numeric(omega)[1])
  if (!is.finite(omega) || omega < 0 || omega > 1) {
    info$reason <- "invalid_derivative_minimum_probability"
    return(list(
      "riseThresholdX" = NA_real_,
      "risingFastX" = NA_real_,
      "info" = info
    ))
  }

  psi <- suppressWarnings(as.numeric(psi)[1])
  if (!is.finite(psi) || psi <= 0 || psi > 1) {
    info$reason <- "invalid_derivative_rise_fraction"
    return(list(
      "riseThresholdX" = NA_real_,
      "risingFastX" = NA_real_,
      "info" = info
    ))
  }

  peakIdx <- peakObj$index
  peakHeight <- peakObj$data$deriv[peakIdx]
  riseHeight <- psi * peakHeight

  # First identify the earliest point at or to the right of the selected peak
  # whose response probability reaches omega.
  omegaIdx <- seq.int(peakIdx, nrow(peakObj$data))
  omegaIdx <- omegaIdx[
    is.finite(peakObj$data$prob[omegaIdx]) &
      peakObj$data$prob[omegaIdx] >= omega
  ]

  info$omega <- omega
  info$peakProbMin <- omega
  info$psi <- psi
  info$riseFrac <- psi
  info$riseHeight <- riseHeight

  if (length(omegaIdx) == 0L) {
    info$reason <- "no_point_at_or_after_peak_met_minimum_probability"
    return(list(
      "riseThresholdX" = NA_real_,
      "risingFastX" = NA_real_,
      "info" = info
    ))
  }

  omegaIdx <- min(omegaIdx)
  info$omegaIdx <- omegaIdx
  info$omegaX <- peakObj$data$x[omegaIdx]
  info$omegaProb <- peakObj$data$prob[omegaIdx]

  if (omegaIdx != peakIdx) {
    # When the selected peak itself has insufficient response probability, do
    # not move the threshold left. Use the first later point reaching omega.
    riseIdxCandidate <- omegaIdx
    info$thresholdBasis <- "minimum_probability_right_of_peak"
  } else {
    # Otherwise move left to the first point where the derivative reaches psi
    # times the selected peak height.
    riseIdx <- seq_len(peakIdx)
    riseIdx <- riseIdx[peakObj$data$deriv[riseIdx] >= riseHeight]
    if (length(riseIdx) == 0L) {
      info$reason <- "derivative_never_reached_fraction_of_selected_peak"
      return(list(
        "riseThresholdX" = NA_real_,
        "risingFastX" = NA_real_,
        "info" = info
      ))
    }
    riseIdxCandidate <- min(riseIdx)
    info$thresholdBasis <- "left_rise_fraction_of_selected_peak"
  }

  riseThresholdXCandidate <- peakObj$data$x[riseIdxCandidate]
  riseThresholdProbCandidate <- peakObj$data$prob[riseIdxCandidate]

  thresholdProbMin <- suppressWarnings(as.numeric(thresholdProbMin)[1])
  if (
    !is.finite(thresholdProbMin) ||
      thresholdProbMin < 0 ||
      thresholdProbMin > 1
  ) {
    thresholdProbMin <- 0
  }

  info$riseThresholdIdxCandidate <- riseIdxCandidate
  info$riseThresholdXCandidate <- riseThresholdXCandidate
  info$riseThresholdProbCandidate <- riseThresholdProbCandidate
  info$thresholdProbMin <- thresholdProbMin

  # This additional, optional constraint applies to the final threshold rather
  # than to peak selection. If needed, shift the threshold rightwards until the
  # fitted response probability reaches the requested minimum.
  laterIdx <- seq.int(riseIdxCandidate, nrow(peakObj$data))
  laterIdx <- laterIdx[
    is.finite(peakObj$data$prob[laterIdx]) &
      peakObj$data$prob[laterIdx] >= thresholdProbMin
  ]

  if (length(laterIdx) == 0L) {
    info$reason <- "no_later_rise_threshold_met_minimum_probability"
    return(list(
      "riseThresholdX" = NA_real_,
      "risingFastX" = NA_real_,
      "info" = info
    ))
  }

  riseIdx <- min(laterIdx)
  riseThresholdX <- peakObj$data$x[riseIdx]
  riseThresholdProb <- peakObj$data$prob[riseIdx]

  info$reason <- "identified_probability_rise_threshold"
  info$riseThresholdIdx <- riseIdx
  info$riseThresholdX <- riseThresholdX
  info$riseThresholdProb <- riseThresholdProb
  info$riseThreshold <- riseHeight
  info$shiftedRightForProbability <- riseIdx > riseIdxCandidate
  info$risingFastIdx <- riseIdx
  info$risingFastX <- riseThresholdX

  list(
    "riseThresholdX" = riseThresholdX,
    "risingFastX" = riseThresholdX,
    "info" = info
  )
}

#' Obtain a probability-rise threshold from a StimGate model object
#' @keywords internal
.getCpUnsLocGetCpTrimRiseThreshold <- function(
  dataMod,
  chnlSettings,
  probCol,
  alpha = NULL,
  omega = NULL,
  psi = NULL,
  riseFrac = NULL,
  thresholdProbMin = 0
) {
  info <- list(reason = "no_valid_probability_derivative_peak")

  derivTbl <- .getCpUnsLocGetCpTrimDerivativeTbl(
    dataMod = dataMod,
    probCol = probCol
  )
  if (is.null(derivTbl) || nrow(derivTbl) < 3L) {
    info$reason <- "probability_derivative_unavailable"
    return(list(
      "riseThresholdX" = NA_real_,
      "risingFastX" = NA_real_,
      "info" = info
    ))
  }

  # `riseFrac` is retained as a backwards-compatible alias for `psi`.
  if (is.null(psi)) {
    psi <- riseFrac
  }
  if (is.null(alpha)) {
    alpha <- .getCpUnsLocGetCpTrimSetting(
      chnlSettings,
      "locDerivAlpha",
      .getCpUnsLocGetCpTrimSetting(
        chnlSettings,
        "locDerivPeakMinRel",
        0.75
      )
    )
  }
  if (is.null(omega)) {
    omega <- .getCpUnsLocGetCpTrimSetting(
      chnlSettings,
      "locDerivOmega",
      .getCpUnsLocGetCpTrimSetting(
        chnlSettings,
        "locDerivPeakProbMin",
        0.15
      )
    )
  }
  if (is.null(psi)) {
    psi <- .getCpUnsLocGetCpTrimSetting(
      chnlSettings,
      "locDerivPsi",
      0.75
    )
  }

  riseObj <- .getCpUnsLocDerivativeRiseThreshold(
    x = derivTbl$x,
    prob = derivTbl$prob,
    deriv = derivTbl$deriv,
    alpha = alpha,
    omega = omega,
    psi = psi,
    thresholdProbMin = thresholdProbMin
  )
  riseObj$info$derivSource <- unique(derivTbl$source)[1]
  riseObj
}

# Antimode-specific derivative threshold using appendix defaults.
#' @keywords internal
.getCpUnsLocGetCpTrimAntimodeRise <- function(
  dataMod,
  chnlSettings,
  probCol
) {
  params <- .getCpUnsLocGetCpTrimStageDerivParams(
    chnlSettings = chnlSettings,
    stage = "antimode"
  )
  thresholdProbMin <- .getCpUnsLocGetCpTrimRiseProbMin(
    chnlSettings = chnlSettings,
    settingName = "locAntimodeDerivRiseProbMin",
    default = 0
  )

  .getCpUnsLocGetCpTrimRiseThreshold(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol,
    alpha = params$alpha,
    omega = params$omega,
    psi = params$psi,
    thresholdProbMin = thresholdProbMin
  )
}

#' @keywords internal
.getCpUnsLocGetCpTrimAntimode <- function(
  dataMod,
  chnlSettings,
  probCol,
  riseObj = NULL
) {
  info <- list(
    applied = FALSE,
    reason = "antimode_trim_not_applied"
  )

  if (!is.data.frame(dataMod) || nrow(dataMod) < 5L) {
    info$reason <- "too_few_model_values_for_antimode_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  if (is.null(riseObj)) {
    riseObj <- .getCpUnsLocGetCpTrimAntimodeRise(
      dataMod = dataMod,
      chnlSettings = chnlSettings,
      probCol = probCol
    )
  }
  info$rise <- riseObj$info

  riseThresholdX <- riseObj$riseThresholdX %||% riseObj$risingFastX
  if (!is.finite(riseThresholdX)) {
    info$reason <- riseObj$info$reason %||%
      "no_valid_probability_derivative_peak"
    return(list("dataMod" = dataMod, "info" = info))
  }

  # Fit the antimode density only to values in the dubious-response region,
  # including the derivative-based upper endpoint.
  exprVec <- suppressWarnings(as.numeric(.getCut(dataMod)))
  exprVec <- exprVec[
    is.finite(exprVec) &
      exprVec <= riseThresholdX
  ]
  if (length(exprVec) < 5L || length(unique(exprVec)) < 3L) {
    info$reason <- "too_few_expression_values_below_antimode_upper_limit"
    return(list("dataMod" = dataMod, "info" = info))
  }

  densObj <- .getCpUnsLocGetCpTrimDensity(
    exprVec = exprVec,
    chnlSettings = chnlSettings,
    densityBw = attr(dataMod, "locDensityBw")
  )
  if (is.null(densObj)) {
    info$reason <- "antimode_density_failed"
    return(list("dataMod" = dataMod, "info" = info))
  }

  antimodeX <- .getCpUnsLocGetCpTrimAntimodes(densObj)
  antimodeLeftX <- antimodeX[
    is.finite(antimodeX) & antimodeX <= riseThresholdX
  ]

  info$densityBwType <- attr(densObj, "locBwType")
  info$densityBwFraction <- attr(densObj, "locBwFraction")
  info$densityBwBase <- attr(densObj, "locBwBaseSummary")
  info$densityBwUsed <- attr(densObj, "locBwUsedSummary")
  info$antimodeX <- antimodeX
  info$antimodeLeftX <- antimodeLeftX
  info$riseThresholdX <- riseThresholdX
  info$risingFastX <- riseThresholdX
  info$nExpressionValuesForDensity <- length(exprVec)

  if (length(antimodeLeftX) == 0L) {
    info$reason <- "no_antimode_in_dubious_response_region"
    return(list("dataMod" = dataMod, "info" = info))
  }

  # Any antimode is accepted. Use the right-most one in the dubious-response
  # region as the filtering point.
  filterX <- max(antimodeLeftX, na.rm = TRUE)
  x <- suppressWarnings(as.numeric(.getCut(dataMod)))
  keep <- is.finite(x) & x >= filterX

  info$filterX <- filterX
  info$nDropped <- sum(!keep, na.rm = TRUE)
  info$applied <- info$nDropped > 0L
  info$reason <- if (isTRUE(info$applied)) {
    "dropped_values_left_of_rightmost_eligible_antimode"
  } else {
    "rightmost_eligible_antimode_kept_all_values"
  }

  list(
    "dataMod" = .getCpUnsLocGetCpTrimSubset(dataMod, keep),
    "info" = info
  )
}

#' @keywords internal
.getCpUnsLocGetCpTrimDerivativeTbl <- function(
  dataMod,
  probCol
) {
  xRange <- range(
    suppressWarnings(as.numeric(.getCut(dataMod))),
    na.rm = TRUE
  )
  if (length(xRange) != 2L || any(!is.finite(xRange)) || diff(xRange) <= 0) {
    return(NULL)
  }

  # Preferred path: the finite-difference derivative evaluated across the fitted
  # monotone response-probability curve when the smoother was fitted.
  derivTbl <- attr(dataMod, "locProbDerivTbl")
  if (
    identical(probCol, "pred") &&
      is.data.frame(derivTbl) &&
      all(c("x", "pred", "deriv") %in% names(derivTbl))
  ) {
    derivTbl <- derivTbl |>
      dplyr::transmute(
        x = suppressWarnings(as.numeric(.data$x)),
        prob = pmin(1, pmax(0, suppressWarnings(as.numeric(.data$pred)))),
        deriv = pmax(0, suppressWarnings(as.numeric(.data$deriv))),
        source = "fitted_probability_finite_difference"
      ) |>
      dplyr::filter(
        is.finite(.data$x),
        is.finite(.data$prob),
        is.finite(.data$deriv),
        .data$x >= .env$xRange[1],
        .data$x <= .env$xRange[2]
      ) |>
      dplyr::arrange(.data$x)

    if (nrow(derivTbl) >= 3L) {
      return(derivTbl)
    }
  }

  # Fallback for an unavailable fitted derivative: finite differences of the
  # selected probability column over the filtered model values.
  x <- suppressWarnings(as.numeric(.getCut(dataMod)))
  prob <- .getCpUnsLocGetCpTrimProbVec(dataMod, probCol)
  probTbl <- tibble::tibble(x = x, prob = prob) |>
    dplyr::filter(is.finite(.data$x), is.finite(.data$prob)) |>
    dplyr::group_by(.data$x) |>
    dplyr::summarise(prob = mean(.data$prob), .groups = "drop") |>
    dplyr::arrange(.data$x)

  if (nrow(probTbl) < 4L || diff(range(probTbl$x)) <= 0) {
    return(NULL)
  }

  dx <- diff(probTbl$x)
  deriv <- diff(probTbl$prob) / dx
  deriv[!is.finite(deriv)] <- 0
  deriv <- pmax(0, deriv)

  tibble::tibble(
    x = (probTbl$x[-1L] + probTbl$x[-nrow(probTbl)]) / 2,
    prob = (probTbl$prob[-1L] + probTbl$prob[-nrow(probTbl)]) / 2,
    deriv = deriv,
    source = "model_probability_finite_difference"
  )
}

# Retain the previous name for compatibility with direct internal calls.
#' @keywords internal
.getCpUnsLocGetCpTrimAntimodeDerivativeTbl <- function(
  dataMod,
  probCol
) {
  .getCpUnsLocGetCpTrimDerivativeTbl(
    dataMod = dataMod,
    probCol = probCol
  )
}

#' @keywords internal
.getCpUnsLocGetCpTrimDensity <- function(
  exprVec,
  chnlSettings,
  densityBw = NULL
) {
  bwFraction <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locAntimodeBwFrac",
    1 / 2
  )
  bwFraction <- suppressWarnings(as.numeric(bwFraction[1]))
  if (!is.finite(bwFraction) || bwFraction <= 0) {
    bwFraction <- 1 / 2
  }

  exprRange <- range(exprVec, na.rm = TRUE)

  # If the original density used an adaptive bandwidth curve, halve that curve
  # pointwise and estimate the antimode density on the same grid.
  if (
    is.list(densityBw) &&
      isTRUE(densityBw$adaptive) &&
      !is.null(densityBw$grid) &&
      !is.null(densityBw$sharedGrid)
  ) {
    grid <- suppressWarnings(as.numeric(densityBw$grid))
    bwBase <- suppressWarnings(as.numeric(densityBw$sharedGrid))
    keepGrid <- is.finite(grid) &
      is.finite(bwBase) &
      bwBase > 0 &
      grid >= exprRange[1] &
      grid <= exprRange[2]
    grid <- grid[keepGrid]
    bwBase <- bwBase[keepGrid]

    if (length(grid) >= 3L && length(grid) == length(bwBase)) {
      densObj <- .getCpUnsLocDensityAdaptiveGrid(
        x = exprVec,
        grid = grid,
        bwGrid = bwBase * bwFraction,
        normalise = TRUE
      )
      if (!is.null(densObj)) {
        attr(densObj, "locBwType") <- "adaptive"
        attr(densObj, "locBwFraction") <- bwFraction
        attr(densObj, "locBwBaseSummary") <- .getCpUnsLocGetCpTrimBwSummary(
          bwBase
        )
        attr(densObj, "locBwUsedSummary") <- .getCpUnsLocGetCpTrimBwSummary(
          bwBase * bwFraction
        )
        return(densObj)
      }
    }
  }

  bwBase <- .getCpUnsLocGetCpTrimDensityBw(
    exprVec = exprVec,
    chnlSettings = chnlSettings,
    densityBw = densityBw
  )
  if (!is.finite(bwBase) || bwBase <= 0) {
    return(NULL)
  }

  bwUsed <- bwBase * bwFraction
  densObj <- try(
    suppressWarnings(stats::density(
      exprVec,
      bw = bwUsed,
      n = 512L,
      from = exprRange[1],
      to = exprRange[2]
    )),
    silent = TRUE
  )
  if (inherits(densObj, "try-error")) {
    return(NULL)
  }

  attr(densObj, "locBwType") <- "fixed"
  attr(densObj, "locBwFraction") <- bwFraction
  attr(densObj, "locBwBaseSummary") <- .getCpUnsLocGetCpTrimBwSummary(bwBase)
  attr(densObj, "locBwUsedSummary") <- .getCpUnsLocGetCpTrimBwSummary(bwUsed)
  densObj
}

#' @keywords internal
.getCpUnsLocGetCpTrimDensityBw <- function(
  exprVec,
  chnlSettings,
  densityBw = NULL
) {
  bwBase <- if (is.list(densityBw)) {
    NA_real_
  } else {
    .getCpUnsLocGetCpTrimAsFiniteNumeric(
      densityBw,
      positive = TRUE
    )
  }
  if (is.finite(bwBase)) {
    return(bwBase)
  }

  # This fallback is mainly for direct calls that bypass the normal pipeline.
  # In ordinary use, locDensityBw is attached when the first density is fitted.
  settingNames <- c("bw", "bwCluster", "bwFallback")
  for (nm in settingNames) {
    bwBase <- .getCpUnsLocGetCpTrimAsFiniteNumeric(
      .getCpUnsLocGetCpTrimSetting(chnlSettings, nm, NA_real_),
      positive = TRUE
    )
    if (is.finite(bwBase)) {
      return(bwBase)
    }
  }

  hpiBw <- try(
    suppressWarnings(ks::hpi(x = exprVec)),
    silent = TRUE
  )
  if (inherits(hpiBw, "try-error")) {
    return(NA_real_)
  }
  .getCpUnsLocGetCpTrimAsFiniteNumeric(hpiBw, positive = TRUE)
}

#' @keywords internal
.getCpUnsLocGetCpTrimBwSummary <- function(bw) {
  bw <- suppressWarnings(as.numeric(bw))
  bw <- bw[is.finite(bw) & bw > 0]
  if (length(bw) == 0L) {
    return(c(min = NA_real_, median = NA_real_, max = NA_real_))
  }
  c(
    min = min(bw),
    median = stats::median(bw),
    max = max(bw)
  )
}

#' @keywords internal
.getCpUnsLocGetCpTrimAntimodes <- function(densObj) {
  x <- suppressWarnings(as.numeric(densObj$x))
  y <- suppressWarnings(as.numeric(densObj$y))
  if (
    length(x) != length(y) ||
      length(y) < 3L ||
      all(!is.finite(y))
  ) {
    return(numeric(0L))
  }

  y[!is.finite(y)] <- Inf
  minIdx <- .getLocalMinimaIdx(y)
  sort(unique(x[minIdx]))
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

  dataMod <- .getCpUnsLocGetCpTrimSubset(
    dataMod = dataMod,
    keep = order(.getCut(dataMod))
  )
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

  gridObj <- .getCpUnsLocGetCpTrimMarginalBreaks(
    dataMod = dataMod,
    startX = startX
  )
  if (is.null(gridObj)) {
    info$reason <- "insufficient_bins_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  breaks <- gridObj$breaks
  refIndex <- gridObj$refIndex
  refNBin <- length(breaks) - refIndex
  if (refNBin < 1L) {
    info$reason <- "no_right_reference_bins_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  rightMask <- dm$x >= startX
  refNCell <- sum(rightMask, na.rm = TRUE)
  refExpectedResp <- sum(dm$prob[rightMask], na.rm = TRUE)
  if (refNCell == 0L || !is.finite(refExpectedResp)) {
    info$reason <- "no_right_reference_cells_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  refCellsPerBin <- refNCell / refNBin
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

  leftBinIndex <- seq_len(refIndex - 1L)
  if (length(leftBinIndex) == 0L) {
    info$reason <- "no_left_bins_for_marginal_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }
  leftBinIndex <- rev(leftBinIndex)

  currentCut <- startX
  scanList <- vector("list", length(leftBinIndex))
  stopReason <- "accepted_all_left_bins"

  for (j in seq_along(leftBinIndex)) {
    i <- leftBinIndex[[j]]
    left <- breaks[[i]]
    right <- breaks[[i + 1L]]
    mask <- dm$x >= left & dm$x < right
    nCell <- sum(mask, na.rm = TRUE)
    expectedResp <- sum(dm$prob[mask], na.rm = TRUE)
    purity <- if (nCell > 0L) expectedResp / nCell else NA_real_

    acceptBin <- nCell == 0L ||
      (is.finite(purity) &&
        nCell <= maxLeftBinCells &&
        purity >= minLeftBinPurity)

    scanList[[j]] <- tibble::tibble(
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
      relCellsPerBin = ifelse(
        refCellsPerBin > 0,
        nCell / refCellsPerBin,
        Inf
      ),
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
  info$gridShift <- gridObj$shift
  info$gridSpacing <- gridObj$spacing
  info$refIndex <- refIndex
  info$refNBin <- refNBin
  info$refNCell <- refNCell
  info$refExpectedResp <- refExpectedResp
  info$refCellsPerBin <- refCellsPerBin
  info$refPurity <- refPurity
  info$cellBinRatio <- cellBinRatio
  info$purityRel <- purityRel
  info$maxLeftBinCells <- maxLeftBinCells
  info$minLeftBinPurity <- minLeftBinPurity
  info$scanTbl <- scanTbl

  list(
    "dataMod" = .getCpUnsLocGetCpTrimSubset(dataMod, keep),
    "info" = info
  )
}

#' @keywords internal
.getCpUnsLocGetCpTrimMarginalBreaks <- function(
  dataMod,
  startX
) {
  binVec <- attr(dataMod, "binVec")
  if (is.null(binVec)) {
    x <- suppressWarnings(as.numeric(.getCut(dataMod)))
    x <- x[is.finite(x)]
    if (length(x) < 2L || diff(range(x)) <= 0) {
      return(NULL)
    }
    binVec <- seq(
      min(x),
      max(x),
      length.out = 512L
    )
  }

  binVec <- suppressWarnings(as.numeric(binVec))
  binVec <- sort(unique(binVec[is.finite(binVec)]))
  if (length(binVec) < 2L || !is.finite(startX)) {
    return(NULL)
  }

  spacing <- stats::median(diff(binVec), na.rm = TRUE)
  if (!is.finite(spacing) || spacing <= 0) {
    return(NULL)
  }

  tol <- max(
    sqrt(.Machine$double.eps) * max(abs(c(binVec, startX)), 1),
    spacing * 1e-10
  )
  eligibleAnchor <- which(binVec <= startX + tol)
  if (length(eligibleAnchor) == 0L) {
    return(NULL)
  }

  # The largest original grid point not exceeding x_ref gives the smallest
  # non-negative shift that makes a shifted grid point exactly equal x_ref.
  refIndex <- max(eligibleAnchor)
  shift <- startX - binVec[[refIndex]]
  if (!is.finite(shift) || shift < -tol) {
    return(NULL)
  }
  shift <- max(0, shift)

  breaks <- binVec + shift
  breaks[[refIndex]] <- startX
  if (any(!is.finite(breaks)) || any(diff(breaks) <= 0)) {
    return(NULL)
  }

  list(
    breaks = breaks,
    refIndex = refIndex,
    shift = shift,
    spacing = spacing
  )
}

#' Apply global and marginal filtering using separate rise thresholds
#' @keywords internal
.getCpUnsLocGetCpTrimLowProbLeft <- function(
  dataMod,
  chnlSettings,
  probCol,
  riseObj = NULL
) {
  info <- list(
    applied = FALSE,
    reason = "low_probability_left_trim_not_run"
  )

  if (!is.data.frame(dataMod) || nrow(dataMod) < 4L) {
    info$reason <- "too_few_cells_for_low_probability_left_trim"
    return(list("dataMod" = dataMod, "info" = info))
  }

  # Retained only for compatibility with the previous internal interface.
  force(riseObj)

  globalParams <- .getCpUnsLocGetCpTrimStageDerivParams(
    chnlSettings = chnlSettings,
    stage = "global"
  )
  globalThresholdProbMin <- .getCpUnsLocGetCpTrimRiseProbMin(
    chnlSettings = chnlSettings,
    settingName = "locGlobalDerivRiseProbMin",
    default = 0
  )
  globalRiseObj <- .getCpUnsLocGetCpTrimRiseThreshold(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol,
    alpha = globalParams$alpha,
    omega = globalParams$omega,
    psi = globalParams$psi,
    thresholdProbMin = globalThresholdProbMin
  )
  info$globalRise <- globalRiseObj$info
  info$globalParams <- globalParams
  info$globalThresholdProbMin <- globalThresholdProbMin

  globalStartX <- globalRiseObj$riseThresholdX %||%
    globalRiseObj$risingFastX
  x <- suppressWarnings(as.numeric(.getCut(dataMod)))

  if (is.finite(globalStartX)) {
    globalKeep <- is.finite(x) & x >= globalStartX
    nGlobalDropped <- sum(!globalKeep, na.rm = TRUE)
    dataModGlobal <- .getCpUnsLocGetCpTrimSubset(dataMod, globalKeep)
  } else {
    nGlobalDropped <- 0L
    dataModGlobal <- dataMod
  }

  info$globalStartX <- globalStartX
  info$nGlobalDropped <- nGlobalDropped
  info$globalDefined <- is.finite(globalStartX)
  # Compatibility with earlier field names.
  info$globalRiseFrac <- globalParams$psi
  info$hardStartX <- globalStartX
  info$hardDerivFrac <- globalParams$psi
  info$hardDerivThreshold <- globalRiseObj$info$riseHeight
  info$peakDeriv <- globalRiseObj$info$peakDeriv
  info$derivSource <- globalRiseObj$info$derivSource

  if (!is.data.frame(dataModGlobal) || nrow(dataModGlobal) == 0L) {
    info$applied <- nGlobalDropped > 0L
    info$reason <- "all_cells_removed_by_global_derivative_trim"
    return(list("dataMod" = dataModGlobal, "info" = info))
  }

  # Calculate the marginal reference threshold only after applying a defined
  # global threshold. An undefined global threshold does not block this stage.
  marginalParams <- .getCpUnsLocGetCpTrimStageDerivParams(
    chnlSettings = chnlSettings,
    stage = "marginal"
  )
  marginalThresholdProbMin <- .getCpUnsLocGetCpTrimRiseProbMin(
    chnlSettings = chnlSettings,
    settingName = "locMarginalDerivRiseProbMin",
    default = 0
  )
  marginalRiseObj <- .getCpUnsLocGetCpTrimRiseThreshold(
    dataMod = dataModGlobal,
    chnlSettings = chnlSettings,
    probCol = probCol,
    alpha = marginalParams$alpha,
    omega = marginalParams$omega,
    psi = marginalParams$psi,
    thresholdProbMin = marginalThresholdProbMin
  )
  info$marginalRise <- marginalRiseObj$info
  info$marginalParams <- marginalParams
  info$marginalThresholdProbMin <- marginalThresholdProbMin

  marginalStartX <- marginalRiseObj$riseThresholdX %||%
    marginalRiseObj$risingFastX
  info$marginalStartX <- marginalStartX
  info$marginalDefined <- is.finite(marginalStartX)
  # Compatibility with earlier field names.
  info$marginalRiseFrac <- marginalParams$psi
  info$startX <- marginalStartX
  info$derivFrac <- marginalParams$psi
  info$derivThreshold <- marginalRiseObj$info$riseHeight

  if (!is.finite(marginalStartX)) {
    info$applied <- nGlobalDropped > 0L
    info$reason <- if (isTRUE(info$applied)) {
      "applied_global_trim_but_no_valid_marginal_rise_threshold"
    } else {
      "no_valid_global_or_marginal_probability_rise_threshold"
    }
    return(list("dataMod" = dataModGlobal, "info" = info))
  }

  marginalFilterObj <- .getCpUnsLocGetCpTrimMarginalLeft(
    dataMod = dataModGlobal,
    chnlSettings = chnlSettings,
    probCol = probCol,
    startX = marginalStartX
  )

  globalApplied <- nGlobalDropped > 0L
  marginalApplied <- isTRUE(marginalFilterObj$info$applied)

  info$applied <- globalApplied || marginalApplied
  info$reason <- if (isTRUE(info$applied)) {
    "dropped_left_region_by_global_or_marginal_trim"
  } else {
    "global_and_marginal_trim_kept_all_cells"
  }
  info$marginal <- marginalFilterObj$info

  list("dataMod" = marginalFilterObj$dataMod, "info" = info)
}

# Retain the old name for callers outside this file.
#' @keywords internal
.getCpUnsLocGetCpTrimFlatLeft <- function(
  dataMod,
  exTblStimNoMin = NULL,
  chnlSettings,
  probCol,
  riseObj = NULL
) {
  force(exTblStimNoMin)
  .getCpUnsLocGetCpTrimLowProbLeft(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol,
    riseObj = riseObj
  )
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
  # Remove the lower-margin values retained only to anchor the smoother at the
  # point where the final response proportion is calculated, not while the
  # filtering thresholds are identified.
  dataModEstimate <- .getCpUnsLocGetCpDataThresholdExcludeMargin(dataMod)

  dataCount <- .getCpUnsLocGetCpDataThresholdCount(dataModEstimate)
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

#' Exclude values retained only as a lower smoothing margin
#' @keywords internal
.getCpUnsLocGetCpDataThresholdExcludeMargin <- function(dataMod) {
  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    return(dataMod)
  }

  minProbXPos <- attr(dataMod, "minProbXPos")
  minProbXPos <- suppressWarnings(as.numeric(minProbXPos)[1])
  if (!is.finite(minProbXPos)) {
    return(dataMod)
  }

  x <- suppressWarnings(as.numeric(.getCut(dataMod)))
  .getCpUnsLocGetCpTrimSubset(
    dataMod = dataMod,
    keep = is.finite(x) & x >= minProbXPos
  )
}

#' @keywords internal
.getCpUnsLocGetCpDataThresholdCount <- function(dataMod) {
  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    return(dataMod)
  }

  if (nrow(dataMod) == 1L) {
    minVal <- min(.getCut(dataMod)) - 1
  } else {
    minVal <- min(.getCut(dataMod))
  }
  dataMod <- dataMod[.getCut(dataMod) > minVal, , drop = FALSE]
  dataMod <- dataMod[order(.getCut(dataMod)), , drop = FALSE]
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
