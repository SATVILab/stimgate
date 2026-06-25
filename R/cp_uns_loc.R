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
    ind = names(exList)[length(exList)],
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
    exListNoMinStim = exListNoMin[-length(exListNoMin)],
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
  indUns <- names(exList)[length(exList)]
  exTblStim <- exList[seq_len(length(exList) - 1)] |>
    dplyr::bind_rows()
  exTblStim <- exTblStim[order(.getCut(exTblStim)), ]
  list(
    exTblStim,
    exList[[length(exList)]]
  ) |>
    stats::setNames(c(names(exList)[1], indUns))
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
  .combineCp(
    cp = cpUnsListNonjoin[["loc"]],
    gateCombn = nonPrejoinCombnVec
  )
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
      exTblUnsOrig <- exListOrig[[length(exListOrig)]]
      exTblStimOrig <- exListOrig[[i]]
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

      .getCpUnsLocInd(
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
    indUns = names(exListOrig)[length(exListOrig)],
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
  objOut <- .getCpUnsLocIndCheckOut(
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
.getCpUnsLocInd <- function(
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
    objOut <- .getCpUnsLocIndTooFew(
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
    .getCpUnsLocIndMaxDensX(exTblStimNoMin)
  )
  exTblUnsThreshold <- .getCpUnsLocSetMaxExpr(
    exTblUnsBias,
    .getCpUnsLocIndMaxDensX(exTblStimNoMin)
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
    pathProject = pathProject
  )
}

.getCpUnsLocIndTooFew <- function(
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
  objOut <- .getCpUnsLocIndCheckOut(
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
.getCpUnsLocIndCheckNCell <- function(exTblStimNoMin, minCell) {
  nrow(exTblStimNoMin) < minCell
}

#' @keywords internal
.getCpUnsLocIndMaxDensX <- function(exTblStimNoMin) {
  max(.getCut(exTblStimNoMin)) -
    0.05 * (diff(range(.getCut(exTblStimNoMin))))
}

#' @keywords internal
.getCpUnsLocIndCheckMaxX <- function(exTblStimNoMin, cpMin) {
  .getCpUnsLocIndMaxDensX(exTblStimNoMin) <= cpMin
}

#' @keywords internal
.getCpUnsLocCheckEarly <- function(exTblStimNoMin, minCell, cpMin) {
  .getCpUnsLocIndCheckNCell(exTblStimNoMin, minCell) ||
    .getCpUnsLocIndCheckMaxX(exTblStimNoMin, cpMin)
}

#' @keywords internal
.getCpUnsLocIndCheckOut <- function(
  cpMin,
  exTblStimNoMin,
  exTblUnsBias,
  stage,
  msg
) {
  .debug(msg) # nolint
  list(
    cp = .getCpUnsLocIndCpNonLoc(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias
    ),
    pList = .getCpUnsLocPListEmpty()
  )
}

#' @keywords internal
.getCpUnsLocIndCpNonLoc <- function(
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
    sdX <- abs(.data) / 1.5
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

  probTbl |>
    dplyr::filter(
      cumsum(probStimNorm >= 0.025) > 0
    )

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
    max(.getCut(exTblStimOrig)) <
      .getCpUnsLocGetMinProbX(probTblPos)
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
    return(.getCpUnsLocIndCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "No responding cells" # nolint
    ))
  }
  margin <- getCpUnsLocGetDataModMargin(
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsNoMin = exTblUnsThreshold
  )

  dataMod <- exTblStimThreshold
  dataMod <- dataMod[
    .getCut(dataMod) >=
      (min(.getCpUnsLocGetMinProbX(probTblList$pos) - margin)),
  ]
  if (nrow(dataMod) == 0L) {
    return(.getCpUnsLocIndCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "No responding cells" # nolint
    ))
  }
  probVec <- approx(
    x = probTblList$pos$xStim,
    y = probTblList$pos$probStimNorm,
    xout = dataMod[[1]],
    method = "linear",
    f = 0.5,
    rule = 2
  )$y

  dataMod |>
    dplyr::mutate(probSmooth = probVec)
}

getCpUnsLocGetDataModMargin <- function(
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

#' @keywords internal
.getCpUnsLocGetProbSmoothActualFirst <- function(dataMod, stage) {
  .debug("Smoothing I") # nolint
  try(
    {
      fml <- as.formula(paste0(
        "probSmooth ~ s(xStim), bs = 'mpi')"
      ))
      scam::scam(
        fml, # nolint
        family = "binomial",
        .data = dataMod |>
          dplyr::mutate(
            probSmooth = pmin(probSmooth, 0.999),
            probSmooth = pmax(probSmooth, 0.001)
          ),
        control = scam::scam.control(
          print.warn = FALSE,
          trace = FALSE,
          devtol.fit = 0.5,
          steptol.fit = 1e-1,
          maxHalf = 5,
          bfgs = list(steptol.bfgs = 1e-1),
          maxit = 1e1
        )
      )
    },
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
  predVec <- predict(fit, type = "response")
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
  try(
    {
      fml <- as.formula(paste0(
        "probSmooth ~ s(xStim), bs = 'micv')"
      ))
      scam::scam(
        fml,
        family = "binomial",
        .data = dataMod,
        control = scam::scam.control(
          print.warn = FALSE,
          trace = FALSE,
          devtol.fit = 0.01
        )
      )
    },
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
  pathProject
) {
  ind <- .getInd(exTblStimNoMin)
  chnl <- .getCpUnsLocGetChnl(exTblStimNoMin)
  stageChnl <- file.path(stage, chnl)
  if (!is.data.frame(dataMod)) {
    .intSaveNm("noDataModDf", NULL, ind, stageChnl, pathProject)
    .intSaveNm("cpInd", dataMod, ind, stageChnl, pathProject)
    return(dataMod)
  }

  dataThreshold <- .getCpUnsLocGetCpDataThreshold(
    dataMod = dataMod,
    exTblStimOrig = exTblStimOrig,
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsBias = exTblUnsBias,
    exTblUnsOrig = exTblUnsOrig,
    bias = bias
  )
  .intSave(ind, stageChnl, pathProject, dataThreshold)
  cpInd <- .getCpUnsLocGetCpActual(
    dataThreshold = dataThreshold,
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsBias = exTblUnsBias,
    cpMin = cpMin,
    stage = stage
  )
  .intSave(ind, stageChnl, pathProject, cpInd)
  .debug("Completed loc gate for single sample") # nolint
  list("cp" = cpInd, "pList" = .getCpUnsLocPListEmpty())
}

#' @keywords internal
.getCpUnsLocGetCpDataThreshold <- function(
  dataMod,
  exTblStimOrig,
  exTblStimNoMin,
  exTblUnsBias,
  exTblUnsOrig,
  bias
) {
  dataCount <- .getCpUnsLocGetCpDataThresholdCount(dataMod)
  probBsEst <- .getCpUnsLocGetCpDataThresholdPropBsEst(
    dataCount = dataCount,
    exTblStimOrig = exTblStimOrig
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
      propUns = propUnsVec,
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
    return(.getCpUnsLocIndCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "Too few responding cells"
    )[["cp"]])
  }
  dataThreshold <- dataThreshold |>
    dplyr::filter(abs(propBsDiff) == min(abs(propBsDiff))) |> # nolint
    dplyr::slice(1) |>
    .getCut()
  dataThreshold
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
  .intSaveNm(
    "cpVecBeforeRep",
    cpVec,
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
  } else {
    .intSaveNm(
      "individualCpUsed",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpVec <- stats::setNames(cpVec, indStim)
  }
  .intSaveNm(
    "cpVecAfterRep",
    cpVec,
    indCombined,
    stageChnl,
    pathProject
  )

  if (!all(purrr::map_lgl(cpVec, is.na))) {
    .intSaveNm(
      "addingUnsCp",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpVec <- c(cpVec, stats::setNames(mean(cpVec, na.rm = TRUE), indUns))
  } else {
    .intSaveNm(
      "addNaForUnsCp",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpVec <- c(cpVec, NA)
  }
  .intSaveNm(
    "cpVecAfterUnsAdded",
    cpVec,
    indCombined,
    stageChnl,
    pathProject
  )

  cpVec
}
