# Local-FDR response estimate and final threshold
#
# Runs post-smoothing filtering, estimates the background-subtracted response
# proportion, and selects the empirical expression threshold that reproduces
# that estimate.

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

  trimObj <- .getCpUnsLocFilterAfterSmoothing(
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

.getCpUnsLocGetCpTrimSubset <- function(dataMod, keep) {
  .getCpUnsLocSubsetRows(dataMod, keep)
}
