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
    chnlSettings = list()) {
  ind <- .getInd(exTblStimNoMin)
  chnl <- .getCpUnsLocGetChnl(exTblStimNoMin)
  stageChnl <- file.path(stage, chnl)
  dataThreshold <- NULL
  densityBw <- attr(dataMod, "locDensityBw")

  if (!is.data.frame(dataMod)) {
    .intSaveNm("noDataModDf", NULL, ind, stageChnl, pathProject)
    .intSaveNm("cpInd", dataMod, ind, stageChnl, pathProject)
    if (is.list(dataMod) && "cp" %in% names(dataMod)) {
      cpObj <- dataMod
      cpObj$locGenerated <- cpObj$locGenerated %||% FALSE
      cpObj$locGeneratedDirect <- cpObj$locGeneratedDirect %||% FALSE
      cpObj$locSource <- cpObj$locSource %||% "not_calculated"
      cpObj$locReason <- cpObj$locReason %||% "data_mod_not_available"
    } else {
      cpObj <- .getCpUnsLocConditionOut(
        cp = .getCpUnsLocConditionCpNonLoc(
          cpMin = cpMin,
          exTblStimNoMin = exTblStimNoMin,
          exTblUnsBias = exTblUnsBias
        ),
        locGenerated = FALSE,
        locGeneratedDirect = FALSE,
        locSource = "not_calculated",
        locReason = "data_mod_not_available"
      )
    }
  } else {
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
      cpObj <- .getCpUnsLocConditionOut(
        cp = cpInd,
        locGenerated = FALSE,
        locGeneratedDirect = FALSE,
        locSource = "not_calculated",
        locReason = trimObj$info$reason %||% "trim_returned_non_local_cutpoint"
      )
    } else {
      dataMod <- trimObj$dataMod
      if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
        cpInd <- .getCpUnsLocConditionCpNonLoc(
          cpMin = cpMin,
          exTblStimNoMin = exTblStimNoMin,
          exTblUnsBias = exTblUnsBias
        )
        .intSave(ind, stageChnl, pathProject, cpInd)
        .debug("Completed loc gate for single sample") # nolint
        cpObj <- .getCpUnsLocConditionOut(
          cp = cpInd,
          locGenerated = FALSE,
          locGeneratedDirect = FALSE,
          locSource = "not_calculated",
          locReason = "empty_data_mod_after_trimming"
        )
      } else {
        dataThreshold <- .getCpUnsLocGetCpDataThreshold(
          dataMod = dataMod,
          exTblStimOrig = exTblStimOrig,
          exTblStimNoMin = exTblStimNoMin,
          exTblUnsOrig = exTblUnsOrig,
          stage = stage,
          pathProject = pathProject
        )
        .intSave(ind, stageChnl, pathProject, dataThreshold)
        cpObj <- .getCpUnsLocGetCpActual(
          dataThreshold = dataThreshold,
          exTblStimNoMin = exTblStimNoMin,
          exTblUnsBias = exTblUnsBias,
          cpMin = cpMin,
          stage = stage,
          exTblStimOrig = exTblStimOrig,
          exTblUnsOrig = exTblUnsOrig,
          densityBw = densityBw
        )
        .intSave(ind, stageChnl, pathProject, cpObj$cp)
        .debug("Completed loc gate for single sample") # nolint
      }
    }
  }

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
  cpObj
}


#' @keywords internal
.getCpUnsLocGetCpDataThreshold <- function(
    dataMod,
    exTblStimOrig,
    exTblStimNoMin,
    exTblUnsOrig,
    pathProject,
    stage) {
  # Remove the lower-margin values retained only to anchor the smoother at the
  # point where the final response proportion is calculated, not while the
  # filtering thresholds are identified.
  dataModEstimate <- .getCpUnsLocGetCpDataThresholdExcludeMargin(dataMod)

  dataCount <- .getCpUnsLocGetCpDataThresholdCount(dataModEstimate)
  probBsEst <- sum(dataCount$pred) / nrow(exTblStimOrig)
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
    exTblUnsOrig = exTblUnsOrig
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
  .getCpUnsLocSubsetRows(
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

.getCpUnsLocTailPropAtThresholds <- function(x, thresholds, denominator) {
  x <- as.numeric(x)
  thresholds <- as.numeric(thresholds)

  # Preserve the current behaviour if the expression vector contains NA/NaN:
  # sum(x >= threshold) without na.rm = TRUE would return NA.
  if (anyNA(x)) {
    return(rep(NA_real_, length(thresholds)))
  }

  out <- rep(NA_real_, length(thresholds))
  ok <- !is.na(thresholds)

  if (!any(ok)) {
    return(out)
  }

  x <- sort(x)

  # left.open = TRUE gives the number of x values STRICTLY below
  # each threshold. Therefore:
  #
  #   length(x) - n_below
  #
  # is the number >= threshold.
  n_below <- findInterval(
    thresholds[ok],
    x,
    left.open = TRUE
  )

  out[ok] <- (length(x) - n_below) / denominator
  out
}


.getCpUnsLocGetCpDataThresholdActual <- function(
    dataCount,
    propBsEst,
    exTblStimOrig,
    exTblUnsOrig) {
  thresholds <- .getCut(dataCount)

  propStimVec <- .getCpUnsLocTailPropAtThresholds(
    x = .getCut(exTblStimOrig),
    thresholds = thresholds,
    denominator = nrow(exTblStimOrig)
  )

  propUnsVec <- .getCpUnsLocTailPropAtThresholds(
    x = .getCut(exTblUnsOrig),
    thresholds = thresholds,
    denominator = nrow(exTblUnsOrig)
  )

  dataCount |>
    dplyr::mutate(
      propStim = propStimVec,
      propUns = propUnsVec
    ) |>
    dplyr::mutate(
      propBs = propStim - propUns,
      propBsDiff = propBs - propBsEst
    )
}

#' Select the local-FDR threshold and place the gate below it
#'
#' The selected threshold is the candidate cell value whose tail frequency
#' (cells at or above it) best matches the estimate. Gates are applied with a
#' strict `x > gate`, so the gate is placed below that cell, within the gap
#' to the next lower stimulated or unstimulated cell: it moves down by the
#' smaller of twice the density bandwidth and half that gap. The applied gate
#' then counts exactly the cells the selection counted.
#'
#' @param exTblStimOrig,exTblUnsOrig data.frame or NULL Expression used to
#'   find the gap below the selected cell; NULL keeps the gate at the cell.
#' @param densityBw numeric, list or NULL Local-FDR density bandwidth; an
#'   adaptive bandwidth uses the shared bandwidth at the selected cell.
#' @return list from `.getCpUnsLocConditionOut()`, with the selected cell value
#'   as attribute `cpSelected`.
#' @keywords internal
.getCpUnsLocGetCpActual <- function(
    dataThreshold,
    exTblStimNoMin,
    exTblUnsBias,
    cpMin,
    stage,
    exTblStimOrig = NULL,
    exTblUnsOrig = NULL,
    densityBw = NULL) {
  if (nrow(dataThreshold) == 0L) {
    return(.getCpUnsLocConditionCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "Too few responding cells"
    ))
  }
  bestIdx <- which.min(abs(dataThreshold$propBsDiff))
  cpVal <- .getCut(dataThreshold)[bestIdx]

  if (!is.finite(cpVal)) {
    return(.getCpUnsLocConditionCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "No finite local-FDR threshold"
    ))
  }

  cpObj <- .getCpUnsLocConditionOut(
    cp = .getCpUnsLocGateBelowCell(
      cp = cpVal,
      x = c(.getCut(exTblStimOrig), .getCut(exTblUnsOrig)),
      densityBw = densityBw
    ),
    locGenerated = TRUE,
    locGeneratedDirect = TRUE,
    locSource = "direct",
    locReason = "local_fdr_threshold_selected"
  )
  attr(cpObj, "cpSelected") <- cpVal
  cpObj
}

#' Place a gate in the gap below a selected cell value
#'
#' @param cp numeric Selected cell value.
#' @param x numeric Stimulated and unstimulated expression; NULL keeps `cp`.
#' @param densityBw numeric, list or NULL Density bandwidth. An adaptive
#'   bandwidth object is interpolated at `cp`; an unavailable bandwidth leaves
#'   only the half-gap limit.
#' @return numeric Gate strictly below `cp` and above every lower value of `x`.
#' @keywords internal
.getCpUnsLocGateBelowCell <- function(cp, x, densityBw = NULL) {
  if (is.null(x) || length(x) == 0L) {
    return(cp)
  }
  below <- x[is.finite(x) & x < cp]
  halfGap <- if (length(below) > 0L) (cp - max(below)) / 2 else Inf
  if (is.list(densityBw)) {
    grid <- densityBw$grid
    bwGrid <- densityBw$sharedGrid
    valid <- if (length(grid) == length(bwGrid)) {
      is.finite(grid) & is.finite(bwGrid) & bwGrid > 0
    } else {
      logical(0)
    }
    bw <- if (sum(valid) > 1L) {
      stats::approx(grid[valid], bwGrid[valid], xout = cp, rule = 2)$y
    } else {
      NA_real_
    }
  } else {
    bw <- suppressWarnings(as.numeric(densityBw)[1L])
  }
  step <- min(halfGap, if (is.finite(bw) && bw > 0) 2 * bw else Inf)
  if (!is.finite(step)) {
    # No lower cell and no bandwidth: move just below the cell.
    step <- sqrt(.Machine$double.eps) * max(1, abs(cp))
  }
  cp - step
}
