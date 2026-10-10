# Local-FDR response estimate and final threshold
#
# Runs post-smoothing filtering, then sets the condition-level gate at the
# lower boundary of the filtered region (locThresholdMethod = "region"), at the
# empirical expression threshold whose background-subtracted frequency matches
# the sum of fitted response probabilities ("match"), or at the region boundary
# unless its frequency exceeds that sum by more than a factor ("cap").

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
  dataThreshold <- NULL
  densityBw <- attr(dataMod, "locDensityBw")
  method <- .getCpUnsLocThresholdMethod(chnlSettings)
  regionX <- NA_real_
  shiftedPeakRef <- attr(dataMod, "locShiftedPeakRef")

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
        regionX <- suppressWarnings(
          as.numeric(trimObj$info$final$xSum %||% NA_real_)[1L]
        )
        cpObj <- if (identical(method, "cap")) {
          .getCpUnsLocGetCpCap(
            dataThreshold = dataThreshold,
            regionX = regionX,
            cap = chnlSettings$locThresholdCap %||% 1.3,
            exTblStimNoMin = exTblStimNoMin,
            exTblUnsBias = exTblUnsBias,
            cpMin = cpMin,
            stage = stage,
            exTblStimOrig = exTblStimOrig,
            exTblUnsOrig = exTblUnsOrig,
            densityBw = densityBw
          )
        } else if (identical(method, "match")) {
          .getCpUnsLocGetCpActual(
            dataThreshold = dataThreshold,
            exTblStimNoMin = exTblStimNoMin,
            exTblUnsBias = exTblUnsBias,
            cpMin = cpMin,
            stage = stage,
            exTblStimOrig = exTblStimOrig,
            exTblUnsOrig = exTblUnsOrig,
            densityBw = densityBw
          )
        } else {
          .getCpUnsLocGetCpRegion(
            dataThreshold = dataThreshold,
            regionX = regionX,
            exTblStimNoMin = exTblStimNoMin,
            exTblUnsBias = exTblUnsBias,
            cpMin = cpMin,
            stage = stage
          )
        }
        .intSave(ind, stageChnl, pathProject, cpObj$cp)
        .debug("Completed loc gate for single sample") # nolint
      }
    }
  }

  attr(cpObj, "locThresholdMethod") <- method
  attr(cpObj, "locRegionX") <- regionX
  # Only when the shifted-peak rule was requested, so default outputs are
  # unchanged.
  if (!is.null(shiftedPeakRef)) {
    cpObj$locShiftedPeakRef <- isTRUE(shiftedPeakRef$applied)
    attr(cpObj, "locShiftedPeakInfo") <- shiftedPeakRef
  }
  # Carried with the gate to limit how far shared gates may lower it.
  cpObj$propBsEst <- .getCpUnsLocProbBsEst(dataThreshold)
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
  stage
) {
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
  exTblUnsOrig
) {
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
  densityBw = NULL
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

#' Resolve the local-FDR threshold method from channel settings
#'
#' Settings completed before the option existed have no value and use
#' "region", the default.
#' @keywords internal
.getCpUnsLocThresholdMethod <- function(chnlSettings) {
  method <- chnlSettings[["locThresholdMethod"]]
  if (.verifyIsNullOrNa(method)) "region" else method
}

#' Use the lower boundary of the filtered region as the gate
#'
#' The gate is the filtering boundary itself (`xSum`); cells strictly above it
#' are positive. The probability-sum frequency estimate is not matched, no
#' empirical cell value is selected and the gate is not moved below a cell.
#' The no-response and non-finite cases fall back exactly as in matching.
#'
#' @param regionX numeric Lower boundary of the filtered region.
#' @return list from `.getCpUnsLocConditionOut()`.
#' @keywords internal
.getCpUnsLocGetCpRegion <- function(
  dataThreshold,
  regionX,
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
  if (!is.finite(regionX)) {
    return(.getCpUnsLocConditionCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "No finite local-FDR region boundary"
    ))
  }
  .getCpUnsLocConditionOut(
    cp = regionX,
    locGenerated = TRUE,
    locGeneratedDirect = TRUE,
    locSource = "direct",
    locReason = "local_fdr_region_boundary_selected"
  )
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

#' Use the region boundary unless its frequency exceeds the capped estimate
#'
#' Keeps the region boundary (`regionX`, as under "region") when the
#' background-subtracted frequency strictly above it is at most `cap` times the
#' probability-sum estimate. Otherwise the gate is placed below the lowest
#' candidate cell value at or above the boundary whose frequency (cells at or
#' above it) is within the cap, as matching places its gate below the selected
#' cell. When no candidate is within the cap, the matching choice is used.
#'
#' @param cap numeric Largest allowed ratio of frequency to estimate.
#' @return list from `.getCpUnsLocConditionOut()`. A selected cell is attribute
#'   `cpSelected`; `locCapExceededAbove` records whether any higher candidate's
#'   frequency exceeds the cap again.
#' @keywords internal
.getCpUnsLocGetCpCap <- function(
  dataThreshold,
  regionX,
  cap,
  exTblStimNoMin,
  exTblUnsBias,
  cpMin,
  stage,
  exTblStimOrig,
  exTblUnsOrig,
  densityBw = NULL
) {
  cpRegion <- .getCpUnsLocGetCpRegion(
    dataThreshold = dataThreshold,
    regionX = regionX,
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsBias = exTblUnsBias,
    cpMin = cpMin,
    stage = stage
  )
  if (!isTRUE(cpRegion$locGenerated)) {
    return(cpRegion)
  }

  limit <- cap * .getCpUnsLocProbBsEst(dataThreshold)
  xStim <- .getCut(exTblStimOrig)
  xUns <- .getCut(exTblUnsOrig)
  propBsRegion <- mean(xStim > regionX) - mean(xUns > regionX)
  x <- .getCut(dataThreshold)
  propBs <- dataThreshold$propBs
  if (is.finite(propBsRegion) && propBsRegion <= limit) {
    cpRegion$locReason <- "local_fdr_cap_region_boundary_selected"
    attr(cpRegion, "locCapExceededAbove") <- any(
      propBs[x > regionX] > limit,
      na.rm = TRUE
    )
    return(cpRegion)
  }

  within <- which(is.finite(x) & x >= regionX & propBs <= limit)
  if (length(within) == 0L) {
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
    if (isTRUE(cpObj$locGenerated)) {
      cpObj$locReason <- "local_fdr_cap_match_fallback"
    }
    attr(cpObj, "locCapExceededAbove") <- NA
    return(cpObj)
  }
  cpVal <- min(x[within])

  cpObj <- .getCpUnsLocConditionOut(
    cp = .getCpUnsLocGateBelowCell(
      cp = cpVal,
      x = c(xStim, xUns),
      densityBw = densityBw
    ),
    locGenerated = TRUE,
    locGeneratedDirect = TRUE,
    locSource = "direct",
    locReason = "local_fdr_cap_threshold_selected"
  )
  attr(cpObj, "cpSelected") <- cpVal
  attr(cpObj, "locCapExceededAbove") <- any(propBs[x > cpVal] > limit, na.rm = TRUE)
  cpObj
}
