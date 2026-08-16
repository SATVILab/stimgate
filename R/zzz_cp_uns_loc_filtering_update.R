# Current post-smoothing filtering logic
#
# This file overrides the ordinary post-smoothing filtering entry point defined
# in cp_uns_loc_filtering.R. It is loaded last so the existing shape-enforced
# route can be retained unchanged while the ordinary route follows the current
# method description: clear-response boundary, contiguous density dominance,
# marginal quality, and optional taut-string antimodality. The old global
# derivative filter is not used.

# Keep the pre-existing implementation available for the separate
# locEnforceShapeThreshold route, which is intentionally not changed here.
.getCpUnsLocFilterAfterSmoothingLegacy <- .getCpUnsLocFilterAfterSmoothing

#' Apply current post-smoothing filtering for the ordinary local-FDR route
#' @keywords internal
.getCpUnsLocFilterAfterSmoothing <- function(
  dataMod,
  exTblStimNoMin,
  exTblUnsBias,
  cpMin,
  stage,
  chnlSettings
) {
  force(stage)

  # Leave the separate shape-enforced arm exactly as it was.
  if (isTRUE(attr(dataMod, "locShapeThresholdRequested"))) {
    return(.getCpUnsLocFilterAfterSmoothingLegacy(
      dataMod = dataMod,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      cpMin = cpMin,
      stage = stage,
      chnlSettings = chnlSettings
    ))
  }

  info <- list(applied = FALSE, reason = "not_filtered")

  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    info$reason <- "no_data_mod"
    return(list(dataMod = dataMod, cp = NULL, info = info))
  }

  dataMod <- .getCpUnsLocSubsetRows(dataMod, order(.getCut(dataMod)))
  dataModFull <- dataMod

  # Values below minProbXPos were retained only as a lower smoothing margin.
  # The preliminary filtering stage remains the hard lower domain for the
  # ordinary route.
  preliminaryLowerBoundX <- suppressWarnings(
    as.numeric(attr(dataModFull, "minProbXPos"))[1L]
  )
  info$minProbXPos <- preliminaryLowerBoundX
  info$clampingCellsRetainedForThresholdSelection <- TRUE

  probCol <- .getCpUnsLocProbabilityColumn(dataModFull, chnlSettings)
  prob <- .getCpUnsLocProbability(dataModFull, probCol)
  maxProb <- suppressWarnings(max(prob, na.rm = TRUE))
  minPeakProb <- .getCpUnsLocSetting(chnlSettings, "locMinPeakProb", 0.25)
  minPeakProb <- .getCpUnsLocUnitValue(minPeakProb, 0.25, allowZero = TRUE)

  info$probCol <- probCol
  info$maxProb <- maxProb
  info$minPeakProb <- minPeakProb

  if (!is.finite(maxProb) || maxProb < minPeakProb) {
    info$applied <- TRUE
    info$reason <- "max_response_probability_below_minimum"
    return(list(
      dataMod = dataModFull[0, , drop = FALSE],
      cp = .getCpUnsLocConditionCpNonLoc(
        cpMin = cpMin,
        exTblStimNoMin = exTblStimNoMin,
        exTblUnsBias = exTblUnsBias
      ),
      info = info
    ))
  }

  # x_clear_init: conservative response-probability boundary. Keep the current
  # derivative machinery and its 85%-of-maximum fallback.
  thresholdClearInit <- .getCpUnsLocStageThreshold(
    dataModFull,
    chnlSettings,
    probCol,
    stage = "marginal"
  )
  if (!is.finite(thresholdClearInit$thresholdX)) {
    thresholdClearInit <- .getCpUnsLocHighProbabilityReference(
      dataMod = dataModFull,
      probCol = probCol,
      fraction = 0.85,
      derivativeInfo = thresholdClearInit$info,
      shapeApplied = FALSE
    )
  }

  xClearInit <- suppressWarnings(
    as.numeric(thresholdClearInit$thresholdX)[1L]
  )
  if (!is.finite(xClearInit)) {
    info$applied <- TRUE
    info$reason <- "no_informative_clear_response_reference"
    info$thresholdClass <- "undefined"
    info$shareable <- FALSE
    return(list(
      dataMod = dataModFull,
      cp = .getCpUnsLocConditionCpNonLoc(
        cpMin = cpMin,
        exTblStimNoMin = exTblStimNoMin,
        exTblUnsBias = exTblUnsBias
      ),
      info = info
    ))
  }

  # Density dominance can move x_clear_init lower. The technical implementation
  # retains the existing tail collapse and half-bandwidth smoothing, but the
  # dominance region is explicitly the contiguous dominant region containing
  # x_clear_init. The peak is searched only between its left onset and
  # x_clear_init.
  dominance <- .getCpUnsLocDominanceBoundaryCurrent(
    density = attr(dataModFull, "locDensityComparison"),
    startX = xClearInit,
    densityBw = attr(dataModFull, "locDensityBw"),
    dominanceRatio = 2,
    onsetWeight = 2 / 3,
    lowerBoundX = preliminaryLowerBoundX
  )
  xDom <- suppressWarnings(as.numeric(dominance$startX)[1L])
  xClear <- .getCpUnsLocFiniteMin(c(xClearInit, xDom))

  info$clear <- list(
    xClearInit = xClearInit,
    derivativeThreshold = thresholdClearInit$info,
    xDom = xDom,
    densityDominance = dominance$info,
    xClear = xClear
  )
  # Retain the older diagnostic name for downstream inspection of saved
  # intermediates, but there is no global filter in this route.
  info$marginalReference <- list(
    referenceX = xClear,
    derivativeThreshold = thresholdClearInit$info,
    densityDominance = dominance$info,
    globalFilterApplied = FALSE
  )

  if (!is.finite(xClear)) {
    info$applied <- TRUE
    info$reason <- "no_informative_clear_response_reference"
    info$thresholdClass <- "undefined"
    info$shareable <- FALSE
    return(list(
      dataMod = dataModFull,
      cp = .getCpUnsLocConditionCpNonLoc(
        cpMin = cpMin,
        exTblStimNoMin = exTblStimNoMin,
        exTblUnsBias = exTblUnsBias
      ),
      info = info
    ))
  }

  # x_qual: scan bins leftward from x_clear. The old global derivative floor is
  # deliberately absent. The only lower bound is the preliminary filter that
  # was already applied before smoothing.
  quality <- .getCpUnsLocQualityBoundaryCurrent(
    dataMod = dataModFull,
    chnlSettings = chnlSettings,
    probCol = probCol,
    xClear = xClear,
    lowerBoundX = preliminaryLowerBoundX
  )
  xQual <- suppressWarnings(as.numeric(quality$thresholdX)[1L])
  if (!is.finite(xQual)) {
    xQual <- xClear
  }
  info$marginal <- quality$info

  xBase <- .getCpUnsLocFiniteMin(c(xClear, xQual))

  # Antimodality is an additional criterion for moving the already-supported
  # boundary lower. Candidate antimodes must lie below xBase. The reference
  # mode is the highest taut-string mode to the right of xBase, including a
  # mode between xQual and xClear when xQual < xClear. No additional mode
  # validity criterion is applied. A candidate trough must fall to <95% of the
  # reference-mode height so a numerically flat density cannot trigger a move.
  antimode <- .getCpUnsLocAntimodeBoundaryCurrent(
    dataMod = dataModFull,
    chnlSettings = chnlSettings,
    xBase = xBase,
    lowerBoundX = preliminaryLowerBoundX,
    heightFrac = .getCpUnsLocSetting(
      chnlSettings,
      "locAntimodeReferenceHeightFrac",
      0.95
    )
  )
  xAntimode <- suppressWarnings(as.numeric(antimode$thresholdX)[1L])
  info$antimode <- antimode$info

  xSum <- .getCpUnsLocFiniteMin(c(xClear, xQual, xAntimode))
  info$final <- list(
    xClearInit = xClearInit,
    xDom = xDom,
    xClear = xClear,
    xQual = xQual,
    xBase = xBase,
    xAntimode = xAntimode,
    xSum = xSum,
    globalFilterApplied = FALSE
  )

  if (!is.finite(xSum)) {
    info$applied <- TRUE
    info$reason <- "no_final_post_smoothing_boundary"
    info$thresholdClass <- "undefined"
    info$shareable <- FALSE
    return(list(
      dataMod = dataModFull,
      cp = .getCpUnsLocConditionCpNonLoc(
        cpMin = cpMin,
        exTblStimNoMin = exTblStimNoMin,
        exTblUnsBias = exTblUnsBias
      ),
      info = info
    ))
  }

  x <- suppressWarnings(as.numeric(.getCut(dataModFull)))
  keep <- is.finite(x) & x >= xSum
  dataModOut <- .getCpUnsLocSubsetRows(dataModFull, keep)
  attr(dataModOut, "locFinalFilterX") <- xSum

  nDropped <- sum(!keep)
  info$applied <- nDropped > 0L
  info$nDropped <- nDropped
  info$reason <- if (info$applied) {
    "filtered_at_lowest_supported_post_smoothing_boundary"
  } else {
    "post_smoothing_boundary_kept_all_values"
  }
  info$lowProbLeft <- list(
    applied = info$applied,
    global = list(
      applied = FALSE,
      reason = "global_filter_removed_from_current_method"
    ),
    marginal = quality$info
  )
  info$flatLeft <- info$lowProbLeft

  if (!is.data.frame(dataModOut) || nrow(dataModOut) == 0L) {
    info$reason <- "all_cells_removed_by_post_smoothing_filtering"
    return(.getCpUnsLocEmptyFilterResult(
      dataMod = dataModOut,
      info = info,
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias
    ))
  }

  list(dataMod = dataModOut, cp = NULL, info = info)
}

#' Identify the current contiguous density-dominance boundary
#'
#' The stimulated and unstimulated densities are first stabilised in the same
#' way as the earlier implementation: each density is collapsed to its mean to
#' the right of x_clear_init, then both are Gaussian-kernel smoothed with half
#' the original fixed density bandwidth. The 2:1 dominance region used here is
#' specifically the contiguous region containing x_clear_init. Its left onset
#' is tempered by the strongest dominance point between that onset and
#' x_clear_init.
#' @keywords internal
.getCpUnsLocDominanceBoundaryCurrent <- function(
  density,
  startX,
  densityBw = NULL,
  dominanceRatio = 2,
  onsetWeight = 2 / 3,
  lowerBoundX = NA_real_
) {
  info <- list(
    applied = FALSE,
    reason = "density_dominance_boundary_unavailable",
    clearInitX = startX,
    dominanceRatio = dominanceRatio,
    onsetWeight = onsetWeight,
    lowerBoundX = lowerBoundX
  )

  if (
    !is.data.frame(density) ||
      !all(c("x", "stim", "unstim") %in% names(density)) ||
      !is.finite(startX) ||
      !is.finite(dominanceRatio) ||
      dominanceRatio <= 0 ||
      !is.finite(onsetWeight) ||
      onsetWeight < 0 ||
      onsetWeight > 1
  ) {
    info$reason <- "invalid_density_dominance_input"
    return(list(startX = NA_real_, info = info))
  }

  if (is.list(densityBw)) {
    info$reason <- "adaptive_density_bandwidth_not_supported_for_dominance"
    return(list(startX = NA_real_, info = info))
  }
  densityBw <- suppressWarnings(as.numeric(densityBw))
  densityBw <- densityBw[is.finite(densityBw) & densityBw > 0]
  if (length(densityBw) == 0L) {
    info$reason <- "fixed_density_bandwidth_unavailable_for_dominance"
    return(list(startX = NA_real_, info = info))
  }
  smoothBw <- 0.5 * densityBw[[1L]]

  x <- suppressWarnings(as.numeric(density$x))
  stim <- suppressWarnings(as.numeric(density$stim))
  unstim <- suppressWarnings(as.numeric(density$unstim))
  keep <- is.finite(x) & is.finite(stim) & is.finite(unstim)
  x <- x[keep]
  stim <- pmax(0, stim[keep])
  unstim <- pmax(0, unstim[keep])
  ord <- order(x)
  x <- x[ord]
  stim <- stim[ord]
  unstim <- unstim[ord]

  lowerBoundX <- suppressWarnings(as.numeric(lowerBoundX)[1L])
  if (is.finite(lowerBoundX)) {
    inDomain <- x >= lowerBoundX
    x <- x[inDomain]
    stim <- stim[inDomain]
    unstim <- unstim[inDomain]
  }

  if (length(x) < 3L || anyDuplicated(x) || diff(range(x)) <= 0) {
    info$reason <- "insufficient_density_dominance_grid"
    return(list(startX = NA_real_, info = info))
  }

  # Preserve the previous technical stabilisation of the already-clear tail.
  tailIdx <- which(x >= startX)
  if (length(tailIdx) > 0L) {
    stim[tailIdx] <- mean(stim[tailIdx])
    unstim[tailIdx] <- mean(unstim[tailIdx])
  }

  z <- outer(x, x, "-") / smoothBw
  weights <- exp(-0.5 * z^2)
  weightSum <- rowSums(weights)
  stimSmooth <- as.numeric(weights %*% stim) / weightSum
  unstimSmooth <- as.numeric(weights %*% unstim) / weightSum
  totalSmooth <- stimSmooth + unstimSmooth
  score <- stimSmooth / totalSmooth
  supportMin <- sqrt(.Machine$double.eps) * max(totalSmooth, na.rm = TRUE)
  validScore <- is.finite(score) &
    is.finite(totalSmooth) &
    totalSmooth > supportMin
  scoreThreshold <- dominanceRatio / (1 + dominanceRatio)
  dominant <- validScore & score >= scoreThreshold

  # x_clear_init anchors the contiguous region. Use the grid point immediately
  # to its left when possible so the search never begins beyond the stated
  # conservative boundary.
  refCandidate <- which(x <= startX)
  if (length(refCandidate) == 0L) {
    info$reason <- "clear_response_reference_left_of_density_grid"
    return(list(startX = NA_real_, info = info))
  }
  refIdx <- max(refCandidate)
  if (!isTRUE(dominant[[refIdx]])) {
    info$reason <- "density_not_dominant_at_clear_response_reference"
    info$smoothBw <- smoothBw
    info$scoreThreshold <- scoreThreshold
    info$referenceGridX <- x[[refIdx]]
    info$referenceScore <- score[[refIdx]]
    return(list(startX = NA_real_, info = info))
  }

  # Move left only through the contiguous dominant run that contains x_clear_init.
  leftIdx <- refIdx
  while (
    leftIdx > 1L &&
      isTRUE(validScore[[leftIdx - 1L]]) &&
      isTRUE(dominant[[leftIdx - 1L]])
  ) {
    leftIdx <- leftIdx - 1L
  }

  if (leftIdx == 1L || !isTRUE(validScore[[leftIdx - 1L]])) {
    info$reason <- "density_dominance_has_no_observed_left_onset"
    info$smoothBw <- smoothBw
    info$scoreThreshold <- scoreThreshold
    return(list(startX = NA_real_, info = info))
  }

  nonDomIdx <- leftIdx - 1L
  scoreChange <- score[[leftIdx]] - score[[nonDomIdx]]
  onsetX <- if (!is.finite(scoreChange) || scoreChange == 0) {
    x[[leftIdx]]
  } else {
    x[[nonDomIdx]] +
      (scoreThreshold - score[[nonDomIdx]]) *
        (x[[leftIdx]] - x[[nonDomIdx]]) / scoreChange
  }

  # The robustness peak is explicitly constrained to the interval from onset
  # through x_clear_init, as described in the current method.
  peakCandidate <- seq.int(leftIdx, refIdx)
  peakCandidate <- peakCandidate[validScore[peakCandidate]]
  if (length(peakCandidate) == 0L) {
    info$reason <- "density_dominance_peak_unavailable_before_clear_reference"
    return(list(startX = NA_real_, info = info))
  }
  peakIdx <- peakCandidate[[which.max(score[peakCandidate])]]
  peakX <- x[[peakIdx]]
  dominanceX <- onsetWeight * onsetX + (1 - onsetWeight) * peakX

  info$applied <- is.finite(dominanceX)
  info$reason <- "identified_contiguous_density_dominance_boundary"
  info$tailCollapsed <- length(tailIdx) > 0L
  info$tailCollapseStartX <- if (length(tailIdx) > 0L) {
    x[[min(tailIdx)]]
  } else {
    NA_real_
  }
  info$smoothBw <- smoothBw
  info$scoreThreshold <- scoreThreshold
  info$referenceGridX <- x[[refIdx]]
  info$onsetX <- onsetX
  info$peakX <- peakX
  info$peakScore <- score[[peakIdx]]
  info$dominanceStartX <- dominanceX
  info$scoreTbl <- data.frame(
    x = x,
    stim = stimSmooth,
    unstim = unstimSmooth,
    score = score,
    dominant = dominant
  )

  list(startX = dominanceX, info = info)
}

#' Obtain the quality-based lower boundary starting at x_clear
#' @keywords internal
.getCpUnsLocQualityBoundaryCurrent <- function(
  dataMod,
  chnlSettings,
  probCol,
  xClear,
  lowerBoundX = NA_real_
) {
  out <- .getCpUnsLocFilterMarginalBins(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol,
    startX = xClear,
    lowerBoundX = lowerBoundX
  )

  xQual <- suppressWarnings(as.numeric(out$info$finalStartX)[1L])
  if (!is.finite(xQual)) {
    xQual <- xClear
  }

  out$info$thresholdX <- xClear
  out$info$purityStartX <- xClear
  out$info$referenceBasis <- "x_clear"
  out$info$thresholdClass <- "informative"
  out$info$shareable <- TRUE
  out$info$preliminaryLowerBoundX <- lowerBoundX
  out$info$xQual <- xQual

  list(thresholdX = xQual, dataMod = out$dataMod, info = out$info)
}

#' Select a taut-string antimode that can move xBase lower
#' @keywords internal
.getCpUnsLocAntimodeBoundaryCurrent <- function(
  dataMod,
  chnlSettings,
  xBase,
  lowerBoundX = NA_real_,
  heightFrac = 0.95
) {
  info <- list(
    applied = FALSE,
    reason = "antimode_boundary_not_applied",
    xBase = xBase,
    lowerBoundX = lowerBoundX,
    heightFrac = heightFrac
  )

  if (!is.data.frame(dataMod) || nrow(dataMod) < 5L || !is.finite(xBase)) {
    info$reason <- "insufficient_data_for_antimode_boundary"
    return(list(thresholdX = NA_real_, info = info))
  }

  heightFrac <- suppressWarnings(as.numeric(heightFrac)[1L])
  if (!is.finite(heightFrac) || heightFrac <= 0 || heightFrac >= 1) {
    heightFrac <- 0.95
  }
  info$heightFrac <- heightFrac

  expr <- suppressWarnings(as.numeric(.getCut(dataMod)))
  keepExpr <- is.finite(expr)
  lowerBoundX <- suppressWarnings(as.numeric(lowerBoundX)[1L])
  if (is.finite(lowerBoundX)) {
    keepExpr <- keepExpr & expr >= lowerBoundX
  }
  exprForDensity <- expr[keepExpr]

  if (length(exprForDensity) < 5L || length(unique(exprForDensity)) < 3L) {
    info$reason <- "too_few_values_for_antimode_density"
    return(list(thresholdX = NA_real_, info = info))
  }

  density <- .getCpUnsLocAntimodeDensity(
    expr = exprForDensity,
    chnlSettings = chnlSettings,
    originalBw = attr(dataMod, "locDensityBw"),
    mtd = "taut_string"
  )
  if (is.null(density)) {
    info$reason <- "antimode_density_failed"
    return(list(thresholdX = NA_real_, info = info))
  }

  extrema <- .getCpUnsLocTautStringExtremaCurrent(density)
  modes <- extrema$modes
  antimodes <- extrema$antimodes

  info$nExpressionValuesForDensity <- length(exprForDensity)
  info$densityMtd <- attr(density, "locDensityMtd")
  info$allModes <- modes
  info$allAntimodes <- antimodes

  if (nrow(antimodes) == 0L) {
    info$reason <- "no_taut_string_antimodes"
    return(list(thresholdX = NA_real_, info = info))
  }

  candidate <- antimodes[antimodes$x < xBase, , drop = FALSE]
  if (is.finite(lowerBoundX)) {
    candidate <- candidate[candidate$x >= lowerBoundX, , drop = FALSE]
  }
  info$candidateAntimodes <- candidate

  if (nrow(candidate) == 0L) {
    info$reason <- "no_antimode_below_supported_boundary"
    return(list(thresholdX = NA_real_, info = info))
  }

  referenceModes <- modes[modes$x > xBase, , drop = FALSE]
  info$referenceModes <- referenceModes
  if (nrow(referenceModes) == 0L) {
    info$reason <- "no_taut_string_mode_right_of_supported_boundary"
    return(list(thresholdX = NA_real_, info = info))
  }

  refIdx <- which.max(referenceModes$height)
  referenceMode <- referenceModes[refIdx, , drop = FALSE]
  referenceHeight <- referenceMode$height[[1L]]
  heightLimit <- heightFrac * referenceHeight

  passing <- candidate[
    is.finite(candidate$height) & candidate$height < heightLimit,
    ,
    drop = FALSE
  ]

  info$referenceMode <- referenceMode
  info$referenceModeHeight <- referenceHeight
  info$heightLimit <- heightLimit
  info$passingAntimodes <- passing

  if (nrow(passing) == 0L) {
    info$reason <- "no_antimode_dropped_below_reference_mode_height"
    return(list(thresholdX = NA_real_, info = info))
  }

  selected <- passing[which.max(passing$x), , drop = FALSE]
  thresholdX <- selected$x[[1L]]

  info$applied <- is.finite(thresholdX)
  info$reason <- "selected_rightmost_valid_antimode"
  info$selectedAntimode <- selected
  info$thresholdX <- thresholdX

  list(thresholdX = thresholdX, info = info)
}

#' Return modes and antimodes of a piecewise-constant taut-string density
#' @keywords internal
.getCpUnsLocTautStringExtremaCurrent <- function(density) {
  empty <- data.frame(x = numeric(0L), height = numeric(0L))
  if (is.null(density) || !identical(density$method, "taut_string")) {
    return(list(modes = empty, antimodes = empty))
  }

  x <- suppressWarnings(as.numeric(density$x))
  y <- suppressWarnings(as.numeric(density$y))
  finite <- is.finite(x) & is.finite(y)
  x <- x[finite]
  y <- y[finite]

  if (length(x) != length(y) || length(y) < 3L) {
    return(list(modes = empty, antimodes = empty))
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  runId <- cumsum(c(
    TRUE,
    !dplyr::near(y[-1L], y[-length(y)])
  ))
  runs <- split(seq_along(y), runId)
  runY <- vapply(runs, function(i) y[i[[1L]]], numeric(1L))
  runLeft <- vapply(runs, function(i) min(x[i]), numeric(1L))
  runRight <- vapply(runs, function(i) max(x[i]), numeric(1L))
  runX <- (runLeft + runRight) / 2

  if (length(runY) < 3L) {
    return(list(modes = empty, antimodes = empty))
  }

  internal <- seq.int(2L, length(runY) - 1L)
  modeIdx <- internal[
    runY[internal] > runY[internal - 1L] &
      runY[internal] > runY[internal + 1L]
  ]
  antimodeIdx <- internal[
    runY[internal] < runY[internal - 1L] &
      runY[internal] < runY[internal + 1L]
  ]

  modes <- data.frame(
    x = runX[modeIdx],
    height = runY[modeIdx]
  )
  antimodes <- data.frame(
    x = runX[antimodeIdx],
    height = runY[antimodeIdx]
  )

  list(
    modes = modes[order(modes$x), , drop = FALSE],
    antimodes = antimodes[order(antimodes$x), , drop = FALSE]
  )
}
