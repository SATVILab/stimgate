# Appendix-aligned post-smoothing filtering

# Local-FDR filtering after probability smoothing
#
# The names in this file follow the appendix terminology:
#   alpha: minimum relative height for an eligible derivative peak
#   omega: minimum response probability associated with that peak
#   psi: fraction of the selected peak used to locate the rising edge
#
# Density shape, when requested, is enforced before smoothing. Post-smoothing
# filtering is therefore sequentially global and marginal. A missing global
# threshold does not prevent the marginal stage from being attempted.

#' Apply all filtering steps after smoothing
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
  info <- list(applied = FALSE, reason = "not_filtered")

  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    info$reason <- "no_data_mod"
    return(list(dataMod = dataMod, cp = NULL, info = info))
  }

  # These cells anchor the smoother and remain available while the filtering
  # thresholds are identified. They are removed only when the final response
  # proportion is calculated.
  info$minProbXPos <- attr(dataMod, "minProbXPos")
  info$clampingCellsRetainedForThresholdSelection <- TRUE
  shapeRequested <- isTRUE(attr(dataMod, "locShapeThresholdRequested"))
  shapeApplied <- isTRUE(attr(dataMod, "locShapeThresholdApplied"))
  shapeLowerBoundX <- suppressWarnings(
    as.numeric(attr(dataMod, "locShapeThresholdX"))[1L]
  )
  shapeTailgateX <- suppressWarnings(
    as.numeric(attr(dataMod, "locShapeTailgateX"))[1L]
  )
  info$shape <- list(
    requested = shapeRequested,
    applied = shapeApplied,
    lowerBoundX = shapeLowerBoundX,
    adjustedTailgateX = shapeTailgateX,
    detail = attr(dataMod, "locShapeThresholdInfo")
  )

  if (shapeRequested && !shapeApplied) {
    info$applied <- TRUE
    info$reason <- "requested_shape_threshold_unavailable"
    info$thresholdClass <- "undefined"
    info$shareable <- FALSE
    return(list(
      dataMod = dataMod,
      cp = .getCpUnsLocConditionCpNonLoc(
        cpMin = cpMin,
        exTblStimNoMin = exTblStimNoMin,
        exTblUnsBias = exTblUnsBias
      ),
      info = info
    ))
  }

  dataMod <- .getCpUnsLocSubsetRows(dataMod, order(.getCut(dataMod)))
  probCol <- .getCpUnsLocProbabilityColumn(dataMod, chnlSettings)
  prob <- .getCpUnsLocProbability(dataMod, probCol)
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
      dataMod = dataMod[0, , drop = FALSE],
      cp = .getCpUnsLocConditionCpNonLoc(
        cpMin = cpMin,
        exTblStimNoMin = exTblStimNoMin,
        exTblUnsBias = exTblUnsBias
      ),
      info = info
    ))
  }
  thresholdGlobal <- .getCpUnsLocStageThreshold(
    dataMod,
    chnlSettings,
    probCol,
    stage = "global"
  )
  thresholdMarginal <- .getCpUnsLocStageThreshold(
    dataMod,
    chnlSettings,
    probCol,
    stage = "marginal"
  )
  if (!is.finite(thresholdMarginal$thresholdX)) {
    thresholdMarginal <- .getCpUnsLocHighProbabilityReference(
      dataMod = dataMod,
      probCol = probCol,
      fraction = 0.85,
      derivativeInfo = thresholdMarginal$info,
      shapeApplied = shapeApplied
    )
  }

  # Dominance is measured before global filtering on the density comparison
  # associated with the selected probability fit. Under shape enforcement,
  # both that density and the fitted probability curve use only observations
  # to the right of the shape threshold. The adjusted tailgate is the earliest
  # point at which dominance may begin; the derivative reference still
  # controls the right-tail collapse.
  dominanceLowerBoundX <- if (shapeApplied) {
    candidate <- c(shapeLowerBoundX, shapeTailgateX)
    candidate <- candidate[is.finite(candidate)]
    if (length(candidate) == 0L) NA_real_ else max(candidate)
  } else {
    NA_real_
  }
  dominanceMarginal <- .getCpUnsLocMarginalDominanceStart(
    density = attr(dataMod, "locDensityComparison"),
    startX = thresholdMarginal$thresholdX,
    densityBw = attr(dataMod, "locDensityBw"),
    dominanceRatio = 2,
    onsetWeight = 2 / 3,
    lowerBoundX = dominanceLowerBoundX
  )
  marginalReferenceX <- .getCpUnsLocMarginalReferenceX(
    derivativeX = thresholdMarginal$thresholdX,
    dominanceX = dominanceMarginal$startX
  )
  info$marginalReference <- list(
    referenceX = marginalReferenceX,
    derivativeThreshold = thresholdMarginal$info,
    densityDominance = dominanceMarginal$info,
    shapeLowerBoundX = shapeLowerBoundX,
    dominanceLowerBoundX = dominanceLowerBoundX
  )

  if (!is.finite(marginalReferenceX)) {
    info$applied <- TRUE
    info$reason <- "no_informative_marginal_reference"
    info$thresholdClass <- "undefined"
    info$shareable <- FALSE
    return(list(
      dataMod = dataMod,
      cp = .getCpUnsLocConditionCpNonLoc(
        cpMin = cpMin,
        exTblStimNoMin = exTblStimNoMin,
        exTblUnsBias = exTblUnsBias
      ),
      info = info
    ))
  }

  # The global stage is only a soft preprocessing filter. Record the first
  # remaining expression value before marginal filtering so it remains a
  # distinct lower constraint in the marginal scan.
  global <- .getCpUnsLocFilterGlobal(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol,
    threshold = thresholdGlobal
  )
  dataMod <- global$dataMod
  info$global <- global$info

  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    info$applied <- TRUE
    info$reason <- "all_cells_removed_by_global_filter"
    return(.getCpUnsLocEmptyFilterResult(
      dataMod = dataMod,
      info = info,
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias
    ))
  }

  globalLowerBoundX <- suppressWarnings(
    min(as.numeric(.getCut(dataMod)), na.rm = TRUE)
  )

  info$antimode <- list(
    applied = FALSE,
    reason = "post_smoothing_antimode_filter_disabled",
    shapeThresholdAppliedBeforeSmoothing = shapeApplied
  )

  marginal <- .getCpUnsLocFilterMarginal(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol,
    threshold = thresholdMarginal,
    dominance = dominanceMarginal,
    globalLowerBoundX = globalLowerBoundX,
    shapeLowerBoundX = shapeLowerBoundX
  )
  dataMod <- marginal$dataMod
  info$marginal <- marginal$info

  # Retain the older diagnostic grouping while exposing the clearer stage names.
  info$lowProbLeft <- list(
    applied = isTRUE(global$info$applied) || isTRUE(marginal$info$applied),
    global = global$info,
    marginal = marginal$info
  )
  info$flatLeft <- info$lowProbLeft

  info$applied <- isTRUE(global$info$applied) ||
    isTRUE(marginal$info$applied)
  info$reason <- if (isTRUE(info$applied)) {
    "filtered_before_threshold"
  } else {
    "all_filtering_stages_kept_all_cells"
  }

  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    info$reason <- "all_cells_removed_by_post_smoothing_filtering"
    return(.getCpUnsLocEmptyFilterResult(
      dataMod = dataMod,
      info = info,
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias
    ))
  }

  list(dataMod = dataMod, cp = NULL, info = info)
}

#' Return the standard non-local result after a filter removes every cell
#' @keywords internal
.getCpUnsLocEmptyFilterResult <- function(
  dataMod,
  info,
  cpMin,
  exTblStimNoMin,
  exTblUnsBias
) {
  list(
    dataMod = dataMod,
    cp = .getCpUnsLocConditionCpNonLoc(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias
    ),
    info = info
  )
}

#' Use a high fitted-probability point when the derivative reference fails
#'
#' The first point reaching a fixed fraction of the maximum fitted probability
#' supplies a conservative right-hand reference when no derivative peak is
#' identifiable. It is calculated from whichever ordinary or shape-restricted
#' probability fit the user selected.
#' @keywords internal
.getCpUnsLocHighProbabilityReference <- function(
  dataMod,
  probCol,
  fraction = 0.85,
  derivativeInfo = NULL,
  shapeApplied = FALSE
) {
  info <- list(
    reason = "high_probability_reference_unavailable",
    stage = "marginal",
    fraction = fraction,
    derivativeFailure = derivativeInfo,
    fallback = TRUE,
    shapeThresholdApplied = isTRUE(shapeApplied)
  )
  derivTbl <- attr(dataMod, "locProbDerivTbl")
  if (
    is.data.frame(derivTbl) &&
      all(c("x", "pred") %in% names(derivTbl))
  ) {
    x <- suppressWarnings(as.numeric(derivTbl$x))
    prob <- suppressWarnings(as.numeric(derivTbl$pred))
    source <- "smoothed_probability_grid"
  } else {
    x <- suppressWarnings(as.numeric(.getCut(dataMod)))
    prob <- .getCpUnsLocProbability(dataMod, probCol)
    source <- "model_data"
  }

  keep <- is.finite(x) & is.finite(prob)
  x <- x[keep]
  prob <- pmin(1, pmax(0, prob[keep]))
  ord <- order(x)
  x <- x[ord]
  prob <- prob[ord]

  if (
    length(x) < 2L ||
      !is.finite(fraction) ||
      fraction <= 0 ||
      fraction > 1
  ) {
    return(list(thresholdX = NA_real_, info = info))
  }

  peakProb <- max(prob, na.rm = TRUE)
  targetProb <- fraction * peakProb
  reached <- which(prob >= targetProb)
  if (length(reached) == 0L || !is.finite(peakProb) || peakProb <= 0) {
    return(list(thresholdX = NA_real_, info = info))
  }

  rightIdx <- min(reached)
  thresholdX <- x[rightIdx]
  if (rightIdx > 1L && prob[rightIdx] != prob[rightIdx - 1L]) {
    leftIdx <- rightIdx - 1L
    thresholdX <- x[leftIdx] +
      (targetProb - prob[leftIdx]) *
        (x[rightIdx] - x[leftIdx]) /
        (prob[rightIdx] - prob[leftIdx])
  }

  info$reason <- "used_probability_85pct_reference"
  info$source <- source
  info$peakProb <- peakProb
  info$targetProb <- targetProb
  info$thresholdX <- thresholdX
  list(thresholdX = thresholdX, info = info)
}

# Antimode filter ------------------------------------------------------------

#' Remove values left of the right-most antimode in the dubious-response region
#' @keywords internal
.getCpUnsLocFilterAntimode <- function(
  dataMod,
  chnlSettings,
  probCol,
  threshold = NULL,
  marginalReferenceX = Inf
) {
  info <- list(applied = FALSE, reason = "antimode_filter_not_applied")
  if (!is.data.frame(dataMod) || nrow(dataMod) < 5L) {
    info$reason <- "too_few_model_values_for_antimode_filter"
    return(list(dataMod = dataMod, info = info))
  }

  if (is.null(threshold)) {
    threshold <- .getCpUnsLocStageThreshold(
      dataMod,
      chnlSettings,
      probCol,
      stage = "antimode"
    )
  }
  upper <- threshold
  info$derivativeThreshold <- upper$info
  upperX <- upper$thresholdX
  if (!is.finite(upperX)) {
    info$reason <- upper$info$reason
    return(list(dataMod = dataMod, info = info))
  }

  expr <- suppressWarnings(as.numeric(.getCut(dataMod)))
  exprForDensity <- expr[is.finite(expr)]
  if (length(exprForDensity) < 5L || length(unique(exprForDensity)) < 3L) {
    info$reason <- "too_few_values_below_antimode_upper_limit"
    return(list(dataMod = dataMod, info = info))
  }

  density <- .getCpUnsLocAntimodeDensity(
    expr = exprForDensity,
    chnlSettings = chnlSettings,
    originalBw = attr(dataMod, "locDensityBw")
  )
  if (is.null(density)) {
    info$reason <- "antimode_density_failed"
    return(list(dataMod = dataMod, info = info))
  }

  antimodes <- .getCpUnsLocAntimodes(density)
  marginalReferenceX <- suppressWarnings(
    as.numeric(marginalReferenceX)[1L]
  )
  if (!is.finite(marginalReferenceX)) {
    marginalReferenceX <- Inf
  }
  eligibleUpperX <- min(upperX, marginalReferenceX)
  eligible <- antimodes[
    is.finite(antimodes) &
      antimodes <= upperX &
      antimodes < marginalReferenceX
  ]

  info$upperX <- upperX
  info$marginalReferenceX <- marginalReferenceX
  info$eligibleUpperX <- eligibleUpperX
  info$rise <- upper$info
  info$riseThresholdX <- upperX
  info$risingFastX <- upperX
  info$antimodes <- antimodes
  info$antimodeX <- antimodes
  info$eligibleAntimodes <- eligible
  info$antimodeLeftX <- eligible
  info$nExpressionValuesForDensity <- length(exprForDensity)
  info$densityMtd <- attr(density, "locDensityMtd")
  info$densityBwType <- attr(density, "locBwType")
  info$densityBwFraction <- attr(density, "locBwFraction")
  info$densityBwBase <- attr(density, "locBwBaseSummary")
  info$densityBwUsed <- attr(density, "locBwUsedSummary")

  if (length(eligible) == 0L) {
    info$reason <- "no_antimode_in_dubious_response_region"
    return(list(dataMod = dataMod, info = info))
  }

  filterX <- max(eligible)
  keep <- is.finite(expr) & expr >= filterX
  info$filterX <- filterX
  info$nDropped <- sum(!keep)
  info$applied <- info$nDropped > 0L
  info$reason <- if (info$applied) {
    "dropped_values_left_of_rightmost_eligible_antimode"
  } else {
    "rightmost_eligible_antimode_kept_all_values"
  }

  list(dataMod = .getCpUnsLocSubsetRows(dataMod, keep), info = info)
}

#' Piecewise-constant taut-string density via cytoUtils
#'
#' Returns a numeric vector of length `n - 1` (one value per interval between
#' consecutive sorted data points) representing a piecewise-constant density
#' consistent with the gate structure returned by `cytoUtils::tautstring()`.
#' Heights are computed as a histogram within each gate region so that mode
#' regions have higher values than antimode regions, enabling downstream
#' antimode detection.
#'
#' @keywords internal
.tautStringPmden <- function(x_sorted) {
  n <- length(x_sorted)
  if (n < 3L) {
    return(list(y = rep(0.0, max(n - 1L, 0L))))
  }

  gates <- try(
    suppressWarnings(cytoUtils::tautstring(x_sorted)),
    silent = TRUE
  )
  if (inherits(gates, "try-error") || length(gates) < 2L) {
    return(list(y = rep(0.0, n - 1L)))
  }

  y <- numeric(n - 1L)

  if (length(gates) <= 2L) {
    total_width <- x_sorted[n] - x_sorted[1L]
    if (total_width > 0) y[] <- 1.0 / total_width
    return(list(y = y))
  }

  seg_counts <- tabulate(
    findInterval(x_sorted, gates, rightmost.closed = TRUE),
    nbins = length(gates) - 1L
  )
  seg_widths <- diff(gates)
  seg_density <- ifelse(
    seg_widths > 0 & seg_counts > 0,
    seg_counts / (n * seg_widths),
    0.0
  )

  x_mid <- (x_sorted[-n] + x_sorted[-1L]) / 2
  mid_seg <- findInterval(x_mid, gates, rightmost.closed = TRUE)
  mid_seg <- pmax(1L, pmin(mid_seg, length(gates) - 1L))
  y <- seg_density[mid_seg]
  list(y = y)
}

#' Fit the density used to identify antimodes
#' @keywords internal
.getCpUnsLocAntimodeDensity <- function(
  expr,
  chnlSettings,
  originalBw = NULL,
  mtd = c("taut_string", "kde")
) {
  mtd <- match.arg(mtd)
  expr <- suppressWarnings(as.numeric(expr))
  expr <- expr[is.finite(expr)]

  if (length(expr) < 5L || length(unique(expr)) < 3L) {
    return(NULL)
  }

  if (mtd == "taut_string") {
    exprSorted <- sort(expr)
    tautFit <- try(
      suppressWarnings(.tautStringPmden(exprSorted)),
      silent = TRUE
    )
    if (inherits(tautFit, "try-error")) {
      return(NULL)
    }

    y <- suppressWarnings(as.numeric(tautFit$y))
    x <- (exprSorted[-1L] + exprSorted[-length(exprSorted)]) / 2
    if (length(x) != length(y) || length(y) < 3L || all(!is.finite(y))) {
      return(NULL)
    }

    out <- list(
      x = x,
      y = y,
      fit = tautFit,
      method = "taut_string"
    )
    attr(out, "locDensityMtd") <- "taut_string"
    attr(out, "locBwType") <- "not_applicable"
    attr(out, "locBwFraction") <- NA_real_
    attr(out, "locBwBaseSummary") <- .getCpUnsLocBwSummary(numeric(0L))
    attr(out, "locBwUsedSummary") <- .getCpUnsLocBwSummary(numeric(0L))
    return(out)
  }

  bwFraction <- suppressWarnings(as.numeric(
    .getCpUnsLocSetting(chnlSettings, "locAntimodeBwFrac", 1 / 2)
  )[1])
  if (!is.finite(bwFraction) || bwFraction <= 0) {
    bwFraction <- 1 / 2
  }

  exprRange <- range(expr, na.rm = TRUE)
  if (
    is.list(originalBw) &&
      isTRUE(originalBw$adaptive) &&
      !is.null(originalBw$grid) &&
      !is.null(originalBw$sharedGrid)
  ) {
    grid <- suppressWarnings(as.numeric(originalBw$grid))
    bw <- suppressWarnings(as.numeric(originalBw$sharedGrid))
    keep <- is.finite(grid) &
      is.finite(bw) &
      bw > 0 &
      grid >= exprRange[1] &
      grid <= exprRange[2]
    grid <- grid[keep]
    bw <- bw[keep]

    if (length(grid) >= 3L && length(grid) == length(bw)) {
      out <- .getCpUnsLocDensityAdaptiveGrid(
        x = expr,
        grid = grid,
        bwGrid = bw * bwFraction,
        normalise = TRUE
      )
      if (!is.null(out)) {
        attr(out, "locDensityMtd") <- "kde"
        attr(out, "locBwType") <- "adaptive"
        attr(out, "locBwFraction") <- bwFraction
        attr(out, "locBwBaseSummary") <- .getCpUnsLocBwSummary(bw)
        attr(out, "locBwUsedSummary") <- .getCpUnsLocBwSummary(
          bw * bwFraction
        )
        return(out)
      }
    }
  }

  bw <- .getCpUnsLocAntimodeBw(expr, chnlSettings, originalBw)
  if (!is.finite(bw) || bw <= 0) {
    return(NULL)
  }

  usedBw <- bw * bwFraction
  out <- try(
    suppressWarnings(stats::density(
      expr,
      bw = usedBw,
      n = 512L,
      from = exprRange[1],
      to = exprRange[2]
    )),
    silent = TRUE
  )
  if (inherits(out, "try-error")) {
    return(NULL)
  }

  attr(out, "locDensityMtd") <- "kde"
  attr(out, "locBwType") <- "fixed"
  attr(out, "locBwFraction") <- bwFraction
  attr(out, "locBwBaseSummary") <- .getCpUnsLocBwSummary(bw)
  attr(out, "locBwUsedSummary") <- .getCpUnsLocBwSummary(usedBw)
  out
}

#' Resolve the original fixed bandwidth used by the local-FDR densities
#' @keywords internal
.getCpUnsLocAntimodeBw <- function(expr, chnlSettings, originalBw = NULL) {
  if (!is.list(originalBw)) {
    bw <- .getCpUnsLocFiniteMin(originalBw, positive = TRUE)
    if (is.finite(bw)) {
      return(bw)
    }
  }

  for (name in c("bw", "bwCluster", "bwFallback")) {
    bw <- .getCpUnsLocFiniteMin(
      .getCpUnsLocSetting(chnlSettings, name, NA_real_),
      positive = TRUE
    )
    if (is.finite(bw)) {
      return(bw)
    }
  }

  bw <- try(suppressWarnings(ks::hpi(x = expr)), silent = TRUE)
  if (inherits(bw, "try-error")) {
    return(NA_real_)
  }
  .getCpUnsLocFiniteMin(bw, positive = TRUE)
}

#' Summarise a scalar or adaptive bandwidth
#' @keywords internal
.getCpUnsLocBwSummary <- function(bw) {
  bw <- suppressWarnings(as.numeric(bw))
  bw <- bw[is.finite(bw) & bw > 0]
  if (length(bw) == 0L) {
    return(c(min = NA_real_, median = NA_real_, max = NA_real_))
  }
  c(min = min(bw), median = stats::median(bw), max = max(bw))
}

#' Locate all antimodes in a density object
#' @keywords internal
.getCpUnsLocAntimodes <- function(density) {
  x <- suppressWarnings(as.numeric(density$x))
  y <- suppressWarnings(as.numeric(density$y))
  if (length(x) != length(y) || length(y) < 3L || all(!is.finite(y))) {
    return(numeric(0L))
  }

  if (identical(density$method, "taut_string")) {
    return(.getCpUnsLocPiecewiseConstantAntimodes(x, y))
  }

  y[!is.finite(y)] <- Inf
  sort(unique(x[.getLocalMinimaIdx(y)]))
}

#' Locate antimodes in a piecewise-constant taut-string density
#' @keywords internal
.getCpUnsLocPiecewiseConstantAntimodes <- function(x, y) {
  finite <- is.finite(x) & is.finite(y)
  x <- x[finite]
  y <- y[finite]
  if (length(x) != length(y) || length(y) < 3L) {
    return(numeric(0L))
  }

  runId <- cumsum(c(
    TRUE,
    !dplyr::near(y[-1L], y[-length(y)])
  ))
  runs <- split(seq_along(y), runId)
  runY <- vapply(runs, function(i) y[i[[1L]]], numeric(1L))
  if (length(runY) < 3L) {
    return(numeric(0L))
  }

  runLeft <- vapply(runs, function(i) min(x[i]), numeric(1L))
  runRight <- vapply(runs, function(i) max(x[i]), numeric(1L))
  internal <- seq.int(2L, length(runY) - 1L)
  minima <- internal[
    runY[internal] < runY[internal - 1L] &
      runY[internal] < runY[internal + 1L]
  ]

  sort(unique((runLeft[minima] + runRight[minima]) / 2))
}

# Global filter --------------------------------------------------------------

#' Apply the global derivative threshold
#' @keywords internal
.getCpUnsLocFilterGlobal <- function(
  dataMod,
  chnlSettings,
  probCol,
  threshold = NULL
) {
  if (is.null(threshold)) {
    threshold <- .getCpUnsLocStageThreshold(
      dataMod,
      chnlSettings,
      probCol,
      stage = "global"
    )
  }
  info <- list(
    applied = FALSE,
    reason = threshold$info$reason,
    thresholdX = threshold$thresholdX,
    derivativeThreshold = threshold$info
  )

  if (!is.finite(threshold$thresholdX)) {
    info$reason <- "global_derivative_threshold_undefined"
    return(list(dataMod = dataMod, info = info))
  }

  x <- suppressWarnings(as.numeric(.getCut(dataMod)))
  keep <- is.finite(x) & x >= threshold$thresholdX
  info$nDropped <- sum(!keep)
  info$applied <- info$nDropped > 0L
  info$reason <- if (info$applied) {
    "dropped_values_left_of_global_derivative_threshold"
  } else {
    "global_derivative_threshold_kept_all_values"
  }

  list(dataMod = .getCpUnsLocSubsetRows(dataMod, keep), info = info)
}

# Marginal filter ------------------------------------------------------------

#' Find the marginal reference threshold and scan bins to its left
#' @keywords internal
.getCpUnsLocFilterMarginal <- function(
  dataMod,
  chnlSettings,
  probCol,
  antimodeX = NULL,
  threshold = NULL,
  dominance = NULL,
  globalLowerBoundX = NA_real_,
  shapeLowerBoundX = NA_real_
) {
  if (is.null(threshold)) {
    threshold <- .getCpUnsLocStageThreshold(
      dataMod,
      chnlSettings,
      probCol,
      stage = "marginal"
    )
  }

  if (is.null(dominance)) {
    dominance <- .getCpUnsLocMarginalDominanceStart(
      density = attr(dataMod, "locDensityComparison"),
      startX = threshold$thresholdX,
      densityBw = attr(dataMod, "locDensityBw"),
      dominanceRatio = 2,
      onsetWeight = 2 / 3,
      lowerBoundX = shapeLowerBoundX
    )
  }

  derivativeX <- suppressWarnings(as.numeric(threshold$thresholdX)[1L])
  dominanceX <- suppressWarnings(as.numeric(dominance$startX)[1L])
  referenceX <- .getCpUnsLocMarginalReferenceX(
    derivativeX = derivativeX,
    dominanceX = dominanceX
  )

  if (!is.finite(referenceX)) {
    return(list(
      dataMod = dataMod,
      info = list(
        applied = FALSE,
        reason = "marginal_reference_undefined",
        thresholdX = NA_real_,
        thresholdClass = "undefined",
        shareable = FALSE,
        derivativeThreshold = threshold$info,
        densityDominance = dominance$info
      )
    ))
  }

  referenceBasis <- if (
    is.finite(dominanceX) &&
      (!is.finite(derivativeX) || dominanceX < derivativeX)
  ) {
    "density_dominance_rise"
  } else {
    "marginal_probability_derivative"
  }

  shapeLowerBoundX <- suppressWarnings(as.numeric(shapeLowerBoundX)[1L])
  densityLowerBound <- list(
    lowerBoundX = shapeLowerBoundX,
    info = list(
      applied = is.finite(shapeLowerBoundX),
      reason = if (is.finite(shapeLowerBoundX)) {
        "enforced_prefit_shape_threshold"
      } else {
        "shape_threshold_not_enforced"
      },
      lowerBoundX = shapeLowerBoundX
    )
  )
  globalLowerBoundX <- suppressWarnings(
    as.numeric(globalLowerBoundX)[1L]
  )
  lowerBoundCandidates <- c(globalLowerBoundX, shapeLowerBoundX)
  lowerBoundCandidates <- lowerBoundCandidates[
    is.finite(lowerBoundCandidates)
  ]
  lowerBoundX <- if (length(lowerBoundCandidates) == 0L) {
    NA_real_
  } else {
    max(lowerBoundCandidates)
  }

  out <- .getCpUnsLocFilterMarginalBins(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol,
    startX = referenceX,
    lowerBoundX = lowerBoundX
  )
  out$info$thresholdX <- referenceX
  out$info$purityStartX <- referenceX
  out$info$referenceBasis <- referenceBasis
  out$info$thresholdClass <- "informative"
  out$info$shareable <- TRUE
  out$info$derivativeThreshold <- threshold$info
  out$info$densityDominance <- dominance$info
  out$info$derivativeThresholdX <- derivativeX
  out$info$dominanceThresholdX <- dominanceX
  out$info$globalLowerBoundX <- globalLowerBoundX
  out$info$shapeLowerBoundX <- shapeLowerBoundX
  out$info$densityLowerBoundX <- lowerBoundX
  out$info$densityLowerBound <- densityLowerBound$info
  out
}

#' Take the lower informative marginal reference
#' @keywords internal
.getCpUnsLocMarginalReferenceX <- function(derivativeX, dominanceX) {
  candidate <- suppressWarnings(as.numeric(c(derivativeX, dominanceX)))
  candidate <- candidate[is.finite(candidate)]
  if (length(candidate) == 0L) {
    return(NA_real_)
  }
  min(candidate)
}

#' Adjust and validate the tailgate floor against the marginal reference
#' @keywords internal
.getCpUnsLocSelectTailgateLowerBound <- function(
  rawTailgateX,
  windowWidth,
  referenceX,
  adjustmentFraction = 1 / 4
) {
  rawTailgateX <- suppressWarnings(as.numeric(rawTailgateX)[1L])
  windowWidth <- suppressWarnings(as.numeric(windowWidth)[1L])
  referenceX <- suppressWarnings(as.numeric(referenceX)[1L])
  adjustmentFraction <- suppressWarnings(
    as.numeric(adjustmentFraction)[1L]
  )
  info <- list(
    lowerBoundXRaw = rawTailgateX,
    windowWidth = windowWidth,
    marginalReferenceX = referenceX,
    adjustmentFraction = adjustmentFraction,
    adjustmentMethod = "none_missing_window_width",
    adjustmentRolledBack = FALSE
  )

  if (!is.finite(rawTailgateX) || !is.finite(referenceX)) {
    info$selectionReason <- "invalid_tailgate_or_marginal_reference"
    info$lowerBoundX <- NA_real_
    return(list(lowerBoundX = NA_real_, info = info))
  }

  adjustedTailgateX <- rawTailgateX
  if (
    is.finite(windowWidth) &&
      windowWidth > 0 &&
      is.finite(adjustmentFraction) &&
      adjustmentFraction >= 0
  ) {
    adjustedTailgateX <-
      rawTailgateX + adjustmentFraction * windowWidth
    info$adjustmentMethod <- "window_width"
  }
  info$lowerBoundXAdjusted <- adjustedTailgateX

  # The adjustment is only a conservative nudge. If it would overtake the
  # informative right-hand reference, retain the raw tailgate instead.
  if (adjustedTailgateX >= referenceX && rawTailgateX < referenceX) {
    selectedTailgateX <- rawTailgateX
    info$adjustmentRolledBack <- TRUE
    info$selectionReason <- "adjustment_crossed_marginal_reference"
  } else {
    selectedTailgateX <- adjustedTailgateX
    info$selectionReason <- "selected_adjusted_tailgate"
  }

  # Shape evidence to the right of the response reference is ignored rather
  # than allowed to replace or invalidate that informative reference.
  if (selectedTailgateX >= referenceX) {
    info$selectionReason <- "tailgate_not_below_marginal_reference"
    info$ignoredLowerBoundX <- selectedTailgateX
    selectedTailgateX <- NA_real_
  }
  info$lowerBoundX <- selectedTailgateX

  list(lowerBoundX = selectedTailgateX, info = info)
}

#' Identify a smoothed left-to-right rise in density dominance
#'
#' Stimulated and unstimulated densities are smoothed by Gaussian kernel
#' regression using half the fixed density bandwidth. When a derivative
#' threshold exists, densities to its right are first collapsed to their
#' respective means. A dominance threshold requires an observed transition
#' from non-dominance to dominance. Its location is two-thirds of the onset
#' location plus one-third of the subsequent dominance-score peak location.
#' @keywords internal
.getCpUnsLocMarginalDominanceStart <- function(
  density,
  startX = NA_real_,
  densityBw = NULL,
  dominanceRatio = 2,
  onsetWeight = 2 / 3,
  lowerBoundX = NA_real_
) {
  info <- list(
    applied = FALSE,
    reason = "density_dominance_rise_unavailable",
    derivativeStartX = startX,
    dominanceRatio = dominanceRatio,
    onsetWeight = onsetWeight,
    lowerBoundX = lowerBoundX
  )

  if (
    !is.data.frame(density) ||
      !all(c("x", "stim", "unstim") %in% names(density)) ||
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
    info$reason <- if (is.finite(lowerBoundX)) {
      "insufficient_density_dominance_grid_right_of_lower_bound"
    } else {
      "insufficient_density_dominance_grid"
    }
    return(list(startX = NA_real_, info = info))
  }

  derivativeStartX <- suppressWarnings(as.numeric(startX)[1L])
  tailIdx <- integer(0L)
  if (is.finite(derivativeStartX)) {
    tailIdx <- which(x >= derivativeStartX)
    if (length(tailIdx) > 0L) {
      stim[tailIdx] <- mean(stim[tailIdx])
      unstim[tailIdx] <- mean(unstim[tailIdx])
    }
  }

  # Gaussian Nadaraya-Watson regression on each density, using half the fixed
  # bandwidth that generated the original density comparison.
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

  transitionIdx <- which(
    validScore[-length(validScore)] &
      validScore[-1L] &
      !dominant[-length(dominant)] &
      dominant[-1L]
  ) + 1L
  if (length(transitionIdx) == 0L) {
    info$reason <- "density_dominance_has_no_observed_left_rise"
    info$tailCollapsed <- length(tailIdx) > 0L
    info$smoothBw <- smoothBw
    info$scoreThreshold <- scoreThreshold
    info$scoreTbl <- data.frame(
      x = x,
      stim = stimSmooth,
      unstim = unstimSmooth,
      score = score,
      dominant = dominant
    )
    return(list(startX = NA_real_, info = info))
  }

  onsetIdx <- transitionIdx[[1L]]
  leftIdx <- onsetIdx - 1L
  scoreChange <- score[onsetIdx] - score[leftIdx]
  onsetX <- if (!is.finite(scoreChange) || scoreChange == 0) {
    x[onsetIdx]
  } else {
    x[leftIdx] +
      (scoreThreshold - score[leftIdx]) *
        (x[onsetIdx] - x[leftIdx]) / scoreChange
  }

  peakCandidate <- seq.int(onsetIdx, length(x))
  peakCandidate <- peakCandidate[validScore[peakCandidate]]
  if (length(peakCandidate) == 0L) {
    info$reason <- "density_dominance_peak_unavailable_after_onset"
    return(list(startX = NA_real_, info = info))
  }
  peakIdx <- peakCandidate[[which.max(score[peakCandidate])]]
  peakX <- x[peakIdx]
  dominanceStartX <- onsetWeight * onsetX + (1 - onsetWeight) * peakX

  info$applied <- is.finite(dominanceStartX)
  info$reason <- "identified_smoothed_density_dominance_rise"
  info$tailCollapsed <- length(tailIdx) > 0L
  info$tailCollapseStartX <- if (length(tailIdx) > 0L) {
    x[min(tailIdx)]
  } else {
    NA_real_
  }
  info$smoothBw <- smoothBw
  info$scoreThreshold <- scoreThreshold
  info$onsetX <- onsetX
  info$peakX <- peakX
  info$peakScore <- score[peakIdx]
  info$dominanceStartX <- dominanceStartX
  info$scoreTbl <- data.frame(
    x = x,
    stim = stimSmooth,
    unstim = unstimSmooth,
    score = score,
    dominant = dominant
  )

  list(startX = dominanceStartX, info = info)
}

#' Bound the marginal scan by the stimulated peak's descending shoulder
#'
#' This does not itself define a filtering threshold. When no antimode was
#' identified, it prevents the marginal scan from considering bins below the
#' point where the negative stimulated-density derivative has flattened to the
#' requested fraction of its most negative value on the peak's right shoulder.
#'
#' @keywords internal
.getCpUnsLocMarginalDensityLowerBound <- function(
  density,
  peakX,
  fraction = 1 / 200
) {
  info <- list(
    applied = FALSE,
    reason = "stimulated_density_lower_bound_undefined",
    fraction = fraction,
    peakX = suppressWarnings(as.numeric(peakX)[1L])
  )

  if (
    !is.data.frame(density) ||
      !all(c("x", "y") %in% names(density)) ||
      !is.finite(info$peakX) ||
      !is.finite(fraction) ||
      fraction <= 0 ||
      fraction >= 1
  ) {
    info$reason <- "invalid_stimulated_density_or_peak"
    return(list(lowerBoundX = NA_real_, info = info))
  }

  x <- suppressWarnings(as.numeric(density$x))
  y <- suppressWarnings(as.numeric(density$y))
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- pmax(0, y[keep])
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  if (length(x) < 3L || anyDuplicated(x) || diff(range(x)) <= 0) {
    info$reason <- "insufficient_stimulated_density_grid"
    return(list(lowerBoundX = NA_real_, info = info))
  }

  deriv <- .getCpUnsLocDensityDerivative(x, y)
  peakIdx <- which.min(abs(x - info$peakX))
  afterPeak <- seq.int(peakIdx, length(x))
  negativeIdx <- afterPeak[deriv[afterPeak] < 0]
  if (length(negativeIdx) == 0L) {
    info$reason <- "no_negative_derivative_right_of_stimulated_peak"
    return(list(lowerBoundX = NA_real_, info = info))
  }

  shoulderStart <- min(negativeIdx)
  afterShoulderStart <- seq.int(shoulderStart, length(x))
  nonnegativeIdx <- afterShoulderStart[
    deriv[afterShoulderStart] >= 0
  ]
  shoulderEnd <- if (length(nonnegativeIdx) == 0L) {
    length(x)
  } else {
    min(nonnegativeIdx)
  }

  shoulderNegativeIdx <- seq.int(
    shoulderStart,
    max(shoulderStart, shoulderEnd - 1L)
  )
  steepestIdx <- shoulderNegativeIdx[
    which.min(deriv[shoulderNegativeIdx])
  ]
  steepestDeriv <- deriv[steepestIdx]
  targetDeriv <- fraction * steepestDeriv

  laterIdx <- if (steepestIdx < shoulderEnd) {
    seq.int(steepestIdx + 1L, shoulderEnd)
  } else {
    integer(0L)
  }
  crossingIdx <- laterIdx[deriv[laterIdx] >= targetDeriv]
  if (length(crossingIdx) == 0L) {
    info$reason <- "density_derivative_did_not_flatten_to_target"
    info$shoulderPeakX <- x[steepestIdx]
    info$shoulderPeakDeriv <- steepestDeriv
    info$targetDeriv <- targetDeriv
    return(list(lowerBoundX = NA_real_, info = info))
  }

  rightIdx <- min(crossingIdx)
  leftIdx <- rightIdx - 1L
  derivChange <- deriv[rightIdx] - deriv[leftIdx]
  lowerBoundX <- if (!is.finite(derivChange) || derivChange == 0) {
    x[rightIdx]
  } else {
    x[leftIdx] +
      (targetDeriv - deriv[leftIdx]) *
        (x[rightIdx] - x[leftIdx]) /
        derivChange
  }

  info$applied <- is.finite(lowerBoundX)
  info$reason <- if (info$applied) {
    "identified_stimulated_density_shoulder_lower_bound"
  } else {
    "stimulated_density_shoulder_lower_bound_undefined"
  }
  info$shoulderPeakX <- x[steepestIdx]
  info$shoulderPeakDeriv <- steepestDeriv
  info$targetDeriv <- targetDeriv
  info$lowerBoundX <- lowerBoundX

  list(lowerBoundX = lowerBoundX, info = info)
}

#' Calculate the density derivative by finite differences
#' @keywords internal
.getCpUnsLocDensityDerivative <- function(x, y) {
  n <- length(x)
  deriv <- numeric(n)
  deriv[1L] <- (y[2L] - y[1L]) / (x[2L] - x[1L])
  deriv[n] <- (y[n] - y[n - 1L]) / (x[n] - x[n - 1L])
  deriv[2L:(n - 1L)] <-
    (y[3L:n] - y[1L:(n - 2L)]) /
    (x[3L:n] - x[1L:(n - 2L)])
  deriv
}

#' Apply the appendix marginal-bin acceptance rule
#' @keywords internal
.getCpUnsLocFilterMarginalBins <- function(
  dataMod,
  chnlSettings,
  probCol,
  startX,
  lowerBoundX = NA_real_
) {
  info <- list(
    applied = FALSE,
    reason = "marginal_filter_not_run",
    startX = startX,
    lowerBoundX = lowerBoundX
  )
  if (!is.data.frame(dataMod) || nrow(dataMod) < 4L || !is.finite(startX)) {
    info$reason <- "insufficient_data_for_marginal_filter"
    return(list(dataMod = dataMod, info = info))
  }

  dataMod <- .getCpUnsLocSubsetRows(dataMod, order(.getCut(dataMod)))
  x <- suppressWarnings(as.numeric(.getCut(dataMod)))
  dm <- tibble::tibble(
    x = x,
    prob = .getCpUnsLocProbability(dataMod, probCol)
  ) |>
    dplyr::filter(is.finite(.data$x), is.finite(.data$prob)) |>
    dplyr::arrange(.data$x)

  if (nrow(dm) < 4L || diff(range(dm$x)) <= 0) {
    info$reason <- "insufficient_finite_data_for_marginal_filter"
    return(list(dataMod = dataMod, info = info))
  }

  grid <- .getCpUnsLocMarginalBreaks(dataMod, startX, nBin = 50L)
  if (is.null(grid)) {
    info$reason <- "insufficient_bins_for_marginal_filter"
    return(list(dataMod = dataMod, info = info))
  }

  breaks <- grid$breaks
  refIndex <- grid$refIndex
  leftBins <- rev(seq_len(refIndex - 1L))
  info$nLeftBinsBeforeDensityBound <- length(leftBins)
  info$nLeftBinsExcludedByDensityBound <- 0L

  if (is.finite(lowerBoundX) && length(leftBins) > 0L) {
    eligibleBin <- breaks[leftBins] >= lowerBoundX
    info$nLeftBinsExcludedByDensityBound <- sum(!eligibleBin)
    leftBins <- leftBins[eligibleBin]

    if (length(leftBins) == 0L) {
      keep <- is.finite(x) & x >= startX
      info$nDropped <- sum(!keep)
      info$applied <- info$nDropped > 0L
      info$reason <- if (info$applied) {
        "used_marginal_reference_because_density_bound_excluded_all_bins"
      } else {
        "marginal_reference_kept_all_values_after_density_bound"
      }
      info$stopReason <- "density_lower_bound_excluded_all_left_bins"
      info$finalStartX <- startX
      info$gridShift <- grid$shift
      info$gridSpacing <- grid$spacing
      info$refIndex <- refIndex
      info$scanTbl <- tibble::tibble()
      return(list(
        dataMod = .getCpUnsLocSubsetRows(dataMod, keep),
        info = info
      ))
    }
  }

  refNBin <- length(breaks) - refIndex
  right <- dm$x >= startX
  refNCell <- sum(right)
  refExpectedResp <- sum(dm$prob[right])

  if (refNBin < 1L || refNCell < 1L || !is.finite(refExpectedResp)) {
    info$reason <- "invalid_right_reference_region"
    return(list(dataMod = dataMod, info = info))
  }

  refCellsPerBin <- refNCell / refNBin
  refPurity <- refExpectedResp / refNCell
  if (!is.finite(refPurity) || refPurity <= 0) {
    info$reason <- "invalid_right_reference_purity"
    return(list(dataMod = dataMod, info = info))
  }

  cellBinRatio <- suppressWarnings(as.numeric(
    .getCpUnsLocSetting(chnlSettings, "locMarginalCellBinRatio", 1.25)
  )[1])
  if (!is.finite(cellBinRatio) || cellBinRatio <= 0) {
    cellBinRatio <- 1.25
  }
  purityRel <- .getCpUnsLocUnitValue(
    .getCpUnsLocSetting(chnlSettings, "locMarginalPurityRel", 0.85),
    0.85,
    allowZero = TRUE
  )

  maxCells <- cellBinRatio * refCellsPerBin
  minPurity <- purityRel * refPurity
  if (length(leftBins) == 0L) {
    info$reason <- "no_bins_left_of_marginal_reference"
    return(list(dataMod = dataMod, info = info))
  }

  currentCut <- startX
  scan <- vector("list", length(leftBins))
  stopReason <- "accepted_all_left_bins"
  consecutiveRejections <- 0L
  lastAcceptedJ <- 0L

  for (j in seq_along(leftBins)) {
    i <- leftBins[[j]]
    inBin <- dm$x >= breaks[[i]] & dm$x < breaks[[i + 1L]]
    nCell <- sum(inBin)
    expectedResp <- sum(dm$prob[inBin])
    purity <- if (nCell == 0L) NA_real_ else expectedResp / nCell
    accepted <- nCell == 0L ||
      (is.finite(purity) && nCell <= maxCells && purity >= minPurity)

    scan[[j]] <- tibble::tibble(
      left = breaks[[i]],
      right = breaks[[i + 1L]],
      nCell = nCell,
      expectedResp = expectedResp,
      purity = purity,
      refCellsPerBin = refCellsPerBin,
      refPurity = refPurity,
      maxCells = maxCells,
      minPurity = minPurity,
      accepted = accepted
    )

    if (accepted) {
      # An accepted bin also retains either one or two immediately preceding
      # rejected bins, because the cutoff moves to this accepted bin's left edge.
      consecutiveRejections <- 0L
      lastAcceptedJ <- j
      currentCut <- breaks[i]
    } else {
      consecutiveRejections <- consecutiveRejections + 1L

      # Three consecutive rejected bins terminate the scan. None of those three
      # bins is retained, because currentCut still marks the last accepted bin.
      if (consecutiveRejections >= 3L) {
        stopReason <- "three_consecutive_rejections"
        break
      }
    }
  }

  keep <- is.finite(x) & x >= currentCut
  info$applied <- sum(!keep) > 0L
  info$reason <- if (info$applied) {
    "dropped_values_left_of_marginal_boundary"
  } else {
    "marginal_filter_kept_all_values"
  }
  info$stopReason <- stopReason
  info$finalStartX <- currentCut
  info$gridShift <- grid$shift
  info$gridSpacing <- grid$spacing
  info$refIndex <- refIndex
  info$refNBin <- refNBin
  info$refNCell <- refNCell
  info$refExpectedResp <- refExpectedResp
  info$refCellsPerBin <- refCellsPerBin
  info$refPurity <- refPurity
  info$cellBinRatio <- cellBinRatio
  info$purityRel <- purityRel
  info$maxLeftBinCells <- maxCells
  info$minLeftBinPurity <- minPurity
  info$scanTbl <- dplyr::bind_rows(scan)
  if (nrow(info$scanTbl) > 0L) {
    info$scanTbl$retained <- seq_len(nrow(info$scanTbl)) <= lastAcceptedJ
  }

  list(dataMod = .getCpUnsLocSubsetRows(dataMod, keep), info = info)
}

#' Shift the original equal-width grid so one breakpoint equals x_ref
#' @keywords internal
.getCpUnsLocMarginalBreaks <- function(dataMod, startX, nBin = NULL) {
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
      length.out = if (is.null(nBin)) 512L else nBin
    )
  }

  binVec <- sort(unique(suppressWarnings(as.numeric(binVec))))
  binVec <- binVec[is.finite(binVec)]
  if (!is.null(nBin) && length(binVec) > 1L && diff(range(binVec)) > 0) {
    binVec <- seq(
      min(binVec),
      max(binVec),
      length.out = nBin
    )
  }
  if (length(binVec) < 2L || !is.finite(startX)) {
    return(NULL)
  }

  spacing <- stats::median(diff(binVec))
  if (!is.finite(spacing) || spacing <= 0) {
    return(NULL)
  }

  tolerance <- max(
    sqrt(.Machine$double.eps) * max(abs(c(binVec, startX)), 1),
    spacing * 1e-10
  )
  anchors <- which(binVec <= startX + tolerance)
  if (length(anchors) == 0L) {
    return(NULL)
  }

  refIndex <- max(anchors)
  shift <- max(0, startX - binVec[[refIndex]])
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
