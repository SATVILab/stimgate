# Appendix-aligned post-smoothing filtering

# Local-FDR filtering after probability smoothing
#
# The names in this file follow the appendix terminology:
#   alpha: minimum relative height for an eligible derivative peak
#   omega: minimum response probability associated with that peak
#   psi: fraction of the selected peak used to locate the rising edge
#
# Filtering is sequential: antimode, global, then marginal. A missing global
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

  antimode <- .getCpUnsLocFilterAntimode(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol
  )
  dataMod <- antimode$dataMod
  info$antimode <- antimode$info

  if (!is.data.frame(dataMod) || nrow(dataMod) == 0L) {
    info$applied <- TRUE
    info$reason <- "all_cells_removed_by_antimode_filter"
    return(.getCpUnsLocEmptyFilterResult(
      dataMod = dataMod,
      info = info,
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias
    ))
  }

  global <- .getCpUnsLocFilterGlobal(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol
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

  marginal <- .getCpUnsLocFilterMarginal(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol
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

  info$applied <- isTRUE(antimode$info$applied) ||
    isTRUE(global$info$applied) ||
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

# Antimode filter ------------------------------------------------------------

#' Remove values left of the right-most antimode in the dubious-response region
#' @keywords internal
.getCpUnsLocFilterAntimode <- function(
  dataMod,
  chnlSettings,
  probCol,
  threshold = NULL
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
    originalBw = attr(dataMod, "locDensityBw"),
    mtd = "taut_string"
  )
  if (is.null(density)) {
    info$reason <- "antimode_density_failed"
    return(list(dataMod = dataMod, info = info))
  }

  antimodes <- .getCpUnsLocAntimodes(density)
  eligible <- antimodes[is.finite(antimodes) & antimodes <= upperX]

  info$upperX <- upperX
  info$rise <- upper$info
  info$riseThresholdX <- upperX
  info$risingFastX <- upperX
  info$antimodes <- antimodes
  info$antimodeX <- antimodes
  info$eligibleAntimodes <- eligible
  info$antimodeLeftX <- eligible
  info$nExpressionValuesForDensity <- length(exprForDensity)
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

#' Fit the antimode density using half the original bandwidth
#' @keywords internal
.getCpUnsLocAntimodeDensity <- function(
  expr,
  chnlSettings,
  originalBw = NULL,
  mtd = "taut_string"
) {
  if (mtd == "kde") {
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

    attr(out, "locBwType") <- "fixed"
    attr(out, "locBwFraction") <- bwFraction
    attr(out, "locBwBaseSummary") <- .getCpUnsLocBwSummary(bw)
    attr(out, "locBwUsedSummary") <- .getCpUnsLocBwSummary(usedBw)
    out
  } else if (mtd == "taut_string") {}
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
  y[!is.finite(y)] <- Inf
  sort(unique(x[.getLocalMinimaIdx(y)]))
}

# Global filter --------------------------------------------------------------

#' Apply the global derivative threshold
#' @keywords internal
.getCpUnsLocFilterGlobal <- function(dataMod, chnlSettings, probCol) {
  threshold <- .getCpUnsLocStageThreshold(
    dataMod,
    chnlSettings,
    probCol,
    stage = "global"
  )
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
.getCpUnsLocFilterMarginal <- function(dataMod, chnlSettings, probCol) {
  threshold <- .getCpUnsLocStageThreshold(
    dataMod,
    chnlSettings,
    probCol,
    stage = "marginal"
  )
  if (!is.finite(threshold$thresholdX)) {
    return(list(
      dataMod = dataMod,
      info = list(
        applied = FALSE,
        reason = "marginal_derivative_threshold_undefined",
        thresholdX = NA_real_,
        derivativeThreshold = threshold$info
      )
    ))
  }

  out <- .getCpUnsLocFilterMarginalBins(
    dataMod = dataMod,
    chnlSettings = chnlSettings,
    probCol = probCol,
    startX = threshold$thresholdX
  )
  out$info$thresholdX <- threshold$thresholdX
  out$info$derivativeThreshold <- threshold$info
  out
}

#' Apply the appendix marginal-bin acceptance rule
#' @keywords internal
.getCpUnsLocFilterMarginalBins <- function(
  dataMod,
  chnlSettings,
  probCol,
  startX
) {
  info <- list(
    applied = FALSE,
    reason = "marginal_filter_not_run",
    startX = startX
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
  leftBins <- rev(seq_len(refIndex - 1L))
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
        break
      }
    }

    currentCut <- breaks[[i]]
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
