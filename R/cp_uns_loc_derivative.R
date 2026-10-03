# Appendix derivative thresholds for local-FDR filtering
#
# Defines the appendix parameters alpha, omega, and psi; extracts the fitted
# derivative over the current expression region; selects flat or pointed peaks;
# and obtains the stage-specific derivative threshold.

# Settings ------------------------------------------------------------------

#' Read one channel setting
#' @keywords internal
.getCpUnsLocSetting <- function(chnlSettings, name, default = NULL) {
  if (!is.null(chnlSettings) && !is.null(chnlSettings[[name]])) {
    chnlSettings[[name]]
  } else {
    default
  }
}

#' Validate a value constrained to the unit interval
#' @keywords internal
.getCpUnsLocUnitValue <- function(
  value,
  default,
  allowZero = FALSE,
  allowNeg = FALSE
) {
  value <- suppressWarnings(as.numeric(value)[1])
  zeroOk <- if (isTRUE(allowZero)) TRUE else value != 0
  negOk <- if (isTRUE(allowNeg)) TRUE else value >= 0
  if (!is.finite(value) || !zeroOk || !negOk || abs(value) > 1) {
    default
  } else {
    value
  }
}

#' Get appendix parameters (alpha, omega, psi) for one filtering stage
#' @keywords internal
.getCpUnsLocDerivParams <- function(stage) {
  switch(match.arg(stage, c("antimode", "global", "marginal")),
    antimode = list(alpha = 2 / 3, omega = 0.15, psi = -0.2),
    global = list(alpha = 0.05, omega = 0.15, psi = 0.2),
    marginal = list(alpha = 0.50, omega = 0.15, psi = -0.2)
  )
}

# Data helpers ---------------------------------------------------------------

#' Choose the fitted probability column used by the filters
#' @keywords internal
.getCpUnsLocProbabilityColumn <- function(dataMod, chnlSettings) {
  requested <- .getCpUnsLocSetting(chnlSettings, "locProbCol", "pred")
  if (requested %in% names(dataMod)) {
    requested
  } else if ("pred" %in% names(dataMod)) {
    "pred"
  } else {
    "probSmooth"
  }
}

#' Extract finite fitted probabilities on \code{[0, 1]}
#' @keywords internal
.getCpUnsLocProbability <- function(dataMod, probCol) {
  prob <- suppressWarnings(as.numeric(dataMod[[probCol]]))
  prob <- pmin(1, pmax(0, prob))
  prob[!is.finite(prob)] <- NA_real_
  prob
}

#' Subset model data while retaining attributes needed downstream
#' @keywords internal
.getCpUnsLocSubsetRows <- function(dataMod, keep) {
  attrs <- c(
    "chnlCut",
    "ind",
    "indUns",
    "binVec",
    "minProbXPos",
    "locProbDerivTbl",
    "locProbSmoothMethod",
    "locDensityBw",
    "locStimDensity",
    "locDensityComparison",
    "locPeakX",
    "locWindowWidth",
    "locShapeThresholdRequested",
    "locShapeThresholdApplied",
    "locShapeThresholdX",
    "locShapeThresholdInfo",
    "locShapeTailgateX",
    "locShapeAntimodeX",
    "locUnshapedProbCurve"
  )
  values <- stats::setNames(
    lapply(attrs, function(name) attr(dataMod, name)),
    attrs
  )
  out <- dataMod[keep, , drop = FALSE]
  for (name in attrs) {
    if (!is.null(values[[name]])) {
      attr(out, name) <- values[[name]]
    }
  }
  out
}

#' Return the smallest finite numeric value
#' @keywords internal
.getCpUnsLocFiniteMin <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

# Derivative threshold -------------------------------------------------------

#' Derivative table over the current expression region
#' @keywords internal
.getCpUnsLocDerivTbl <- function(dataMod, probCol) {
  x <- suppressWarnings(as.numeric(.getCut(dataMod)))
  xRange <- range(x, na.rm = TRUE)
  if (length(xRange) != 2L || any(!is.finite(xRange)) || diff(xRange) <= 0) {
    return(NULL)
  }

  fitted <- attr(dataMod, "locProbDerivTbl")
  if (
    identical(probCol, "pred") &&
      is.data.frame(fitted) &&
      all(c("x", "pred", "deriv") %in% names(fitted))
  ) {
    fitted <- fitted |>
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

    if (nrow(fitted) >= 3L) {
      return(fitted)
    }
  }

  # Fallback for direct calls or failed derivative storage.
  prob <- .getCpUnsLocProbability(dataMod, probCol)

  keep <- is.finite(x) & is.finite(prob)
  x <- x[keep]
  prob <- prob[keep]

  if (!anyDuplicated(x)) {
    ord <- order(x)

    probTbl <- tibble::tibble(
      x = x[ord],
      prob = prob[ord]
    )
  } else {
    probTbl <- tibble::tibble(
      x = x,
      prob = prob
    ) |>
      dplyr::group_by(.data$x) |>
      dplyr::summarise(
        prob = mean(.data$prob),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$x)
  }

  if (nrow(probTbl) < 4L || diff(range(probTbl$x)) <= 0) {
    return(NULL)
  }

  dx <- diff(probTbl$x)
  deriv <- pmax(0, diff(probTbl$prob) / dx)
  deriv[!is.finite(deriv)] <- 0

  tibble::tibble(
    x = (utils::head(probTbl$x, -1L) + utils::tail(probTbl$x, -1L)) / 2,
    prob = (utils::head(probTbl$prob, -1L) + utils::tail(probTbl$prob, -1L)) /
      2,
    deriv = deriv,
    source = "model_probability_finite_difference"
  )
}

#' Select the left-most derivative peak meeting alpha
#'
#' Flat-topped peaks are represented by the left-most point of the plateau.
#' @keywords internal
.getCpUnsLocDerivPeak <- function(
  x,
  prob,
  deriv,
  alpha = 0.75,
  leftRiseFrac = 0.15
) {
  info <- list(reason = "no_valid_derivative_peak")
  x <- suppressWarnings(as.numeric(x))
  prob <- suppressWarnings(as.numeric(prob))
  deriv <- suppressWarnings(as.numeric(deriv))

  if (length(x) != length(prob) || length(x) != length(deriv)) {
    info$reason <- "derivative_peak_input_lengths_differ"
    return(list(index = NA_integer_, data = NULL, info = info))
  }

  keep <- is.finite(x) & is.finite(prob) & is.finite(deriv)
  peakData <- data.frame(
    x = x[keep],
    prob = pmin(1, pmax(0, prob[keep])),
    deriv = pmax(0, deriv[keep])
  )
  peakData <- peakData[order(peakData$x), , drop = FALSE]

  if (nrow(peakData) < 3L) {
    info$reason <- "too_few_finite_derivative_points"
    return(list(index = NA_integer_, data = peakData, info = info))
  }

  alpha <- .getCpUnsLocUnitValue(alpha, 0.75)
  if (max(peakData$deriv, na.rm = TRUE) <= 0) {
    info$reason <- "no_positive_probability_derivative"
    return(list(index = NA_integer_, data = peakData, info = info))
  }

  runs <- rle(peakData$deriv)
  runEnd <- cumsum(runs$lengths)
  runStart <- runEnd - runs$lengths + 1L
  runValue <- runs$values

  peakRun <- rep(FALSE, length(runValue))
  if (length(runValue) >= 3L) {
    internal <- seq.int(2L, length(runValue) - 1L)
    peakRun[internal] <- runValue[internal] > runValue[internal - 1L] &
      runValue[internal] > runValue[internal + 1L]
  }
  peakIndex <- runStart[peakRun]

  # Retain the global-maximum fallback only when the maximum has an
  # observed, meaningful rising flank to its left.
  usedGlobalFallback <- length(peakIndex) == 0L
  if (usedGlobalFallback) {
    peakIndex <- which.max(peakData$deriv)
  }

  leftRiseFrac <- .getCpUnsLocUnitValue(
    leftRiseFrac,
    0.15,
    allowZero = TRUE
  )
  nPeakData <- nrow(peakData)

  minDerivBefore <- c(
    Inf,
    cummin(peakData$deriv)[seq_len(nPeakData - 1L)]
  )

  hasMeaningfulLeftRise <-
    peakIndex > 1L &
    minDerivBefore[peakIndex] <= leftRiseFrac * peakData$deriv[peakIndex]

  peakIndex <- peakIndex[hasMeaningfulLeftRise]

  info$leftRiseFrac <- leftRiseFrac
  info$usedGlobalMaximumFallback <- usedGlobalFallback

  if (length(peakIndex) == 0L) {
    info$reason <- "no_derivative_peak_with_meaningful_left_rise"
    return(list(
      index = NA_integer_,
      data = peakData,
      info = info
    ))
  }

  maxPeak <- max(peakData$deriv[peakIndex], na.rm = TRUE)
  eligible <- peakIndex[
    peakData$deriv[peakIndex] >= alpha * maxPeak
  ]

  info$alpha <- alpha
  info$globalMaxDeriv <- max(peakData$deriv, na.rm = TRUE)
  info$maxPeakDeriv <- maxPeak
  info$peakSummary <- data.frame(
    index = peakIndex,
    x = peakData$x[peakIndex],
    prob = peakData$prob[peakIndex],
    deriv = peakData$deriv[peakIndex],
    relativeHeight = peakData$deriv[peakIndex] / maxPeak,
    relToGlobal = peakData$deriv[peakIndex] / info$globalMaxDeriv,
    eligible = peakIndex %in% eligible
  )

  if (length(eligible) == 0L) {
    info$reason <- "no_derivative_peak_met_alpha"
    return(list(index = NA_integer_, data = peakData, info = info))
  }

  selected <- min(eligible)
  info$reason <- "identified_leftmost_valid_derivative_peak"
  info$peakIdx <- selected
  info$peakX <- peakData$x[selected]
  info$peakProb <- peakData$prob[selected]
  info$peakDeriv <- peakData$deriv[selected]

  list(index = selected, data = peakData, info = info)
}

#' Locate x_deriv(alpha, omega, psi)
#' @keywords internal
.getCpUnsLocDerivThreshold <- function(
  x,
  prob,
  deriv,
  alpha,
  omega,
  psi,
  capRightWidth = FALSE,
  leftRiseFrac = 0.15,
  stage
) {
  peak <- .getCpUnsLocDerivPeak(
    x = x,
    prob = prob,
    deriv = deriv,
    alpha = alpha,
    leftRiseFrac = leftRiseFrac
  )
  info <- peak$info
  if (is.na(peak$index)) {
    return(list(thresholdX = NA_real_, info = info))
  }

  omega <- .getCpUnsLocUnitValue(omega, 0.15, allowZero = TRUE)
  # `stage` is only needed (and evaluated) when psi is invalid.
  psi <- .getCpUnsLocUnitValue(
    psi,
    .getCpUnsLocDerivParams(stage)$psi,
    allowNeg = TRUE
  )

  iPeak <- peak$index
  peakHeight <- peak$data$deriv[iPeak]
  riseHeight <- abs(psi) * peakHeight

  info$omega <- omega
  info$psi <- psi
  info$riseHeight <- riseHeight

  if (peak$data$prob[iPeak] < omega) {
    candidate <- seq.int(iPeak, nrow(peak$data))
    candidate <- candidate[peak$data$prob[candidate] >= omega]

    if (length(candidate) == 0L) {
      info$reason <- "probability_never_reached_omega_after_peak"
      return(list(thresholdX = NA_real_, info = info))
    }

    if (psi < 0) {
      candidate <- candidate[
        peak$data$deriv[candidate] >= riseHeight
      ]
      if (length(candidate) == 0L) {
        info$reason <- "no_point_met_probability_and_derivative_constraints"
        return(list(thresholdX = NA_real_, info = info))
      }
    }

    iThreshold <- min(candidate)
    info$thresholdBasis <- "first_point_right_of_peak_reaching_omega"
  } else {
    if (psi > 0) {
      candidate <- seq_len(iPeak)
      candidate <- candidate[peak$data$deriv[candidate] >= riseHeight]
      if (length(candidate) == 0L) {
        info$reason <- "derivative_never_reached_psi_times_peak"
        return(list(thresholdX = NA_real_, info = info))
      }
      info$thresholdBasis <- "left_rise_to_psi_times_peak"
    } else {
      # here we find the indices
      # such that the derivative is less than or equal to riseHeight
      candidate <- seq.int(iPeak, nrow(peak$data))
      candidate <- candidate[peak$data$deriv[candidate] <= riseHeight]
      if (length(candidate) == 0L) {
        info$reason <- "derivative_never_fell_below_psi_times_peak"
        return(list(thresholdX = NA_real_, info = info))
      }
      info$thresholdBasis <- "right_fall_to_psi_times_peak"
    }
    iThreshold <- min(candidate)
  }

  info$rightFractionThresholdIdx <- iThreshold
  info$rightFractionThresholdX <- peak$data$x[iThreshold]
  info$rightWidthCapApplied <- FALSE

  if (isTRUE(capRightWidth) && psi < 0 && peak$data$prob[iPeak] >= omega) {
    widthCap <- .getCpUnsLocDerivRightWidthCap(
      peakData = peak$data,
      iPeak = iPeak,
      rightFrac = abs(psi)
    )
    info$rightWidthCap <- widthCap$info

    if (is.finite(widthCap$index) && widthCap$index < iThreshold) {
      iThreshold <- widthCap$index
      info$rightWidthCapApplied <- TRUE
      info$thresholdBasis <- "minimum_of_right_fraction_and_width_cap"
    }
  }

  thresholdX <- peak$data$x[iThreshold]
  info$reason <- "identified_derivative_threshold"
  info$thresholdIdx <- iThreshold
  info$thresholdX <- thresholdX
  info$thresholdProb <- peak$data$prob[iThreshold]

  list(thresholdX = thresholdX, info = info)
}

#' Cap a right-side derivative threshold using the selected peak's left width
#'
#' The left width is measured at half height. The multiplier is chosen so that
#' a Gaussian-shaped peak reaches `rightFrac` at the same standardised distance.
#' @keywords internal
.getCpUnsLocDerivRightWidthCap <- function(
  peakData,
  iPeak,
  rightFrac,
  leftFrac = 0.5
) {
  info <- list(
    reason = "right_width_cap_undefined",
    rightFrac = rightFrac,
    leftFrac = leftFrac
  )

  # The caller guarantees iPeak >= 2 and rightFrac in (0, 1].
  if (rightFrac >= 1) {
    return(list(index = NA_integer_, info = info))
  }

  peakHeight <- peakData$deriv[iPeak]
  leftHeight <- leftFrac * peakHeight
  leftOfPeak <- seq_len(iPeak - 1L)
  below <- leftOfPeak[peakData$deriv[leftOfPeak] < leftHeight]
  if (length(below) == 0L) {
    info$reason <- "left_half_height_crossing_unavailable"
    return(list(index = NA_integer_, info = info))
  }

  # deriv[iLeft] < leftHeight <= deriv[iRight], so the change is positive.
  iLeft <- max(below)
  iRight <- iLeft + 1L
  leftX <- peakData$x[iLeft] +
    (leftHeight - peakData$deriv[iLeft]) *
      (peakData$x[iRight] - peakData$x[iLeft]) /
      (peakData$deriv[iRight] - peakData$deriv[iLeft])

  peakX <- peakData$x[iPeak]
  leftWidth <- peakX - leftX
  if (!is.finite(leftWidth) || leftWidth <= 0) {
    info$reason <- "invalid_right_width_cap"
    return(list(index = NA_integer_, info = info))
  }

  widthRatio <- sqrt(log(1 / rightFrac) / log(1 / leftFrac))
  widthCapX <- peakX + widthRatio * leftWidth
  rightOfPeak <- seq.int(iPeak, nrow(peakData))
  capIndex <- rightOfPeak[
    which.min(abs(peakData$x[rightOfPeak] - widthCapX))
  ]

  info$reason <- "identified_gaussian_matched_right_width_cap"
  info$leftX <- leftX
  info$leftWidth <- leftWidth
  info$widthRatio <- widthRatio
  info$widthCapX <- widthCapX
  info$widthCapGridX <- peakData$x[capIndex]

  list(index = capIndex, info = info)
}

#' Obtain a stage-specific derivative threshold from model data
#' @keywords internal
.getCpUnsLocStageThreshold <- function(
  dataMod,
  chnlSettings,
  probCol,
  stage
) {
  params <- .getCpUnsLocDerivParams(stage)
  derivTbl <- .getCpUnsLocDerivTbl(dataMod, probCol)
  if (is.null(derivTbl)) {
    return(list(
      thresholdX = NA_real_,
      info = list(
        reason = "probability_derivative_unavailable",
        stage = stage,
        params = params
      )
    ))
  }

  out <- .getCpUnsLocDerivThreshold(
    x = derivTbl$x,
    prob = derivTbl$prob,
    deriv = derivTbl$deriv,
    alpha = params$alpha,
    omega = params$omega,
    psi = params$psi,
    capRightWidth = identical(stage, "marginal"),
    leftRiseFrac = if (isTRUE(attr(dataMod, "locShapeThresholdApplied"))) {
      0.5
    } else {
      0.15
    },
    stage = stage
  )
  out$info$stage <- stage
  out$info$params <- params
  out$info$derivSource <- unique(derivTbl$source)[1]
  out
}
