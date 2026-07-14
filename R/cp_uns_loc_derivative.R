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

#' Read the first available setting from a precedence-ordered vector
#' @keywords internal
.getCpUnsLocFirstSetting <- function(chnlSettings, names, default = NULL) {
  for (name in names) {
    value <- .getCpUnsLocSetting(chnlSettings, name)
    if (!is.null(value)) {
      return(value)
    }
  }
  default
}

#' Validate a value constrained to the unit interval
#' @keywords internal
.getCpUnsLocUnitValue <- function(value, default, allowZero = FALSE) {
  value <- suppressWarnings(as.numeric(value)[1])
  lowerOk <- if (isTRUE(allowZero)) value >= 0 else value > 0
  if (!is.finite(value) || !lowerOk || value > 1) {
    default
  } else {
    value
  }
}

#' Get appendix parameters (alpha, omega, psi) for one filtering stage
#' @keywords internal
.getCpUnsLocDerivParams <- function(chnlSettings, stage) {
  stage <- match.arg(stage, c("antimode", "global", "marginal"))
  stageTitle <- switch(
    stage,
    antimode = "Antimode",
    global = "Global",
    marginal = "Marginal"
  )
  defaults <- switch(
    stage,
    antimode = c(alpha = 2 / 3, omega = 0.15, psi = -0.1),
    global = c(alpha = 0.05, omega = 0.15, psi = 0.2),
    marginal = c(alpha = 0.50, omega = 0.15, psi = -0.1)
  )

  alpha <- .getCpUnsLocFirstSetting(
    chnlSettings,
    c(
      paste0("loc", stageTitle, "DerivAlpha"),
      "locDerivAlpha",
      paste0("loc", stageTitle, "DerivPeakMinRel"),
      "locDerivPeakMinRel"
    ),
    defaults[["alpha"]]
  )
  omega <- .getCpUnsLocFirstSetting(
    chnlSettings,
    c(
      paste0("loc", stageTitle, "DerivOmega"),
      "locDerivOmega",
      paste0("loc", stageTitle, "DerivPeakProbMin"),
      "locDerivPeakProbMin"
    ),
    defaults[["omega"]]
  )
  psi <- .getCpUnsLocFirstSetting(
    chnlSettings,
    c(
      paste0("loc", stageTitle, "DerivPsi"),
      "locDerivPsi",
      paste0("loc", stageTitle, "DerivRiseFrac")
    ),
    defaults[["psi"]]
  )

  list(
    alpha = .getCpUnsLocUnitValue(alpha, defaults[["alpha"]]),
    omega = .getCpUnsLocUnitValue(
      omega,
      defaults[["omega"]],
      allowZero = TRUE
    ),
    psi = .getCpUnsLocUnitValue(psi, defaults[["psi"]])
  )
}

#' Get the optional minimum probability at the final derivative threshold
#' @keywords internal
.getCpUnsLocThresholdProbMin <- function(chnlSettings, stage) {
  stage <- match.arg(stage, c("antimode", "global", "marginal"))
  stageTitle <- switch(
    stage,
    antimode = "Antimode",
    global = "Global",
    marginal = "Marginal"
  )
  common <- .getCpUnsLocUnitValue(
    .getCpUnsLocSetting(chnlSettings, "locDerivRiseProbMin", 0),
    0,
    allowZero = TRUE
  )
  .getCpUnsLocUnitValue(
    .getCpUnsLocSetting(
      chnlSettings,
      paste0("loc", stageTitle, "DerivRiseProbMin"),
      common
    ),
    common,
    allowZero = TRUE
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

#' Extract finite fitted probabilities on [0, 1]
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
    "locStimPeakX"
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
.getCpUnsLocFiniteMin <- function(x, positive = FALSE) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (isTRUE(positive)) {
    x <- x[x > 0]
  }
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
  probTbl <- tibble::tibble(
    x = x,
    prob = .getCpUnsLocProbability(dataMod, probCol)
  ) |>
    dplyr::filter(is.finite(.data$x), is.finite(.data$prob)) |>
    dplyr::group_by(.data$x) |>
    dplyr::summarise(prob = mean(.data$prob), .groups = "drop") |>
    dplyr::arrange(.data$x)

  if (nrow(probTbl) < 4L || diff(range(probTbl$x)) <= 0) {
    return(NULL)
  }

  dx <- diff(probTbl$x)
  deriv <- pmax(0, diff(probTbl$prob) / dx)
  deriv[!is.finite(deriv)] <- 0

  tibble::tibble(
    x = (head(probTbl$x, -1L) + tail(probTbl$x, -1L)) / 2,
    prob = (head(probTbl$prob, -1L) + tail(probTbl$prob, -1L)) / 2,
    deriv = deriv,
    source = "model_probability_finite_difference"
  )
}

#' Select the left-most derivative peak meeting alpha
#'
#' Flat-topped peaks are represented by the left-most point of the plateau.
#' @keywords internal
.getCpUnsLocDerivPeak <- function(x, prob, deriv, alpha = 0.75) {
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

  # Retain the established fallback when the derivative has no internal peak.
  usedGlobalFallback <- length(peakIndex) == 0L
  if (usedGlobalFallback) {
    peakIndex <- which.max(peakData$deriv)
  }

  maxPeak <- max(peakData$deriv[peakIndex], na.rm = TRUE)
  eligible <- peakIndex[peakData$deriv[peakIndex] >= alpha * maxPeak]

  info$alpha <- alpha
  info$peakMinRel <- alpha
  info$globalMaxDeriv <- max(peakData$deriv, na.rm = TRUE)
  info$maxPeakDeriv <- maxPeak
  info$usedGlobalMaximumFallback <- usedGlobalFallback
  info$peakSummary <- data.frame(
    index = peakIndex,
    idx = peakIndex,
    x = peakData$x[peakIndex],
    prob = peakData$prob[peakIndex],
    deriv = peakData$deriv[peakIndex],
    relativeHeight = peakData$deriv[peakIndex] / maxPeak,
    relToMaxPeak = peakData$deriv[peakIndex] / maxPeak,
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
  thresholdProbMin = 0
) {
  peak <- .getCpUnsLocDerivPeak(x, prob, deriv, alpha)
  info <- peak$info
  if (is.null(peak$data) || !is.finite(peak$index)) {
    return(list(thresholdX = NA_real_, info = info))
  }

  omega <- .getCpUnsLocUnitValue(omega, 0.15, allowZero = TRUE)
  psi <- .getCpUnsLocUnitValue(psi, 0.75)
  thresholdProbMin <- .getCpUnsLocUnitValue(
    thresholdProbMin,
    0,
    allowZero = TRUE
  )

  iPeak <- peak$index
  peakHeight <- peak$data$deriv[iPeak]
  riseHeight <- abs(psi) * peakHeight

  info$omega <- omega
  info$peakProbMin <- omega
  info$psi <- psi
  info$riseFrac <- psi
  info$riseHeight <- riseHeight
  info$thresholdProbMin <- thresholdProbMin

  if (peak$data$prob[iPeak] < omega) {
    candidate <- seq.int(iPeak, nrow(peak$data))
    candidate <- candidate[peak$data$prob[candidate] >= omega]
    if (length(candidate) == 0L) {
      info$reason <- "probability_never_reached_omega_after_peak"
      return(list(thresholdX = NA_real_, info = info))
    }
    if (psi < 0) {
      candidate <- candidate[peak$data$deriv[candidate] >= riseHeight]
    } else {
      # no need to apply psi here as we're forced to
      # move rightware, and positive psi seeks to move left.
      candidate <- candidate
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

  info$thresholdIdxCandidate <- iThreshold
  info$thresholdXCandidate <- peak$data$x[iThreshold]
  info$thresholdProbCandidate <- peak$data$prob[iThreshold]
  info$riseThresholdIdxCandidate <- iThreshold
  info$riseThresholdXCandidate <- peak$data$x[iThreshold]
  info$riseThresholdProbCandidate <- peak$data$prob[iThreshold]

  # Optional extra constraint on the threshold itself, independent of omega.
  later <- seq.int(iThreshold, nrow(peak$data))
  later <- later[peak$data$prob[later] >= thresholdProbMin]
  if (length(later) == 0L) {
    info$reason <- "no_later_threshold_met_probability_minimum"
    return(list(thresholdX = NA_real_, info = info))
  }

  iThresholdFinal <- min(later)
  thresholdX <- peak$data$x[iThresholdFinal]
  info$reason <- "identified_derivative_threshold"
  info$thresholdIdx <- iThresholdFinal
  info$thresholdX <- thresholdX
  info$thresholdProb <- peak$data$prob[iThresholdFinal]
  info$riseThresholdIdx <- iThresholdFinal
  info$riseThresholdX <- thresholdX
  info$riseThresholdProb <- peak$data$prob[iThresholdFinal]
  info$riseThreshold <- riseHeight
  info$risingFastIdx <- iThresholdFinal
  info$risingFastX <- thresholdX
  info$shiftedRightForProbability <- iThresholdFinal > iThreshold

  list(thresholdX = thresholdX, info = info)
}

#' Obtain a stage-specific derivative threshold from model data
#' @keywords internal
.getCpUnsLocStageThreshold <- function(
  dataMod,
  chnlSettings,
  probCol,
  stage
) {
  params <- .getCpUnsLocDerivParams(chnlSettings, stage)
  derivTbl <- .getCpUnsLocDerivTbl(dataMod, probCol)
  if (is.null(derivTbl) || nrow(derivTbl) < 3L) {
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
    thresholdProbMin = .getCpUnsLocThresholdProbMin(chnlSettings, stage)
  )
  out$info$stage <- stage
  out$info$params <- params
  out$info$derivSource <- unique(derivTbl$source)[1]
  out
}
