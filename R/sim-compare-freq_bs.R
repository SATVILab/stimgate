#' Helpers for comparing background-subtracted frequency estimators
#'
#' The functions in this file compare the fixed-bandwidth StimGate estimate
#' against two alternative gates:
#'   1. an R translation of the fbeta histogram method
#'   2. a tailgate-style density derivative cutpoint
#'
#' @keywords internal

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
}

#' @keywords internal
.simCompareAddMissingColumns <- function(.data, cols) {
  for (nm in names(cols)) {
    if (!nm %in% names(.data)) {
      .data[[nm]] <- cols[[nm]]
    }
  }
  .data
}

#' @keywords internal
.simCompareGetTrans <- function(transformation) {
  if (exists(".simBandwidthGetTrans", mode = "function")) {
    return(.simBandwidthGetTrans(transformation))
  }
  if (exists(".simMiscGetTrans", mode = "function")) {
    return(.simMiscGetTrans(transformation))
  }

  switch(
    transformation,
    "gamma" = calc_gamma,
    "gaussian" = calc_gaussian,
    "skew" = calc_skew,
    stop("Transformation not recognised: ", transformation)
  )
}

#' @keywords internal
.simCompareSampleFromInd <- function(ind, nCondition) {
  if (exists(".simBandwidthSampleFromInd", mode = "function")) {
    return(.simBandwidthSampleFromInd(ind, nCondition))
  }
  ind_num <- suppressWarnings(as.numeric(ind))
  as.character(((ind_num - 1) %/% nCondition) + 1)
}

#' @keywords internal
.simCompareReadLocDetails <- function(pathProject, nSample, nCondition) {
  if (!exists(".simBandwidthReadLocDetails", mode = "function")) {
    stop(
      ".simBandwidthReadLocDetails() must be available before running ",
      "the StimGate comparison. Source sim-bandwidth.R first."
    )
  }
  .simBandwidthReadLocDetails(
    pathProject = pathProject,
    nSample = nSample,
    nCondition = nCondition
  )
}

#' Moving mean with leading NA values, matching fbeta.py
#'
#' @keywords internal
.simCompareMoveMean <- function(x, window) {
  x <- as.numeric(x)
  window <- as.integer(window)
  if (length(x) == 0L || window <= 1L) {
    return(x)
  }
  if (length(x) < window) {
    return(rep(NA_real_, length(x)))
  }
  xs <- cumsum(x)
  x1 <- xs[window:length(x)]
  x2 <- c(0, xs[seq_len(length(x) - window)])
  c(rep(NA_real_, window - 1L), (x1 - x2) / window)
}

#' Calculate the fbeta score from two density vectors
#'
#' This is an R translation of calculate_fscore() in fbeta.py.
#'
#' @keywords internal
.simCompareFbetaCalculateFscore <- function(
  negPdf,
  posPdf,
  beta = 0.8,
  theta = 2
) {
  negPdf <- as.numeric(negPdf)
  posPdf <- as.numeric(posPdf)
  if (length(negPdf) != length(posPdf)) {
    stop("negPdf and posPdf must have the same length.")
  }

  n <- length(negPdf)
  fpos <- ifelse(posPdf > theta * negPdf, posPdf - negPdf, 0)
  fpos[!is.finite(fpos)] <- 0

  tp <- vapply(
    seq_len(n),
    function(i) sum(fpos[i:n], na.rm = FALSE),
    numeric(1)
  )
  fn <- vapply(
    seq_len(n),
    function(i) {
      if (i == 1L) {
        return(0)
      }
      sum(fpos[seq_len(i - 1L)], na.rm = FALSE)
    },
    numeric(1)
  )
  fp <- vapply(
    seq_len(n),
    function(i) sum(negPdf[i:n], na.rm = FALSE),
    numeric(1)
  )

  precision <- tp / (tp + fp)
  precision[tp == 0] <- 0
  recall <- tp / (tp + fn)
  recall[recall == 0] <- 0

  fscores <- (1 + beta * beta) *
    (precision * recall) /
    (beta * beta * precision + recall)
  fscores[!is.finite(fscores)] <- 0
  precision[!is.finite(precision)] <- 0
  recall[!is.finite(recall)] <- 0

  list(
    fscores = fscores,
    precision = precision,
    recall = recall
  )
}

#' @keywords internal
.simCompareHistDensity <- function(x, bins) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(rep(0, length(bins) - 1L))
  }

  # np.histogram ignores values outside the supplied breaks.
  # Match that behaviour before calling stats::hist().
  x <- x[x >= bins[1] & x <= bins[length(bins)]]
  if (length(x) == 0L) {
    return(rep(0, length(bins) - 1L))
  }

  stats::hist(
    x,
    breaks = bins,
    plot = FALSE,
    include.lowest = TRUE,
    right = FALSE
  )$density
}

#' Calculate an fbeta positivity threshold from unstim and stim vectors
#'
#' This mirrors get_positivity_threshold() in fbeta.py, with the typo in the
#' Python function name corrected internally.
#'
#' @keywords internal
.simCompareFbetaThreshold <- function(
  xUns,
  xStim,
  beta = 0.8,
  theta = 2,
  width = 10,
  numBins = NULL
) {
  xUns <- as.numeric(xUns)
  xStim <- as.numeric(xStim)
  xUns <- xUns[is.finite(xUns)]
  xStim <- xStim[is.finite(xStim)]

  if (length(xUns) < 2L || length(xStim) < 2L) {
    return(list(
      threshold = NA_real_,
      thresholdMetric = NA_real_,
      thresholdOrigin = "failed_too_few_cells",
      pdfx = numeric(),
      pdfneg = numeric(),
      pdfpos = numeric(),
      fscores = numeric(),
      precision = numeric(),
      recall = numeric()
    ))
  }

  if (is.null(numBins)) {
    numBins <- floor(sqrt(max(length(xUns), length(xStim))))
  }
  numBins <- max(2L, as.integer(numBins))

  rng <- range(xUns, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) <= 0) {
    return(list(
      threshold = NA_real_,
      thresholdMetric = NA_real_,
      thresholdOrigin = "failed_zero_unstim_range",
      pdfx = numeric(),
      pdfneg = numeric(),
      pdfpos = numeric(),
      fscores = numeric(),
      precision = numeric(),
      recall = numeric()
    ))
  }

  bins <- seq(rng[1], rng[2], length.out = numBins + 1L)
  pdfNeg <- .simCompareHistDensity(xUns, bins)
  pdfPos <- .simCompareHistDensity(xStim, bins)

  pdfNeg <- .simCompareMoveMean(pdfNeg, window = width)
  pdfPos <- .simCompareMoveMean(pdfPos, window = width)

  xs <- (bins[-length(bins)] + bins[-1L]) / 2
  score <- .simCompareFbetaCalculateFscore(
    negPdf = pdfNeg,
    posPdf = pdfPos,
    beta = beta,
    theta = theta
  )

  if (length(score$fscores) == 0L || all(!is.finite(score$fscores))) {
    threshold <- NA_real_
    metric <- NA_real_
    origin <- "failed_no_finite_score"
  } else {
    threshold <- xs[which.max(score$fscores)]
    metric <- max(score$fscores, na.rm = TRUE)
    origin <- if (is.finite(metric) && metric > 0) {
      "calculated"
    } else {
      "calculated_zero_fscore"
    }
  }

  list(
    threshold = threshold,
    thresholdMetric = metric,
    thresholdOrigin = origin,
    pdfx = xs,
    pdfneg = pdfNeg,
    pdfpos = pdfPos,
    fscores = score$fscores,
    precision = score$precision,
    recall = score$recall
  )
}

#' @keywords internal
.simCompareFindPeaks <- function(
  x,
  y = NULL,
  numPeaks = NULL,
  adjust = 2,
  bandwidth = NULL,
  ...
) {
  x <- as.vector(x)

  if (length(x) < 2L) {
    return(data.frame(x = NA_real_, y = NA_real_))
  }

  if (is.null(y)) {
    dens <- if (is.null(bandwidth)) {
      stats::density(x, adjust = adjust, ...)
    } else {
      stats::density(x, bw = bandwidth, adjust = adjust, ...)
    }
  } else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of x and y must be equal.")
    }
    dens <- list(x = x, y = y)
  }

  secondDeriv <- diff(sign(diff(dens$y)))
  whichMaxima <- which(secondDeriv == -2) + 1L
  whichMaxima <- whichMaxima[
    findInterval(dens$x[whichMaxima], range(x, na.rm = TRUE)) == 1L
  ]
  whichMaxima <- whichMaxima[order(dens$y[whichMaxima], decreasing = TRUE)]

  if (length(whichMaxima) > 0L) {
    peaksX <- dens$x[whichMaxima]
    if (is.null(numPeaks) || numPeaks > length(peaksX)) {
      numPeaks <- length(peaksX)
    }
    data.frame(
      x = peaksX[seq_len(numPeaks)],
      y = dens$y[whichMaxima][seq_len(numPeaks)]
    )
  } else {
    data.frame(x = NA_real_, y = NA_real_)
  }
}

#' @keywords internal
.simCompareFindValleys <- function(
  x,
  y = NULL,
  numValleys = NULL,
  adjust = 2,
  bandwidth = NULL,
  ...
) {
  x <- as.vector(x)

  if (length(x) < 2L) {
    return(NA_real_)
  }

  if (is.null(y)) {
    dens <- if (is.null(bandwidth)) {
      stats::density(x, adjust = adjust, ...)
    } else {
      stats::density(x, bw = bandwidth, adjust = adjust, ...)
    }
  } else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of x and y must be equal.")
    }
    dens <- list(x = x, y = y)
  }

  secondDeriv <- diff(sign(diff(dens$y)))
  whichMinima <- which(secondDeriv == 2) + 1L
  whichMinima <- whichMinima[
    findInterval(dens$x[whichMinima], range(x, na.rm = TRUE)) == 1L
  ]
  whichMinima <- whichMinima[order(dens$y[whichMinima], decreasing = FALSE)]

  if (length(whichMinima) > 0L) {
    valleys <- dens$x[whichMinima]
    if (is.null(numValleys) || numValleys > length(valleys)) {
      numValleys <- length(valleys)
    }
    valleys[seq_len(numValleys)]
  } else {
    NA_real_
  }
}

#' @keywords internal
.simCompareDerivDensity <- function(
  x,
  deriv = 1,
  bandwidth = NULL,
  adjust = 1,
  numPoints = 10000,
  ...
) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(list(x = numeric(), y = numeric()))
  }

  if (is.null(bandwidth)) {
    bandwidth <- ks::hpi(x, deriv.order = deriv)
  }

  kdeObj <- ks::kdde(
    x = x,
    deriv.order = deriv,
    h = bandwidth * adjust,
    gridsize = numPoints,
    ...
  )

  list(x = kdeObj$eval.points, y = kdeObj$estimate)
}

#' Tailgate-style cytokine cutpoint with optional fixed bandwidth
#'
#' This mirrors cytoUtils:::.cytokineCutpoint(), using the peak and valley logic
#' from openCyto, but keeps the fixed bandwidth explicit.
#'
#' @keywords internal
.simCompareTailgateThreshold <- function(
  x,
  adjust = 1,
  bandwidth = NULL,
  numPeaks = 1,
  refPeak = 1,
  method = c("firstDeriv", "secondDeriv"),
  tol = 1e-2,
  side = "right",
  strict = FALSE,
  autoTol = FALSE,
  numPoints = 10000
) {
  method <- match.arg(method)
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(list(
      threshold = NA_real_,
      thresholdMetric = NA_real_,
      thresholdOrigin = "failed_too_few_unique_cells"
    ))
  }

  peaks <- .simCompareFindPeaks(
    x,
    numPeaks = numPeaks,
    adjust = adjust,
    bandwidth = bandwidth
  )[, "x"]
  peaks <- sort(peaks[is.finite(peaks)])

  if (length(peaks) == 0L) {
    return(list(
      threshold = NA_real_,
      thresholdMetric = NA_real_,
      thresholdOrigin = "failed_no_peak"
    ))
  }

  if (refPeak > length(peaks)) {
    msg <- paste(
      "The reference peak is larger than the number of peaks found.",
      "Setting the reference peak to the last peak."
    )
    if (strict) {
      stop(msg, call. = FALSE)
    }
    refPeak <- length(peaks)
  }

  if (method == "firstDeriv") {
    derivOut <- .simCompareDerivDensity(
      x = x,
      adjust = adjust,
      bandwidth = bandwidth,
      deriv = 1,
      numPoints = numPoints
    )

    if (length(derivOut$x) == 0L) {
      cutpoint <- NA_real_
    } else {
      if (autoTol) {
        tol <- 0.01 * max(abs(derivOut$y), na.rm = TRUE)
      }

      if (side == "right") {
        derivValleys <- .simCompareFindValleys(
          x = derivOut$x,
          y = derivOut$y,
          adjust = adjust,
          bandwidth = bandwidth
        )
        derivValleys <- derivValleys[derivValleys > peaks[refPeak]]
        derivValleys <- sort(derivValleys)[1]

        if (is.na(derivValleys)) {
          cutpoint <- NA_real_
        } else {
          cutpointCandidates <- derivOut$x[
            derivOut$x > derivValleys & abs(derivOut$y) < tol
          ]
          cutpoint <- if (length(cutpointCandidates) > 0L) {
            cutpointCandidates[1]
          } else {
            NA_real_
          }
        }
      } else if (side == "left") {
        derivOut$y <- -derivOut$y
        derivValleys <- .simCompareFindValleys(
          x = derivOut$x,
          y = derivOut$y,
          adjust = adjust,
          bandwidth = bandwidth
        )
        derivValleys <- derivValleys[derivValleys < peaks[refPeak]]
        derivValleys <- sort(derivValleys, decreasing = TRUE)[1]

        if (is.na(derivValleys)) {
          cutpoint <- NA_real_
        } else {
          cutpointCandidates <- derivOut$x[
            derivOut$x < derivValleys & abs(derivOut$y) < tol
          ]
          cutpoint <- if (length(cutpointCandidates) > 0L) {
            cutpointCandidates[length(cutpointCandidates)]
          } else {
            NA_real_
          }
        }
      } else {
        stop("Unrecognised side: ", side)
      }
    }
  } else {
    derivOut <- .simCompareDerivDensity(
      x = x,
      adjust = adjust,
      bandwidth = bandwidth,
      deriv = 2,
      numPoints = numPoints
    )

    if (length(derivOut$x) == 0L) {
      cutpoint <- NA_real_
    } else if (side == "right") {
      derivPeaks <- .simCompareFindPeaks(
        x = derivOut$x,
        y = derivOut$y,
        adjust = adjust,
        bandwidth = bandwidth
      )[, "x"]
      derivPeaks <- derivPeaks[derivPeaks > peaks[refPeak]]
      cutpoint <- sort(derivPeaks)[1]
    } else if (side == "left") {
      derivOut$y <- -derivOut$y
      derivPeaks <- .simCompareFindPeaks(
        x = derivOut$x,
        y = derivOut$y,
        adjust = adjust,
        bandwidth = bandwidth
      )[, "x"]
      derivPeaks <- derivPeaks[derivPeaks < peaks[refPeak]]
      cutpoint <- sort(derivPeaks, decreasing = TRUE)[length(derivPeaks)]
    } else {
      stop("Unrecognised side: ", side)
    }
  }

  list(
    threshold = as.numeric(cutpoint)[1],
    thresholdMetric = NA_real_,
    thresholdOrigin = if (is.finite(cutpoint)) {
      "calculated"
    } else {
      "failed_no_cutpoint"
    }
  )
}

#' @keywords internal
.simCompareHighValueGate <- function(xStim, xUns, margin = 0.05) {
  x <- c(xStim, xUns)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(Inf)
  }
  rng <- range(x, na.rm = TRUE)
  rng[2] + max(1, diff(rng)) * margin
}

#' @keywords internal
.simCompareEstimateFromThreshold <- function(
  xStim,
  xUns,
  threshold,
  fallbackHighValue = TRUE,
  fallbackMargin = 0.05
) {
  xStim <- as.numeric(xStim)
  xUns <- as.numeric(xUns)
  nCellStim <- length(xStim)
  nCellUns <- length(xUns)

  thresholdUsed <- as.numeric(threshold)[1]
  usedFallback <- FALSE
  if (!is.finite(thresholdUsed) && isTRUE(fallbackHighValue)) {
    thresholdUsed <- .simCompareHighValueGate(
      xStim = xStim,
      xUns = xUns,
      margin = fallbackMargin
    )
    usedFallback <- TRUE
  }

  nPosStim <- if (is.finite(thresholdUsed)) {
    sum(xStim > thresholdUsed, na.rm = TRUE)
  } else {
    NA_integer_
  }
  nPosUns <- if (is.finite(thresholdUsed)) {
    sum(xUns > thresholdUsed, na.rm = TRUE)
  } else {
    NA_integer_
  }

  propStim <- nPosStim / nCellStim
  propUns <- nPosUns / nCellUns
  propRespEst <- propStim - propUns

  list(
    threshold = thresholdUsed,
    thresholdFallbackUsed = usedFallback,
    nCellStim = nCellStim,
    nCellUns = nCellUns,
    nPosStim = nPosStim,
    nPosUns = nPosUns,
    propStim = propStim,
    propUns = propUns,
    propRespEst = propRespEst
  )
}

#' @keywords internal
.simCompareTruthTable <- function(
  labelsList,
  nSample,
  nCondition,
  chnl = "F1"
) {
  purrr::map_df(seq_len(nSample), function(sampleCurr) {
    indUns <- (sampleCurr - 1L) * nCondition + 1L
    indStim <- seq.int(indUns + 1L, sampleCurr * nCondition)
    labelVecUns <- labelsList[[indUns]]
    propUnsTruth <- sum(grepl("^gp$", labelVecUns)) / length(labelVecUns)

    purrr::map_df(indStim, function(ind) {
      labelVecStim <- labelsList[[ind]]
      propStimTruth <- sum(grepl("^gp$", labelVecStim)) /
        length(labelVecStim)
      tibble::tibble(
        sample = as.character(sampleCurr),
        ind = as.character(ind),
        chnl = chnl,
        propStimTruth = propStimTruth,
        propUnsTruth = propUnsTruth,
        propRespTruth = propStimTruth - propUnsTruth
      )
    })
  })
}

#' @keywords internal
.simCompareAlternativeRows <- function(
  flowFrameList,
  labelsList,
  nSample,
  nCondition,
  chnl = "F1",
  biasUns = 0,
  fbetaBeta = 0.8,
  fbetaTheta = 2,
  fbetaWidth = 10,
  fbetaNumBins = NULL,
  tailgateX = c("unstim", "combined", "stim"),
  tailgateAdjust = 1,
  tailgateBandwidth = NULL,
  tailgateNumPeaks = 1,
  tailgateRefPeak = 1,
  tailgateMethod = c("firstDeriv", "secondDeriv"),
  tailgateTol = 1e-2,
  tailgateSide = "right",
  tailgateAutoTol = FALSE,
  tailgateNumPoints = 10000,
  fallbackHighValue = TRUE,
  fallbackMargin = 0.05
) {
  tailgateX <- match.arg(tailgateX)
  tailgateMethod <- match.arg(tailgateMethod)

  truthTbl <- .simCompareTruthTable(
    labelsList = labelsList,
    nSample = nSample,
    nCondition = nCondition,
    chnl = chnl
  )

  out <- purrr::map_df(seq_len(nSample), function(sampleCurr) {
    indUns <- (sampleCurr - 1L) * nCondition + 1L
    indStimVec <- seq.int(indUns + 1L, sampleCurr * nCondition)

    xUnsRaw <- as.numeric(flowCore::exprs(flowFrameList[[indUns]])[, chnl])
    xUnsGate <- xUnsRaw + (biasUns %||% 0)

    purrr::map_df(indStimVec, function(indStim) {
      xStim <- as.numeric(flowCore::exprs(flowFrameList[[indStim]])[, chnl])

      fbetaObj <- tryCatch(
        .simCompareFbetaThreshold(
          xUns = xUnsGate,
          xStim = xStim,
          beta = fbetaBeta,
          theta = fbetaTheta,
          width = fbetaWidth,
          numBins = fbetaNumBins
        ),
        error = function(e) {
          list(
            threshold = NA_real_,
            thresholdMetric = NA_real_,
            thresholdOrigin = paste0("error: ", e$message)
          )
        }
      )

      fbetaEst <- .simCompareEstimateFromThreshold(
        xStim = xStim,
        xUns = xUnsGate,
        threshold = fbetaObj$threshold,
        fallbackHighValue = fallbackHighValue,
        fallbackMargin = fallbackMargin
      )

      xTail <- switch(
        tailgateX,
        "unstim" = xUnsGate,
        "combined" = c(xUnsGate, xStim),
        "stim" = xStim
      )

      tailgateObj <- tryCatch(
        .simCompareTailgateThreshold(
          x = xTail,
          adjust = tailgateAdjust,
          bandwidth = tailgateBandwidth,
          numPeaks = tailgateNumPeaks,
          refPeak = tailgateRefPeak,
          method = tailgateMethod,
          tol = tailgateTol,
          side = tailgateSide,
          strict = FALSE,
          autoTol = tailgateAutoTol,
          numPoints = tailgateNumPoints
        ),
        error = function(e) {
          list(
            threshold = NA_real_,
            thresholdMetric = NA_real_,
            thresholdOrigin = paste0("error: ", e$message)
          )
        }
      )

      tailgateEst <- .simCompareEstimateFromThreshold(
        xStim = xStim,
        xUns = xUnsGate,
        threshold = tailgateObj$threshold,
        fallbackHighValue = fallbackHighValue,
        fallbackMargin = fallbackMargin
      )

      tibble::tibble(
        sample = as.character(sampleCurr),
        ind = as.character(indStim),
        chnl = chnl,
        approach = c("fbeta", "tailgate"),
        method = c("fbeta", "tailgate"),
        threshold = c(fbetaEst$threshold, tailgateEst$threshold),
        thresholdOrigin = c(
          fbetaObj$thresholdOrigin,
          tailgateObj$thresholdOrigin
        ),
        gateReturnPoint = c(
          if (isTRUE(fbetaEst$thresholdFallbackUsed)) {
            "fbeta_fallback_high_value"
          } else {
            "fbeta_calculated"
          },
          if (isTRUE(tailgateEst$thresholdFallbackUsed)) {
            "tailgate_fallback_high_value"
          } else {
            "tailgate_calculated"
          }
        ),
        thresholdMetric = c(
          fbetaObj$thresholdMetric %||% NA_real_,
          tailgateObj$thresholdMetric %||% NA_real_
        ),
        thresholdFallbackUsed = c(
          fbetaEst$thresholdFallbackUsed,
          tailgateEst$thresholdFallbackUsed
        ),
        nCellStim = c(fbetaEst$nCellStim, tailgateEst$nCellStim),
        nCellUns = c(fbetaEst$nCellUns, tailgateEst$nCellUns),
        nPosStim = c(fbetaEst$nPosStim, tailgateEst$nPosStim),
        nPosUns = c(fbetaEst$nPosUns, tailgateEst$nPosUns),
        propStim = c(fbetaEst$propStim, tailgateEst$propStim),
        propUns = c(fbetaEst$propUns, tailgateEst$propUns),
        propRespEst = c(fbetaEst$propRespEst, tailgateEst$propRespEst),
        detailLevel = NA_character_,
        locGenerated = NA,
        locGeneratedDirect = NA,
        locSource = NA_character_,
        locReason = NA_character_,
        error = NA_character_
      )
    })
  })

  out |>
    dplyr::left_join(truthTbl, by = c("sample", "ind", "chnl"))
}

#' @keywords internal
.simCompareStimgateFailureRows <- function(truthTbl, errorMessage) {
  truthTbl |>
    dplyr::mutate(
      approach = "stimgate",
      method = "stimgate_error",
      threshold = NA_real_,
      thresholdOrigin = "error",
      gateReturnPoint = "stimgate_error",
      thresholdMetric = NA_real_,
      thresholdFallbackUsed = NA,
      nCellStim = NA_real_,
      nCellUns = NA_real_,
      nPosStim = NA_integer_,
      nPosUns = NA_integer_,
      propStim = NA_real_,
      propUns = NA_real_,
      propRespEst = NA_real_,
      detailLevel = NA_character_,
      locGenerated = NA,
      locGeneratedDirect = NA,
      locSource = NA_character_,
      locReason = NA_character_,
      error = errorMessage
    )
}

#' @keywords internal
.simCompareStimgateRows <- function(
  gs,
  labelsList,
  pathProject,
  nSample,
  nCondition,
  nMarker,
  biasUns,
  bw,
  bwFallback = bw,
  bwMin = "none",
  bwMax = "none",
  bwMtd = "hpi1",
  bwAdj = 1,
  bwNcellMin = 1e2,
  bwNcellMax = 1e5,
  bwCluster = NULL,
  minCell = 1e2,
  maxPosProbX = Inf,
  gateQuant = c(0.25, 0.75),
  locProbCol = "pred",
  locMinPeakProb = 0.25,
  locDipAlpha = 0.2,
  locAntimodeHeightFrac = 1 / 6,
  locAntimodeLowRel = 0.25,
  locAntimodeLowAbs = 0.15,
  locFlatDerivFrac = 1 / 2,
  locFlatHardDerivFrac = 1 / 4,
  locLeftLowRel = 0.25,
  locLeftLowAbs = 0.15,
  locLeftCellFrac = 0.5,
  locLeftLengthFrac = 0.5,
  locMarginalPurityRel = 0.5,
  locMarginalCellBinRatio = 2,
  locMarginalRefQuantile = 0.75,
  locTolRefPeak = "highest",
  gateCombn = "min",
  tolClust = NULL,
  calcCytPosGates = FALSE,
  calcSinglePosGates = FALSE,
  includeLocCondition = TRUE
) {
  truthTbl <- .simCompareTruthTable(
    labelsList = labelsList,
    nSample = nSample,
    nCondition = nCondition,
    chnl = "F1"
  )

  out <- tryCatch(
    {
      batchList <- lapply(seq_len(nSample), function(i) {
        seq((i - 1L) * nCondition + 1L, i * nCondition)
      })

      oldIntermediate <- Sys.getenv("STIMGATE_INTERMEDIATE", unset = NA)
      on.exit(
        {
          if (is.na(oldIntermediate)) {
            Sys.unsetenv("STIMGATE_INTERMEDIATE")
          } else {
            Sys.setenv("STIMGATE_INTERMEDIATE" = oldIntermediate)
          }
        },
        add = TRUE
      )
      Sys.setenv("STIMGATE_INTERMEDIATE" = "TRUE")

      invisible(gateStim(
        .data = gs,
        pathProject = pathProject,
        popGate = "root",
        batchList = batchList,
        marker = paste0("MarkerF", seq_len(nMarker)),
        calcCytPosGates = calcCytPosGates,
        calcSinglePosGates = calcSinglePosGates,
        biasUns = biasUns,
        bw = bw,
        bwFallback = bwFallback,
        bwMin = bwMin,
        bwMax = bwMax,
        bwMtd = bwMtd,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        bwCluster = bwCluster,
        minCell = minCell,
        maxPosProbX = maxPosProbX,
        gateQuant = gateQuant,
        tolClust = tolClust,
        locProbCol = locProbCol,
        locMinPeakProb = locMinPeakProb,
        locDipAlpha = locDipAlpha,
        locAntimodeHeightFrac = locAntimodeHeightFrac,
        locAntimodeLowRel = locAntimodeLowRel,
        locAntimodeLowAbs = locAntimodeLowAbs,
        locFlatDerivFrac = locFlatDerivFrac,
        locFlatHardDerivFrac = locFlatHardDerivFrac,
        locLeftLowRel = locLeftLowRel,
        locLeftLowAbs = locLeftLowAbs,
        locLeftCellFrac = locLeftCellFrac,
        locLeftLengthFrac = locLeftLengthFrac,
        locMarginalPurityRel = locMarginalPurityRel,
        locMarginalCellBinRatio = locMarginalCellBinRatio,
        locMarginalRefQuantile = locMarginalRefQuantile,
        locTolRefPeak = locTolRefPeak,
        gateCombn = gateCombn
      ))

      detailTbl <- .simCompareReadLocDetails(
        pathProject = pathProject,
        nSample = nSample,
        nCondition = nCondition
      )

      if (!includeLocCondition) {
        detailTbl <- detailTbl |>
          dplyr::filter(.data$detailLevel %in% "sample")
      } else {
        detailTbl <- detailTbl |>
          dplyr::filter(.data$detailLevel %in% c("condition", "sample"))
      }

      detailTbl <- .simCompareAddMissingColumns(
        detailTbl,
        list(
          method = NA_character_,
          propRespEst = NA_real_,
          propBsEst = NA_real_,
          propBs = NA_real_,
          threshold = NA_real_,
          thresholdOrigin = NA_character_,
          gateReturnPoint = NA_character_,
          nCellStim = NA_real_,
          nCellUns = NA_real_,
          nPosStim = NA_integer_,
          nPosUns = NA_integer_,
          propStim = NA_real_,
          propUns = NA_real_,
          detailLevel = NA_character_,
          locGenerated = NA,
          locGeneratedDirect = NA,
          locSource = NA_character_,
          locReason = NA_character_
        )
      )

      detailTbl |>
        dplyr::mutate(
          approach = "stimgate",
          method = paste0("stimgate_", .data$method),
          propRespEst = dplyr::coalesce(
            suppressWarnings(as.numeric(.data$propRespEst)),
            suppressWarnings(as.numeric(.data$propBsEst)),
            suppressWarnings(as.numeric(.data$propBs))
          ),
          thresholdMetric = NA_real_,
          thresholdFallbackUsed = grepl(
            "fallback_high_value",
            .data$gateReturnPoint %||% ""
          ),
          error = NA_character_
        ) |>
        dplyr::left_join(truthTbl, by = c("sample", "ind", "chnl")) |>
        dplyr::select(
          sample,
          ind,
          chnl,
          approach,
          method,
          threshold,
          thresholdOrigin,
          gateReturnPoint,
          thresholdMetric,
          thresholdFallbackUsed,
          nCellStim,
          nCellUns,
          nPosStim,
          nPosUns,
          propStim,
          propUns,
          propRespEst,
          propStimTruth,
          propUnsTruth,
          propRespTruth,
          detailLevel,
          locGenerated,
          locGeneratedDirect,
          locSource,
          locReason,
          error,
          dplyr::everything()
        )
    },
    error = function(e) {
      .simCompareStimgateFailureRows(
        truthTbl = truthTbl,
        errorMessage = e$message
      )
    }
  )

  out
}

#' Compare StimGate, fbeta and tailgate estimates on the same simulated data
#'
#' @keywords internal
.simCompareFreqBs <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nIter,
  biasUns,
  bw,
  bwFallback = bw,
  bwMin = "none",
  bwMax = "none",
  bwMtd = "hpi1",
  bwAdj = 1,
  bwNcellMin = 1e2,
  bwNcellMax = 1e5,
  bwCluster = NULL,
  probExact = FALSE,
  nCellStim,
  probResponse,
  meanPos,
  transformation,
  samplePerturbationSd,
  conditionPerturbationSd,
  clusterPerturbationSd,
  backgroundRelativeToResponse,
  ncellUnsRelativeToStim,
  covEvMin = 1,
  covEvMax = 2,
  tolClust = NULL,
  minCell = 1e2,
  maxPosProbX = Inf,
  gateQuant = c(0.25, 0.75),
  locProbCol = "pred",
  locMinPeakProb = 0.25,
  locDipAlpha = 0.2,
  locAntimodeHeightFrac = 1 / 6,
  locAntimodeLowRel = 0.25,
  locAntimodeLowAbs = 0.15,
  locFlatDerivFrac = 1 / 2,
  locFlatHardDerivFrac = 1 / 4,
  locLeftLowRel = 0.25,
  locLeftLowAbs = 0.15,
  locLeftCellFrac = 0.5,
  locLeftLengthFrac = 0.5,
  locMarginalPurityRel = 0.5,
  locMarginalCellBinRatio = 2,
  locMarginalRefQuantile = 0.75,
  locTolRefPeak = "highest",
  gateCombn = "min",
  calcCytPosGates = FALSE,
  calcSinglePosGates = FALSE,
  includeLocCondition = TRUE,
  fbetaBeta = 0.8,
  fbetaTheta = 2,
  fbetaWidth = 10,
  fbetaNumBins = NULL,
  tailgateX = c("unstim", "combined", "stim"),
  tailgateAdjust = 1,
  tailgateBandwidth = bw,
  tailgateNumPeaks = 1,
  tailgateRefPeak = 1,
  tailgateMethod = c("firstDeriv", "secondDeriv"),
  tailgateTol = 1e-2,
  tailgateSide = "right",
  tailgateAutoTol = FALSE,
  tailgateNumPoints = 10000,
  fallbackHighValue = TRUE,
  fallbackMargin = 0.05
) {
  if (!identical(as.integer(nMarker), 1L)) {
    stop("This comparison helper currently expects nMarker = 1.")
  }
  if (!identical(as.integer(nCondition), 2L)) {
    stop("This comparison helper currently expects nCondition = 2.")
  }
  if (!identical(as.integer(nCluster), 2L)) {
    stop("This comparison helper currently expects nCluster = 2.")
  }

  tailgateX <- match.arg(tailgateX)
  tailgateMethod <- match.arg(tailgateMethod)

  purrr::map_df(seq_len(nIter), function(iterNum) {
    nCellUns <- round(nCellStim * ncellUnsRelativeToStim)
    nCellByCondition <- c(nCellUns, nCellStim)
    transformationFunc <- .simCompareGetTrans(transformation)
    meanExprMat <- matrix(
      c(0, meanPos),
      byrow = TRUE,
      ncol = 1
    )
    clusterLabelVec <- c("gn", "gp")
    probResponseUns <- probResponse * backgroundRelativeToResponse
    probVecUns <- c(1 - probResponseUns, probResponseUns)
    probResponseVecByStimCondition <- list(c(-probResponse, probResponse))

    outListExperiment <- simCytExperiment(
      nSample = nSample,
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nCellByCondition = nCellByCondition,
      transformationFunc = transformationFunc,
      mixtureType = "gaussianOnly",
      meanExprMat = meanExprMat,
      clusterLabelVec = clusterLabelVec,
      probVecUns = probVecUns,
      probExact = probExact,
      probResponseVecByStimCondition = probResponseVecByStimCondition,
      conditionPerturbationSd = conditionPerturbationSd,
      clusterPerturbationSd = clusterPerturbationSd,
      samplePerturbationSd = samplePerturbationSd,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    flowFrameList <- outListExperiment[["flowFrameList"]]
    labelsList <- outListExperiment[["labelsList"]]
    fs <- as(flowFrameList, "flowSet")
    gs <- flowWorkspace::GatingSet(fs)

    pathProject <- file.path(
      tempdir(),
      "stimgate",
      "sim-compare-freq-bs",
      paste0(
        "iter-",
        iterNum,
        "-",
        Sys.getpid(),
        "-",
        format(Sys.time(), "%Y%m%d%H%M%S"),
        "-",
        sample.int(1e7, 1)
      )
    )
    on.exit(
      {
        if (dir.exists(pathProject)) {
          unlink(pathProject, recursive = TRUE)
        }
      },
      add = TRUE
    )
    if (dir.exists(pathProject)) {
      unlink(pathProject, recursive = TRUE)
    }
    dir.create(pathProject, recursive = TRUE, showWarnings = FALSE)

    stimgateTbl <- .simCompareStimgateRows(
      gs = gs,
      labelsList = labelsList,
      pathProject = pathProject,
      nSample = nSample,
      nCondition = nCondition,
      nMarker = nMarker,
      biasUns = biasUns,
      bw = bw,
      bwFallback = bwFallback,
      bwMin = bwMin,
      bwMax = bwMax,
      bwMtd = bwMtd,
      bwAdj = bwAdj,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax,
      bwCluster = bwCluster,
      minCell = minCell,
      maxPosProbX = maxPosProbX,
      gateQuant = gateQuant,
      locProbCol = locProbCol,
      locMinPeakProb = locMinPeakProb,
      locDipAlpha = locDipAlpha,
      locAntimodeHeightFrac = locAntimodeHeightFrac,
      locAntimodeLowRel = locAntimodeLowRel,
      locAntimodeLowAbs = locAntimodeLowAbs,
      locFlatDerivFrac = locFlatDerivFrac,
      locFlatHardDerivFrac = locFlatHardDerivFrac,
      locLeftLowRel = locLeftLowRel,
      locLeftLowAbs = locLeftLowAbs,
      locLeftCellFrac = locLeftCellFrac,
      locLeftLengthFrac = locLeftLengthFrac,
      locMarginalPurityRel = locMarginalPurityRel,
      locMarginalCellBinRatio = locMarginalCellBinRatio,
      locMarginalRefQuantile = locMarginalRefQuantile,
      locTolRefPeak = locTolRefPeak,
      gateCombn = gateCombn,
      tolClust = tolClust,
      calcCytPosGates = calcCytPosGates,
      calcSinglePosGates = calcSinglePosGates,
      includeLocCondition = includeLocCondition
    )

    alternativeTbl <- .simCompareAlternativeRows(
      flowFrameList = flowFrameList,
      labelsList = labelsList,
      nSample = nSample,
      nCondition = nCondition,
      chnl = "F1",
      biasUns = biasUns,
      fbetaBeta = fbetaBeta,
      fbetaTheta = fbetaTheta,
      fbetaWidth = fbetaWidth,
      fbetaNumBins = fbetaNumBins,
      tailgateX = tailgateX,
      tailgateAdjust = tailgateAdjust,
      tailgateBandwidth = tailgateBandwidth,
      tailgateNumPeaks = tailgateNumPeaks,
      tailgateRefPeak = tailgateRefPeak,
      tailgateMethod = tailgateMethod,
      tailgateTol = tailgateTol,
      tailgateSide = tailgateSide,
      tailgateAutoTol = tailgateAutoTol,
      tailgateNumPoints = tailgateNumPoints,
      fallbackHighValue = fallbackHighValue,
      fallbackMargin = fallbackMargin
    )

    dplyr::bind_rows(stimgateTbl, alternativeTbl) |>
      dplyr::mutate(
        iter = iterNum,
        nCellStimSim = nCellStim,
        nCellUnsSim = nCellUns,
        biasUns = biasUns,
        bw = bw,
        bwFallback = bwFallback,
        bwMin = bwMin,
        bwMax = bwMax,
        bwMtd = bwMtd,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        bwCluster = bwCluster %||% NA_real_,
        samplePerturbationSd = samplePerturbationSd,
        conditionPerturbationSd = conditionPerturbationSd,
        clusterPerturbationSd = clusterPerturbationSd,
        backgroundRelativeToResponse = backgroundRelativeToResponse,
        ncellUnsRelativeToStim = ncellUnsRelativeToStim,
        fbetaBeta = fbetaBeta,
        fbetaTheta = fbetaTheta,
        fbetaWidth = fbetaWidth,
        fbetaNumBins = fbetaNumBins %||% NA_integer_,
        tailgateX = tailgateX,
        tailgateAdjust = tailgateAdjust,
        tailgateBandwidth = tailgateBandwidth %||% NA_real_,
        tailgateNumPeaks = tailgateNumPeaks,
        tailgateRefPeak = tailgateRefPeak,
        tailgateMethod = tailgateMethod,
        tailgateTol = tailgateTol,
        tailgateSide = tailgateSide,
        tailgateAutoTol = tailgateAutoTol
      ) |>
      dplyr::select(
        iter,
        chnl,
        sample,
        ind,
        approach,
        method,
        dplyr::everything()
      )
  })
}

#' Run the comparison helper over a simulation grid
#'
#' @keywords internal
.simCompareFreqBsGrid <- function(
  sim_grid,
  nSample = 5,
  nIter = 5,
  nMarker = 1,
  nCondition = 2,
  nCluster = 2,
  probExact = TRUE,
  covEvMin = 2,
  covEvMax = 2,
  tolClust = NULL,
  includeLocCondition = TRUE,
  ...
) {
  purrr::map_df(seq_len(nrow(sim_grid)), function(i) {
    row <- sim_grid[i, , drop = FALSE]

    out <- .simCompareFreqBs(
      nSample = nSample,
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nIter = nIter,
      biasUns = row$bias_uns[[1]],
      bw = row$bw[[1]],
      bwFallback = row$bw[[1]],
      bwMin = "none",
      bwMax = "none",
      probExact = probExact,
      nCellStim = row$n_cell[[1]],
      probResponse = row$prob_response[[1]],
      meanPos = row$mean_pos[[1]],
      transformation = row$transformation[[1]],
      samplePerturbationSd = row$sample_perturbation_sd[[1]],
      conditionPerturbationSd = row$condition_perturbation_sd[[1]],
      clusterPerturbationSd = row$cluster_perturbation_sd[[1]],
      backgroundRelativeToResponse = row$background_relative_to_response[[1]],
      ncellUnsRelativeToStim = row$n_cell_uns_relative_to_stim[[1]],
      covEvMin = covEvMin,
      covEvMax = covEvMax,
      tolClust = tolClust,
      includeLocCondition = includeLocCondition,
      ...
    )

    out <- out |>
      dplyr::select(-dplyr::any_of(names(row)))

    dplyr::bind_cols(
      row[rep(1L, nrow(out)), , drop = FALSE],
      out
    )
  })
}

#' Summarise comparison runs by scenario and method
#'
#' @keywords internal
.simCompareSummariseFreqBs <- function(
  .data,
  scenarioCols = NULL,
  keepMethods = c("stimgate_loc_sample", "fbeta", "tailgate")
) {
  if (!"error" %in% names(.data)) {
    .data$error <- NA_character_
  }

  if (is.null(scenarioCols)) {
    scenarioCols <- intersect(
      c(
        "transformation",
        "prob_response",
        "n_cell",
        "mean_pos",
        "bw",
        "bias_uns",
        "sample_perturbation_sd",
        "condition_perturbation_sd",
        "cluster_perturbation_sd",
        "background_relative_to_response",
        "n_cell_uns_relative_to_stim",
        "approach",
        "method"
      ),
      names(.data)
    )
  }

  .data |>
    dplyr::filter(.data$method %in% keepMethods) |>
    dplyr::mutate(
      run_error = .data$error,
      freq_error = .data$propRespEst - .data$propRespTruth,
      abs_error = abs(.data$freq_error),
      sq_error = .data$freq_error^2,
      rel_error = dplyr::if_else(
        .data$propRespTruth != 0,
        .data$freq_error / .data$propRespTruth,
        NA_real_
      )
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenarioCols))) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_est = sum(is.finite(.data$propRespEst)),
      n_run_error = sum(!is.na(.data$run_error) & .data$run_error != ""),
      n_threshold_fallback = sum(
        .data$thresholdFallbackUsed %in% TRUE,
        na.rm = TRUE
      ),
      propRespTruth_mean = mean(.data$propRespTruth, na.rm = TRUE),
      propRespEst_mean = mean(.data$propRespEst, na.rm = TRUE),
      bias = mean(.data$freq_error, na.rm = TRUE),
      mean_abs_error = mean(.data$abs_error, na.rm = TRUE),
      median_abs_error = stats::median(.data$abs_error, na.rm = TRUE),
      rmse = sqrt(mean(.data$sq_error, na.rm = TRUE)),
      mean_rel_error = mean(.data$rel_error, na.rm = TRUE),
      threshold_mean = mean(.data$threshold, na.rm = TRUE),
      threshold_median = stats::median(.data$threshold, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$n_cell,
      dplyr::desc(.data$prob_response),
      .data$transformation,
      .data$mean_pos,
      .data$method
    )
}
