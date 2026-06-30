# Shared bandwidth helpers for standard and *Norm bandwidth methods.
# Intended method names:
#   nrd0, sj, hpi0, hpi1, hpi2, hpi3
#   nrd0Norm, sjNorm, hpi0Norm, hpi1Norm, hpi2Norm, hpi3Norm

#' @keywords internal
.bwMethodIsNorm <- function(bwMtd) {
  is.character(bwMtd) &&
    length(bwMtd) == 1L &&
    grepl("Norm$", bwMtd)
}

#' @keywords internal
.bwMethodBase <- function(bwMtd) {
  sub("Norm$", "", bwMtd)
}

#' Calculate bandwidth using ordinary or background-normalised methods
#' @keywords internal
.bwCalcOne <- function(
  x,
  bwMtd,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  normPeakFrac = 0.1,
  normPeakMinRel = 0.2,
  normExtraFrac = 0.1,
  normExtraMax = 1000L,
  normExtraJitterFrac = 0.25,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(NA_real_)
  }

  bwMtd <- as.character(bwMtd)[1]
  isNorm <- .bwMethodIsNorm(bwMtd)
  bwMtdBase <- .bwMethodBase(bwMtd)

  if (isTRUE(isNorm)) {
    bwOut <- .bwCalcOneNorm(
      x = x,
      bwMtd = bwMtdBase,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax,
      normPeakFrac = normPeakFrac,
      normPeakMinRel = normPeakMinRel,
      normExtraFrac = normExtraFrac,
      normExtraMax = normExtraMax,
      normExtraJitterFrac = normExtraJitterFrac,
      normLambda = normLambda,
      normDensityN = normDensityN,
      normExcessBwMtd = normExcessBwMtd,
      normExcessNcell = normExcessNcell
    )
  } else {
    xBw <- .bwCalcOneSampleOrdinary(
      x = x,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax
    )

    bwOut <- .bwCalcOneBase(
      x = xBw,
      bwMtd = bwMtdBase
    )
  }

  if (!is.finite(bwOut) || bwOut <= 0) {
    return(NA_real_)
  }

  as.numeric(bwOut)[1] * bwAdj
}

#' @keywords internal
.bwCalcOneSampleOrdinary <- function(
  x,
  bwNcellMin = NULL,
  bwNcellMax = NULL
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(x)
  }

  if (!is.null(bwNcellMin) && length(x) < bwNcellMin) {
    sdX <- .bwRobustSd(x)
    x <- sample(x, replace = TRUE, size = bwNcellMin) +
      stats::rnorm(bwNcellMin, mean = 0, sd = sdX / 10)
  }

  if (!is.null(bwNcellMax) && length(x) > bwNcellMax) {
    x <- sample(x, size = bwNcellMax, replace = FALSE)
  }

  x
}

#' @keywords internal
.bwCalcOneBase <- function(x, bwMtd) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(NA_real_)
  }

  bwCalc <- switch(
    bwMtd,
    "nrd0" = try(stats::bw.nrd0(x), silent = TRUE),
    "sj" = try(stats::bw.SJ(x), silent = TRUE),
    {
      derivOrder <- suppressWarnings(
        as.numeric(gsub("^hpi", "", bwMtd))
      )

      if (!is.finite(derivOrder)) {
        return(NA_real_)
      }

      try(
        suppressWarnings(
          ks::hpi(x, deriv.order = derivOrder)
        ),
        silent = TRUE
      )
    }
  )

  if (
    inherits(bwCalc, "try-error") ||
      !is.finite(bwCalc) ||
      bwCalc <= 0
  ) {
    bwCalc <- try(stats::bw.nrd0(x), silent = TRUE)
  }

  if (
    inherits(bwCalc, "try-error") ||
      !is.finite(bwCalc) ||
      bwCalc <= 0
  ) {
    return(NA_real_)
  }

  as.numeric(bwCalc)[1]
}

#' @keywords internal
.bwCalcOneNorm <- function(
  x,
  bwMtd,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  normPeakFrac = 0.1,
  normPeakMinRel = 0.2,
  normExtraFrac = 0.1,
  normExtraMax = 1000L,
  normExtraJitterFrac = 0.25,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 20L || length(unique(x)) < 5L) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  coreObj <- .bwNormFindBackgroundCore(
    x = x,
    peakFrac = normPeakFrac,
    peakMinRel = normPeakMinRel,
    densityN = normDensityN
  )

  if (is.null(coreObj)) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  xCore <- x[x >= coreObj$xLeft & x <= coreObj$xRight]
  xCore <- xCore[is.finite(xCore)]

  if (length(xCore) < 20L || length(unique(xCore)) < 5L) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  boxObj <- .bwNormChooseBoxCox(
    xCore = xCore,
    lambda = normLambda
  )

  if (is.null(boxObj)) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  xExtra <- .bwNormSampleExcess(
    x = x,
    coreObj = coreObj,
    normExtraFrac = normExtraFrac,
    normExtraMax = normExtraMax,
    normExtraJitterFrac = normExtraJitterFrac,
    densityN = normDensityN,
    normExcessBwMtd = normExcessBwMtd,
    normExcessNcell = normExcessNcell,
    normPeakFrac = normPeakFrac
  )

  nExtra <- length(xExtra)
  nCoreTarget <- .bwNormCoreTargetN(
    nCore = length(xCore),
    nExtra = nExtra,
    nTotal = length(x),
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax
  )

  xCoreBw <- sample(
    xCore,
    size = nCoreTarget,
    replace = length(xCore) < nCoreTarget
  )

  xBw <- c(xCoreBw, xExtra)
  xBw <- xBw[is.finite(xBw)]

  if (length(xBw) < 20L || length(unique(xBw)) < 5L) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  zBw <- .bwBoxCoxTransform(
    x = xBw,
    lambda = boxObj$lambda,
    offset = boxObj$offset,
    eps = boxObj$eps
  )

  bwZ <- .bwCalcOneBase(
    x = zBw,
    bwMtd = bwMtd
  )

  if (!is.finite(bwZ) || bwZ <= 0) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  zCore <- .bwBoxCoxTransform(
    x = xCore,
    lambda = boxObj$lambda,
    offset = boxObj$offset,
    eps = boxObj$eps
  )

  scaleX <- stats::IQR(xCore, na.rm = TRUE)
  scaleZ <- stats::IQR(zCore, na.rm = TRUE)

  if (!is.finite(scaleX) || scaleX <= 0 || !is.finite(scaleZ) || scaleZ <= 0) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  bwZ * scaleX / scaleZ
}

#' @keywords internal
.bwNormFindBackgroundCore <- function(
  x,
  peakFrac = 0.1,
  peakMinRel = 0.2,
  densityN = 512L
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 20L || length(unique(x)) < 5L) {
    return(NULL)
  }

  bwPilot <- try(stats::bw.nrd0(x), silent = TRUE)
  if (
    inherits(bwPilot, "try-error") ||
      !is.finite(bwPilot) ||
      bwPilot <= 0
  ) {
    bwPilot <- stats::IQR(x, na.rm = TRUE) / 20
  }
  if (!is.finite(bwPilot) || bwPilot <= 0) {
    return(NULL)
  }

  dens <- try(
    stats::density(x, bw = bwPilot, n = densityN),
    silent = TRUE
  )

  if (inherits(dens, "try-error")) {
    return(NULL)
  }

  dx <- dens$x
  dy <- dens$y

  if (length(dx) < 5L || all(!is.finite(dy))) {
    return(NULL)
  }

  peakIdxAll <- .bwLocalMaxima(dy)

  if (length(peakIdxAll) == 0L) {
    peakIdx <- which.max(dy)
  } else {
    peakHeightMax <- max(dy, na.rm = TRUE)
    peakIdxMeaningful <- peakIdxAll[
      dy[peakIdxAll] >= peakMinRel * peakHeightMax
    ]

    peakIdx <- if (length(peakIdxMeaningful) > 0L) {
      peakIdxMeaningful[1]
    } else {
      peakIdxAll[1]
    }
  }

  peakHeight <- dy[peakIdx]

  if (!is.finite(peakHeight) || peakHeight <= 0) {
    return(NULL)
  }

  lowHeight <- peakFrac * peakHeight
  minIdxAll <- .bwLocalMinima(dy)

  leftLowIdx <- rev(which(seq_along(dy) < peakIdx & dy <= lowHeight))[1]
  leftIdx <- if (is.finite(leftLowIdx)) leftLowIdx else 1L

  rightAntimodeIdx <- minIdxAll[minIdxAll > peakIdx][1]
  rightLowIdx <- which(seq_along(dy) > peakIdx & dy <= lowHeight)[1]

  rightCandidate <- c(rightAntimodeIdx, rightLowIdx)
  rightCandidate <- rightCandidate[is.finite(rightCandidate)]

  rightIdx <- if (length(rightCandidate) > 0L) {
    min(rightCandidate)
  } else {
    length(dx)
  }

  if (!is.finite(dx[leftIdx]) || !is.finite(dx[rightIdx])) {
    return(NULL)
  }
  if (dx[rightIdx] <= dx[leftIdx]) {
    return(NULL)
  }

  list(
    xLeft = dx[leftIdx],
    xRight = dx[rightIdx],
    xPeak = dx[peakIdx],
    peakHeight = peakHeight,
    lowHeight = lowHeight,
    density = tibble::tibble(x = dx, y = dy)
  )
}

#' @keywords internal
.bwLocalMaxima <- function(y) {
  if (length(y) < 3L) {
    return(integer(0))
  }

  idx <- seq.int(2L, length(y) - 1L)
  idx[y[idx] >= y[idx - 1L] & y[idx] > y[idx + 1L]]
}

#' @keywords internal
.bwLocalMinima <- function(y) {
  if (length(y) < 3L) {
    return(integer(0))
  }

  idx <- seq.int(2L, length(y) - 1L)
  idx[y[idx] <= y[idx - 1L] & y[idx] < y[idx + 1L]]
}

#' @keywords internal
.bwNormChooseBoxCox <- function(
  xCore,
  lambda = seq(-2, 2, length.out = 81)
) {
  xCore <- suppressWarnings(as.numeric(xCore))
  xCore <- xCore[is.finite(xCore)]

  if (length(xCore) < 20L || length(unique(xCore)) < 5L) {
    return(NULL)
  }

  xRange <- diff(range(xCore, na.rm = TRUE))
  eps <- max(.Machine$double.eps, xRange * 1e-8)
  offset <- min(xCore, na.rm = TRUE) - eps

  scoreVec <- purrr::map_dbl(lambda, function(lambdaCurr) {
    z <- .bwBoxCoxTransform(
      x = xCore,
      lambda = lambdaCurr,
      offset = offset,
      eps = eps
    )

    .bwNormalityMomentScore(z)
  })

  if (all(!is.finite(scoreVec))) {
    return(NULL)
  }

  i <- which.min(scoreVec)

  list(
    lambda = lambda[[i]],
    score = scoreVec[[i]],
    offset = offset,
    eps = eps
  )
}

#' @keywords internal
.bwBoxCoxTransform <- function(
  x,
  lambda,
  offset,
  eps
) {
  xShift <- pmax(x - offset, eps)

  if (abs(lambda) < 1e-8) {
    return(log(xShift))
  }

  (xShift^lambda - 1) / lambda
}

#' @keywords internal
.bwNormalityMomentScore <- function(z) {
  z <- suppressWarnings(as.numeric(z))
  z <- z[is.finite(z)]

  if (length(z) < 5L || length(unique(z)) < 3L) {
    return(Inf)
  }

  zMean <- mean(z, na.rm = TRUE)
  zSd <- stats::sd(z, na.rm = TRUE)

  if (!is.finite(zSd) || zSd <= 0) {
    return(Inf)
  }

  zStd <- (z - zMean) / zSd

  skewness <- mean(zStd^3, na.rm = TRUE)
  kurtosis <- mean(zStd^4, na.rm = TRUE)

  if (!is.finite(skewness) || !is.finite(kurtosis)) {
    return(Inf)
  }

  abs(skewness) + abs(kurtosis - 3) / 3
}

#' @keywords internal
.bwNormSampleExcess <- function(
  x,
  coreObj,
  normExtraFrac = 0.1,
  normExtraMax = 1000L,
  normExtraJitterFrac = 0.25,
  densityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normPeakFrac = 0.1,
  normScamK = 30L
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  nExtraTarget <- min(
    ceiling(normExtraFrac * length(x)),
    as.integer(normExtraMax)
  )

  if (!is.finite(nExtraTarget) || nExtraTarget <= 0L) {
    return(numeric(0))
  }

  densObj <- .bwNormExcessDensityDecreasing(
    x = x,
    coreObj = coreObj,
    bwMtd = normExcessBwMtd,
    nCell = normExcessNcell,
    densityN = densityN,
    peakFrac = normPeakFrac,
    scamK = normScamK
  )

  if (is.null(densObj)) {
    return(numeric(0))
  }

  xRightCut <- densObj$xRightCut

  candidate <- is.finite(x) & x > xRightCut
  if (!any(candidate)) {
    return(numeric(0))
  }

  xCand <- x[candidate]

  initAtCand <- stats::approx(
    x = densObj$x,
    y = densObj$yInit,
    xout = xCand,
    rule = 2
  )$y

  decAtCand <- stats::approx(
    x = densObj$x,
    y = densObj$yDec,
    xout = xCand,
    rule = 2
  )$y

  gammaProb <- ifelse(
    is.finite(initAtCand) &
      initAtCand > 0 &
      is.finite(decAtCand) &
      decAtCand >= 0,
    pmin(1, pmax(0, decAtCand / initAtCand)),
    1
  )

  samplingRate <- pmax(0, 1 - gammaProb)

  if (!any(is.finite(samplingRate) & samplingRate > 0)) {
    return(numeric(0))
  }

  xExtra <- .bwNormPreferentialUpsample(
    x = xCand,
    rate = samplingRate,
    nTarget = nExtraTarget
  )

  if (length(xExtra) == 0L) {
    return(numeric(0))
  }

  sdDensity <- .bwDensitySd(
    x = densObj$x,
    y = densObj$yInit
  )

  if (!is.finite(sdDensity) || sdDensity <= 0) {
    sdDensity <- .bwRobustSd(x)
  }

  if (is.finite(sdDensity) && sdDensity > 0) {
    xExtra <- xExtra +
      stats::rnorm(
        length(xExtra),
        mean = 0,
        sd = normExtraJitterFrac * sdDensity
      )
  }

  xExtra
}


#' @keywords internal
.bwNormCoreTargetN <- function(
  nCore,
  nExtra,
  nTotal,
  bwNcellMin = NULL,
  bwNcellMax = NULL
) {
  if (is.null(bwNcellMax)) {
    nTarget <- nCore
  } else {
    nTarget <- max(20L, as.integer(bwNcellMax) - nExtra)
  }

  if (!is.null(bwNcellMin)) {
    nTarget <- max(nTarget, as.integer(bwNcellMin) - nExtra)
  }

  max(20L, nTarget)
}

#' @keywords internal
.bwRobustSd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(.Machine$double.eps)
  }

  iqrX <- diff(stats::quantile(x, c(0.75, 0.25), na.rm = TRUE))
  sdX <- abs(iqrX) / 1.5

  if (!is.finite(sdX) || sdX <= 0) {
    sdX <- stats::sd(x, na.rm = TRUE)
  }

  if (!is.finite(sdX) || sdX <= 0) {
    sdX <- .Machine$double.eps
  }

  sdX
}

#' @keywords internal
.bwWeightedSd <- function(x, w) {
  x <- suppressWarnings(as.numeric(x))
  w <- suppressWarnings(as.numeric(w))
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]
  w <- w[ok]

  if (length(x) < 2L || sum(w) <= 0) {
    return(NA_real_)
  }

  w <- w / sum(w)
  mu <- sum(w * x)
  sqrt(sum(w * (x - mu)^2))
}

#' @keywords internal
.bwNormExcessDensityDecreasing <- function(
  x,
  coreObj,
  bwMtd = "hpi3",
  nCell = 10000L,
  densityN = 512L,
  peakFrac = 0.1,
  scamK = 30L
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 20L || length(unique(x)) < 5L) {
    return(NULL)
  }

  xBw <- .bwCalcOneSampleOrdinary(
    x = x,
    bwNcellMin = nCell,
    bwNcellMax = nCell
  )

  bw <- .bwCalcOneBase(
    x = xBw,
    bwMtd = bwMtd
  )

  if (!is.finite(bw) || bw <= 0) {
    return(NULL)
  }

  densInit <- try(
    stats::density(
      x,
      bw = bw,
      n = densityN,
      from = min(x, na.rm = TRUE),
      to = max(x, na.rm = TRUE)
    ),
    silent = TRUE
  )

  if (inherits(densInit, "try-error")) {
    return(NULL)
  }

  dx <- densInit$x
  dy <- pmax(densInit$y, .Machine$double.eps)

  peakIdx <- which.min(abs(dx - coreObj$xPeak))

  if (!is.finite(peakIdx) || peakIdx < 1L || peakIdx > length(dx)) {
    peakIdx <- which.max(dy)
  }

  xRightCut <- .bwNormRightCutFromDensity(
    dx = dx,
    dy = dy,
    peakIdx = peakIdx,
    peakFrac = peakFrac
  )

  yDec <- .bwNormFitDecreasingDensityScam(
    dx = dx,
    dy = dy,
    peakIdx = peakIdx,
    scamK = scamK
  )

  if (is.null(yDec)) {
    yDec <- .bwNormFitDecreasingDensityIso(
      dx = dx,
      dy = dy,
      peakIdx = peakIdx
    )
  }

  if (is.null(yDec)) {
    return(NULL)
  }

  yDec <- pmin(yDec, dy)
  yDec <- pmax(yDec, 0)

  list(
    x = dx,
    yInit = dy,
    yDec = yDec,
    bw = bw,
    xPeak = dx[peakIdx],
    peakIdx = peakIdx,
    xRightCut = xRightCut
  )
}

#' @keywords internal
.bwNormRightCutFromDensity <- function(
  dx,
  dy,
  peakIdx,
  peakFrac = 0.1
) {
  lowHeight <- peakFrac * dy[peakIdx]

  idx <- seq.int(2L, length(dy) - 1L)
  minIdx <- idx[dy[idx] <= dy[idx - 1L] & dy[idx] < dy[idx + 1L]]

  rightAntimode <- minIdx[minIdx > peakIdx][1]
  rightLow <- which(seq_along(dy) > peakIdx & dy <= lowHeight)[1]

  rightCandidate <- c(rightAntimode, rightLow)
  rightCandidate <- rightCandidate[is.finite(rightCandidate)]

  rightIdx <- if (length(rightCandidate) > 0L) {
    min(rightCandidate)
  } else {
    length(dx)
  }

  dx[rightIdx]
}

#' @keywords internal
.bwNormFitDecreasingDensityScam <- function(
  dx,
  dy,
  peakIdx,
  scamK = 30L
) {
  n <- length(dx)

  if (peakIdx >= n - 3L) {
    return(NULL)
  }

  fitTbl <- tibble::tibble(
    x = dx[seq.int(peakIdx, n)],
    logDens = log(pmax(dy[seq.int(peakIdx, n)], .Machine$double.eps))
  )

  fitTbl <- .bwNormThinDensityGrid(
    fitTbl,
    maxPerBin = 20L
  )

  if (nrow(fitTbl) < 6L) {
    return(NULL)
  }

  k <- min(
    as.integer(scamK),
    max(4L, nrow(fitTbl) - 1L)
  )

  fit <- try(
    scam::scam(
      logDens ~ s(x, bs = "mpd", k = k, m = c(2, 1)),
      data = fitTbl,
      family = gaussian(),
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        maxit = 50
      )
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(NULL)
  }

  predTbl <- tibble::tibble(
    x = dx[seq.int(peakIdx, n)]
  )

  pred <- try(
    stats::predict(fit, newdata = predTbl, type = "response"),
    silent = TRUE
  )

  if (inherits(pred, "try-error") || any(!is.finite(pred))) {
    return(NULL)
  }

  yOut <- dy
  yOut[seq.int(peakIdx, n)] <- exp(pred)

  # Keep the left side as the initial density. The constrained part only matters
  # to the right of the main peak.
  yOut
}

#' @keywords internal
.bwNormFitDecreasingDensityIso <- function(
  dx,
  dy,
  peakIdx
) {
  n <- length(dx)

  if (peakIdx >= n - 3L) {
    return(NULL)
  }

  yRight <- log(pmax(dy[seq.int(peakIdx, n)], .Machine$double.eps))

  iso <- try(
    stats::isoreg(
      seq_along(yRight),
      -yRight
    ),
    silent = TRUE
  )

  if (inherits(iso, "try-error")) {
    return(NULL)
  }

  yRightDec <- exp(-iso$yf)

  yOut <- dy
  yOut[seq.int(peakIdx, n)] <- yRightDec
  yOut
}

#' @keywords internal
.bwNormThinDensityGrid <- function(
  fitTbl,
  maxPerBin = 20L
) {
  if (nrow(fitTbl) <= maxPerBin) {
    return(fitTbl)
  }

  bin <- cut(
    fitTbl$x,
    breaks = pretty(fitTbl$x, n = ceiling(nrow(fitTbl) / maxPerBin)),
    include.lowest = TRUE
  )

  fitTbl$bin <- bin

  fitTbl |>
    dplyr::group_by(.data$bin) |>
    dplyr::slice(seq_len(min(dplyr::n(), maxPerBin))) |>
    dplyr::ungroup() |>
    dplyr::select(-.data$bin)
}

#' @keywords internal
.bwNormPreferentialUpsample <- function(
  x,
  rate,
  nTarget
) {
  ok <- is.finite(x) & is.finite(rate) & rate > 0
  x <- x[ok]
  rate <- rate[ok]

  if (length(x) == 0L || sum(rate) <= 0) {
    return(numeric(0))
  }

  expectedCopies <- rate / sum(rate) * nTarget

  nCopies <- floor(expectedCopies)
  fracCopies <- expectedCopies - nCopies

  nCopies <- nCopies +
    stats::rbinom(
      n = length(nCopies),
      size = 1L,
      prob = pmin(1, pmax(0, fracCopies))
    )

  if (!any(nCopies > 0L)) {
    return(numeric(0))
  }

  rep(x, nCopies)
}

#' @keywords internal
.bwDensitySd <- function(
  x,
  y
) {
  ok <- is.finite(x) & is.finite(y) & y >= 0
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2L || sum(y) <= 0) {
    return(NA_real_)
  }

  area <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

  if (!is.finite(area) || area <= 0) {
    return(NA_real_)
  }

  y <- y / area

  meanX <- sum(diff(x) * (head(x * y, -1) + tail(x * y, -1)) / 2)

  meanX2 <- sum(diff(x) * (head((x^2) * y, -1) + tail((x^2) * y, -1)) / 2)

  varX <- meanX2 - meanX^2

  if (!is.finite(varX) || varX <= 0) {
    return(NA_real_)
  }

  sqrt(varX)
}
