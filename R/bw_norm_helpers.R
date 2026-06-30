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
  normDensityN = 512L
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
      normDensityN = normDensityN
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
  normDensityN = 512L
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
    boxObj = boxObj,
    normExtraFrac = normExtraFrac,
    normExtraMax = normExtraMax,
    normExtraJitterFrac = normExtraJitterFrac,
    densityN = normDensityN
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

  if (
    !is.finite(scaleX) || scaleX <= 0 ||
      !is.finite(scaleZ) || scaleZ <= 0
  ) {
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
  boxObj,
  normExtraFrac = 0.1,
  normExtraMax = 1000L,
  normExtraJitterFrac = 0.25,
  densityN = 512L
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

  zAll <- .bwBoxCoxTransform(
    x = x,
    lambda = boxObj$lambda,
    offset = boxObj$offset,
    eps = boxObj$eps
  )

  xCore <- x[x >= coreObj$xLeft & x <= coreObj$xRight]
  zCore <- .bwBoxCoxTransform(
    x = xCore,
    lambda = boxObj$lambda,
    offset = boxObj$offset,
    eps = boxObj$eps
  )

  if (length(zCore) < 20L || length(unique(zCore)) < 5L) {
    return(numeric(0))
  }

  densFull <- try(stats::density(zAll, n = densityN), silent = TRUE)
  densBg <- try(stats::density(zCore, n = densityN), silent = TRUE)

  if (inherits(densFull, "try-error") || inherits(densBg, "try-error")) {
    return(numeric(0))
  }

  zPeak <- .bwBoxCoxTransform(
    x = coreObj$xPeak,
    lambda = boxObj$lambda,
    offset = boxObj$offset,
    eps = boxObj$eps
  )

  fullAtPeak <- stats::approx(
    x = densFull$x,
    y = densFull$y,
    xout = zPeak,
    rule = 2
  )$y

  bgAtPeak <- stats::approx(
    x = densBg$x,
    y = densBg$y,
    xout = zPeak,
    rule = 2
  )$y

  if (
    !is.finite(fullAtPeak) ||
      fullAtPeak <= 0 ||
      !is.finite(bgAtPeak) ||
      bgAtPeak <= 0
  ) {
    return(numeric(0))
  }

  bgScale <- fullAtPeak / bgAtPeak

  fullAtCell <- stats::approx(
    x = densFull$x,
    y = densFull$y,
    xout = zAll,
    rule = 2
  )$y

  bgAtCell <- stats::approx(
    x = densBg$x,
    y = densBg$y,
    xout = zAll,
    rule = 2
  )$y

  probBg <- ifelse(
    is.finite(fullAtCell) &
      fullAtCell > 0 &
      is.finite(bgAtCell) &
      bgAtCell >= 0,
    pmin(1, bgScale * bgAtCell / fullAtCell),
    1
  )

  probExcess <- pmax(0, 1 - probBg)

  candidate <- x > coreObj$xRight &
    is.finite(x) &
    is.finite(probExcess) &
    probExcess > 0

  if (!any(candidate)) {
    return(numeric(0))
  }

  xCand <- x[candidate]
  wCand <- probExcess[candidate]
  wCand <- wCand / sum(wCand, na.rm = TRUE)

  if (any(!is.finite(wCand)) || sum(wCand, na.rm = TRUE) <= 0) {
    return(numeric(0))
  }

  xExtra <- sample(
    xCand,
    size = nExtraTarget,
    replace = TRUE,
    prob = wCand
  )

  sdExtra <- .bwWeightedSd(xCand, wCand)

  if (!is.finite(sdExtra) || sdExtra <= 0) {
    sdExtra <- .bwRobustSd(xCand)
  }

  if (is.finite(sdExtra) && sdExtra > 0) {
    xExtra <- xExtra +
      stats::rnorm(
        length(xExtra),
        mean = 0,
        sd = normExtraJitterFrac * sdExtra
      )
  }

  pmax(xExtra, boxObj$offset + boxObj$eps)
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
