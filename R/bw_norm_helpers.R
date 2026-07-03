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
  normPeakMinRel = 0.75,
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normExtraJitterFrac = 0.25,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normMtd = "moments"
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
      normExcessNcell = normExcessNcell,
      normMtd = normMtd
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

  bwNcellMin <- .bwAsSafeSampleN(
    bwNcellMin,
    default = NULL,
    lower = 2L
  )

  bwNcellMax <- .bwAsSafeSampleN(
    bwNcellMax,
    default = NULL,
    lower = 2L
  )

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
  normPeakMinRel = 0.75,
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normExtraJitterFrac = 0.25,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normMtd = c("moments", "boxcox")
) {
  normMtd <- match.arg(normMtd)

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

  xCore <- x[x <= coreObj$thresholdX]
  xCore <- xCore[is.finite(xCore)]

  if (length(xCore) < 20L || length(unique(xCore)) < 5L) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  boxObj <- NULL
  if (identical(normMtd, "boxcox")) {
    boxObj <- .bwNormChooseBoxCox(
      xCore = xCore,
      lambda = normLambda
    )

    if (is.null(boxObj)) {
      return(.bwCalcOneBase(x, bwMtd))
    }
  }

  # Add synthetic high-side values only when the observed number above the
  # coreset boundary is below the requested high-side fraction.
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

  if (identical(normMtd, "boxcox")) {
    xBw <- c(x, xExtra)
    xBw <- xBw[is.finite(xBw)]

    if (length(xBw) < 20L || length(unique(xBw)) < 5L) {
      return(.bwCalcOneBase(x, bwMtd))
    }

    zBw <- .bwBoxCoxTransform(
      x = xBw,
      lambda = boxObj$lambda,
      winsoriseMin = boxObj$winsoriseMin
    )
  } else {
    sigmaCore <- sd(xCore)
    sigmaExtra <- sd(xExtra)
    sigmaCoreShrink <- 4 / 5 * sigmaCore + 1 / 5 * sigmaExtra
    sigmaExtraShrink <- 1 / 2 * sigmaCore + 1 / 2 * sigmaExtra
    zBw <- c(
      .bwNormSampleNormalComponent(
        mu = mean(xCore),
        sd = sigmaCoreShrink,
        n = length(x) - length(xExtra),
        fallbackSd = .bwRobustSd(x)
      ),
      .bwNormSampleNormalComponent(
        mu = mean(xExtra),
        sd = sigmaExtraShrink,
        n = length(xExtra),
        fallbackSd = .bwRobustSd(xCore)
      )
    )
  }

  zBw <- zBw[is.finite(zBw)]
  zBw <- .bwCalcOneSampleOrdinary(
    x = zBw,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax
  )

  if (length(zBw) < 20L || length(unique(zBw)) < 5L) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  bwZ <- .bwCalcOneBase(
    x = zBw,
    bwMtd = bwMtd
  )

  if (!is.finite(bwZ) || bwZ <= 0) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  # Moment-normalised values are already on the original expression scale.
  if (identical(normMtd, "moments")) {
    return(bwZ)
  }

  zCore <- .bwBoxCoxTransform(
    x = xCore,
    lambda = boxObj$lambda,
    winsoriseMin = boxObj$winsoriseMin
  )

  scaleX <- stats::IQR(xCore, na.rm = TRUE)
  scaleZ <- stats::IQR(zCore, na.rm = TRUE)

  if (!is.finite(scaleX) || scaleX <= 0 || !is.finite(scaleZ) || scaleZ <= 0) {
    return(.bwCalcOneBase(x, bwMtd))
  }

  bwZ * scaleX / scaleZ
}

#' @keywords internal

#' @keywords internal
.bwNormFindBackgroundCore <- function(
  x,
  peakFrac = 0.1,
  peakMinRel = 0.75,
  densityN = 1024L
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 20L || length(unique(x)) < 5L) {
    return(NULL)
  }

  xPilot <- if (length(x) > 1e5L) {
    sample(x, size = 1e5L, replace = FALSE)
  } else {
    x
  }

  bwPilot <- try(
    suppressWarnings(ks::hpi(xPilot, deriv.order = 0)),
    silent = TRUE
  )
  if (
    inherits(bwPilot, "try-error") ||
      !is.finite(bwPilot) ||
      bwPilot <= 0
  ) {
    bwPilot <- stats::IQR(xPilot, na.rm = TRUE) / 20
  }
  if (!is.finite(bwPilot) || bwPilot <= 0) {
    bwPilot <- .bwRobustSd(xPilot) / 5
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

  dx <- suppressWarnings(as.numeric(dens$x))
  dy <- suppressWarnings(as.numeric(dens$y))
  dy <- pmax(dy, 0)

  if (
    length(dx) < 5L ||
      length(dx) != length(dy) ||
      all(!is.finite(dy))
  ) {
    return(NULL)
  }

  peakMainLeftIdx <- .getPeakMainLeftIdx(
    y = dy,
    peakMinRel = peakMinRel,
    troughMaxRel = peakMinRel
  )

  if (
    length(peakMainLeftIdx) != 1L ||
      !is.finite(peakMainLeftIdx) ||
      peakMainLeftIdx < 1L ||
      peakMainLeftIdx > length(dx)
  ) {
    peakMainLeftIdx <- which.max(dy)
  }

  peakMainLeftX <- dx[peakMainLeftIdx]
  peakHeight <- dy[peakMainLeftIdx]

  if (!is.finite(peakHeight) || peakHeight <= 0) {
    return(NULL)
  }

  thresholdTrough <- .bwNormFindBackgroundCoreThresholdTrough(
    dx = dx,
    dy = dy,
    peakMainLeftIdx = peakMainLeftIdx,
    troughMaxRelMain = peakMinRel,
    troughMaxRelNext = peakMinRel,
    troughMaxRelAbs = peakMinRel
  )

  thresholdFlat <- .bwNormFindBackgroundCoreThresholdFlattened(
    dx = dx,
    dy = dy,
    peakMainLeftIdx = peakMainLeftIdx,
    peakMinRel = peakMinRel
  )

  thresholdVec <- c(thresholdTrough, thresholdFlat)
  thresholdVec <- thresholdVec[
    is.finite(thresholdVec) &
      thresholdVec > peakMainLeftX &
      thresholdVec <= max(dx, na.rm = TRUE)
  ]

  thresholdX <- if (length(thresholdVec) > 0L) {
    min(thresholdVec)
  } else {
    max(dx, na.rm = TRUE)
  }

  thresholdIdx <- which(dx >= thresholdX)[1]
  if (!is.finite(thresholdIdx)) {
    thresholdIdx <- length(dx)
    thresholdX <- dx[thresholdIdx]
  }

  list(
    thresholdX = thresholdX,
    thresholdIdx = thresholdIdx,
    peakX = peakMainLeftX,
    peakHeight = peakHeight,
    lowHeight = peakFrac * peakHeight,
    density = tibble::tibble(x = dx, y = dy)
  )
}

.bwNormFindBackgroundCoreThresholdTrough <- function(
  dx,
  dy,
  peakMainLeftIdx,
  troughMaxRelMain = 0.75,
  troughMaxRelNext = 0.75,
  troughMaxRelAbs = 0.75
) {
  dx <- suppressWarnings(as.numeric(dx))
  dy <- suppressWarnings(as.numeric(dy))
  dy <- pmax(dy, 0)

  if (
    length(dx) != length(dy) ||
      length(dy) < 5L ||
      peakMainLeftIdx >= length(dy) - 1L
  ) {
    return(numeric(0L))
  }

  troughIdxAll <- .getLocalMinimaIdx(dy)
  troughIdxAll <- troughIdxAll[troughIdxAll > peakMainLeftIdx]
  if (length(troughIdxAll) == 0L) {
    return(numeric(0L))
  }

  peakIdxAll <- .getLocalMaximaIdx(dy)
  peakIdxAbove <- peakIdxAll[peakIdxAll > peakMainLeftIdx]
  if (length(peakIdxAbove) == 0L) {
    return(numeric(0L))
  }

  peakHeightMain <- dy[peakMainLeftIdx]
  peakHeightAbs <- max(dy, na.rm = TRUE)

  for (troughIdx in troughIdxAll) {
    peakIdxRight <- peakIdxAbove[peakIdxAbove > troughIdx][1L]
    if (!is.finite(peakIdxRight)) {
      next
    }

    troughHeight <- dy[troughIdx]
    peakHeightRight <- dy[peakIdxRight]

    lowEnoughMain <- troughHeight <= troughMaxRelMain * peakHeightMain
    lowEnoughRight <- troughHeight <= troughMaxRelNext * peakHeightRight
    lowEnoughAbs <- troughHeight <= troughMaxRelAbs * peakHeightAbs

    if (
      isTRUE(lowEnoughMain) && isTRUE(lowEnoughRight) && isTRUE(lowEnoughAbs)
    ) {
      return(dx[troughIdx])
    }
  }

  numeric(0L)
}

.bwNormFindBackgroundCoreThresholdFlattened <- function(
  dx,
  dy,
  peakMainLeftIdx,
  peakMinRel = 0.75,
  autoTol = TRUE,
  tol = 1e-8,
  moveBackFrac = 0.1
) {
  dx <- suppressWarnings(as.numeric(dx))
  dy <- suppressWarnings(as.numeric(dy))
  dy <- pmax(dy, 0)

  if (
    length(dx) != length(dy) ||
      length(dx) < 5L ||
      peakMainLeftIdx >= length(dx) - 2L
  ) {
    return(numeric(0L))
  }

  peakHeight <- dy[peakMainLeftIdx]
  if (!is.finite(peakHeight) || peakHeight <= 0) {
    return(numeric(0L))
  }

  # Only look after the density has dropped enough that a shoulder/local wobble
  # near the peak is not mistaken for a tail flattening point.
  rightDropIdx <- which(
    seq_along(dy) > peakMainLeftIdx &
      dy <= peakMinRel * peakHeight
  )[1L]

  if (!is.finite(rightDropIdx) || rightDropIdx >= length(dx) - 1L) {
    return(numeric(0L))
  }

  deriv <- c(NA_real_, diff(dy) / diff(dx))
  derivRight <- deriv[seq.int(rightDropIdx, length(deriv))]
  xRight <- dx[seq.int(rightDropIdx, length(dx))]

  ok <- is.finite(xRight) & is.finite(derivRight)
  xRight <- xRight[ok]
  derivRight <- derivRight[ok]

  if (length(xRight) < 3L) {
    return(numeric(0L))
  }

  negDeriv <- pmax(0, -derivRight)
  if (all(!is.finite(negDeriv)) || max(negDeriv, na.rm = TRUE) <= 0) {
    return(numeric(0L))
  }

  maxDropIdx <- which.max(negDeriv)
  peakDeriv <- negDeriv[maxDropIdx]

  if (!is.finite(peakDeriv) || peakDeriv <= 0) {
    return(numeric(0L))
  }

  thresholdDeriv <- if (isTRUE(autoTol)) {
    peakDeriv / 100
  } else {
    tol
  }

  flatRelIdx <- which(
    seq_along(negDeriv) > maxDropIdx &
      negDeriv <= thresholdDeriv
  )[1L]

  if (!is.finite(flatRelIdx)) {
    return(numeric(0L))
  }

  xFlat <- xRight[flatRelIdx]

  # Move slightly back towards the peak so the coreset includes the main right
  # tail but not the long flat/excess region.
  peakX <- dx[peakMainLeftIdx]
  peakX + (1 - moveBackFrac) * (xFlat - peakX)
}


.bwNormChooseBoxCox <- function(
  xCore,
  lambda = seq(-2, 2, length.out = 81)
) {
  xCore <- suppressWarnings(as.numeric(xCore))
  xCore <- xCore[is.finite(xCore)]

  if (length(xCore) < 20L || length(unique(xCore)) < 5L) {
    return(NULL)
  }

  xMin <- min(xCore, na.rm = TRUE)
  xCore <- xCore[xCore > xMin]

  if (length(xCore) < 20L || length(unique(xCore)) < 5L) {
    return(NULL)
  }

  xCoreQuantVec <- stats::quantile(
    xCore,
    probs = c(0.01, 0.99),
    na.rm = TRUE,
    names = FALSE
  )

  if (any(!is.finite(xCoreQuantVec))) {
    return(NULL)
  }

  xCore <- .winsorise(xCore, probs = c(0.01, 0.99), na.rm = TRUE)
  winsoriseMin <- max(xCoreQuantVec[[1]], .Machine$double.eps)
  xCore <- pmax(xCore, winsoriseMin)

  if (length(xCore) < 20L || length(unique(xCore)) < 5L) {
    return(NULL)
  }

  scoreVec <- purrr::map_dbl(lambda, function(lambdaCurr) {
    z <- .bwBoxCoxTransform(
      x = xCore,
      lambda = lambdaCurr,
      winsoriseMin = winsoriseMin
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
    winsoriseMin = winsoriseMin
  )
}

#' @keywords internal
.bwBoxCoxTransform <- function(
  x,
  lambda,
  winsoriseMin
) {
  x <- pmax(x, winsoriseMin)
  if (abs(lambda) < 1e-8) {
    return(log(x))
  }

  (x^lambda - 1) / lambda
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
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normExtraJitterFrac = 0.25,
  densityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normPeakFrac = 0.1,
  normScamK = 30L
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 20L || length(unique(x)) < 5L) {
    return(numeric(0L))
  }

  nExtraTarget <- ceiling(normExtraFrac * length(x))

  nExtraTarget <- .bwAsSafeSampleN(
    min(nExtraTarget, normExtraMax),
    default = 0L,
    lower = 0L
  )

  # If the observed data already contain enough high-side cells, do not add
  # synthetic high-side values.
  if (is.null(nExtraTarget) || nExtraTarget <= 0L) {
    return(numeric(0L))
  }

  densObj <- .bwNormExcessDensityDecreasing(
    x = x,
    coreObj = coreObj,
    bwMtd = normExcessBwMtd,
    nCell = length(x),
    densityN = densityN,
    peakFrac = normPeakFrac,
    scamK = normScamK
  )

  if (is.null(densObj)) {
    return(numeric(0L))
  }

  candidate <- is.finite(x) & x > coreObj$thresholdX
  if (!any(candidate)) {
    return(numeric(0L))
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
    return(numeric(0L))
  }

  xExtra <- .bwNormPreferentialUpsample(
    x = xCand,
    rate = samplingRate,
    nTarget = nExtraTarget
  )

  if (length(xExtra) == 0L) {
    return(numeric(0L))
  }

  sdDensity <- .sdSum(
    x[x <= coreObj$thresholdX],
    .winsorise(xExtra, probs = c(0, 0.95))
  )

  if (!is.finite(sdDensity) || sdDensity <= 0) {
    sdDensity <- .bwRobustSd(x)
  }
  if (!is.finite(sdDensity) || sdDensity <= 0) {
    sdDensity <- .Machine$double.eps
  }

  xExtra +
    stats::rnorm(
      length(xExtra),
      mean = 0,
      sd = normExtraJitterFrac * sdDensity
    )
}


#' @keywords internal
.bwNormSampleNormalComponent <- function(
  mu,
  sd,
  n = length(x),
  fallbackSd = NULL
) {
  n <- .bwAsSafeSampleN(n, default = 0L, lower = 0L)
  if (is.null(n) || n <= 0L) {
    return(numeric(0L))
  }

  stats::rnorm(n = n, mean = mu, sd = sd)
}

#' @keywords internal
.bwNormCoreTargetN <- function(
  nCore,
  nExtra,
  nTotal,
  bwNcellMin = NULL,
  bwNcellMax = NULL
) {
  nCore <- .bwAsSafeSampleN(nCore, default = 0L, lower = 0L)
  nExtra <- .bwAsSafeSampleN(nExtra, default = 0L, lower = 0L)
  nTotal <- .bwAsSafeSampleN(nTotal, default = nCore + nExtra, lower = 0L)

  if (!is.finite(nCore) || nCore <= 0L) {
    return(0L)
  }

  if (!is.finite(nExtra) || nExtra < 0L) {
    nExtra <- 0L
  }

  # Start from the actual core size.
  # bwNcellMax is a cap, not a target.
  nTarget <- nCore

  bwNcellMaxSafe <- .bwAsSafeSampleN(
    bwNcellMax,
    default = NULL,
    lower = 20L
  )

  if (!is.null(bwNcellMaxSafe)) {
    nTarget <- min(
      nTarget,
      max(20L, bwNcellMaxSafe - nExtra)
    )
  }

  bwNcellMinSafe <- .bwAsSafeSampleN(
    bwNcellMin,
    default = NULL,
    lower = 0L
  )

  if (!is.null(bwNcellMinSafe)) {
    nTarget <- max(
      nTarget,
      bwNcellMinSafe - nExtra
    )
  }

  nTarget <- max(20L, nTarget)

  .bwAsSafeSampleN(
    nTarget,
    default = 20L,
    lower = 20L
  )
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


.winsorise <- function(x, probs = c(0.05, 0.95), na.rm = TRUE) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L || all(!is.finite(x))) {
    return(x)
  }

  qs <- stats::quantile(
    x,
    probs = probs,
    na.rm = na.rm,
    names = FALSE
  )

  if (length(qs) != 2L || any(!is.finite(qs))) {
    return(x)
  }

  pmin(pmax(x, qs[[1L]]), qs[[2L]])
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

  nCellSafe <- .bwAsSafeSampleN(
    nCell,
    default = length(x),
    lower = 20L
  )

  xBw <- if (!is.null(nCellSafe) && length(x) > nCellSafe) {
    sample(x, size = nCellSafe, replace = FALSE)
  } else {
    x
  }

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

  dx <- suppressWarnings(as.numeric(densInit$x))
  dy <- pmax(suppressWarnings(as.numeric(densInit$y)), .Machine$double.eps)

  peakIdx <- which.min(abs(dx - coreObj$peakX))

  if (!is.finite(peakIdx) || peakIdx < 1L || peakIdx > length(dx)) {
    peakIdx <- which.max(dy)
  }

  yDec <- .bwNormFitDecreasingDensity(
    x = x,
    dx = dx,
    dy = dy,
    thresholdX = coreObj$thresholdX,
    peakX = coreObj$peakX,
    peakIdx = peakIdx,
    scamK = scamK
  )

  if (is.null(yDec) || length(yDec) != length(dx)) {
    return(NULL)
  }

  yDec <- pmin(yDec, max(dy))
  yDec <- pmax(yDec, 0)

  list(
    x = dx,
    yInit = dy,
    yDec = yDec,
    bw = bw,
    peakX = coreObj$peakX,
    thresholdX = coreObj$thresholdX,
    peakIdx = peakIdx
  )
}

.bwNormRightCutFromDensity <- function(
  dx,
  dy,
  peakIdx,
  peakFrac = 0.1
) {
  lowHeight <- peakFrac * dy[peakIdx]

  rightIdx <- which(seq_along(dy) > peakIdx & dy <= lowHeight)[1]

  dx[rightIdx]
}

#' @keywords internal

.bwNormFitDecreasingDensity <- function(
  x,
  dx,
  dy,
  thresholdX,
  peakX,
  peakIdx,
  scamK = 30L
) {
  n <- length(dx)

  if (peakIdx >= n - 3L) {
    return(NULL)
  }

  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  yAtX <- stats::approx(
    x = dx,
    y = dy,
    xout = x,
    rule = 2
  )$y

  fitTblInit <- tibble::tibble(
    x = x,
    logDens = log(pmax(yAtX, 1e2 * .Machine$double.eps))
  ) |>
    dplyr::arrange(.data$x) |>
    dplyr::filter(.data$x >= .env$peakX)

  if (nrow(fitTblInit) < 6L) {
    return(NULL)
  }

  # Values beyond the coreset boundary are allowed to influence where the high
  # points are, but not to pull the decreasing background continuation upwards.
  logDensRepRegion <- fitTblInit$logDens[
    fitTblInit$x >= peakX & fitTblInit$x <= thresholdX
  ]

  logDensRepVal <- min(logDensRepRegion, na.rm = TRUE)
  if (!is.finite(logDensRepVal)) {
    logDensRepVal <- min(fitTblInit$logDens, na.rm = TRUE)
  }
  if (!is.finite(logDensRepVal)) {
    return(NULL)
  }

  fitTblInit <- fitTblInit |>
    dplyr::mutate(
      logDens = dplyr::if_else(
        .data$x > .env$thresholdX,
        pmin(.data$logDens, .env$logDensRepVal),
        .data$logDens
      )
    )

  fitTblThin <- .bwNormThinDensityGrid(
    fitTblInit,
    maxPerBin = 20L,
    dx = dx
  )

  if (nrow(fitTblThin) < 6L) {
    return(NULL)
  }

  k <- min(as.integer(scamK), max(4L, nrow(fitTblThin) - 1L))

  fit <- try(
    scam::scam(
      logDens ~ s(x, bs = "mpd", k = k, m = c(2, 1)),
      data = fitTblThin,
      family = gaussian(),
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        maxit = 50
      )
    ),
    silent = TRUE
  )

  yOut <- dy
  predIdx <- dx >= peakX

  if (!inherits(fit, "try-error")) {
    pred <- try(
      stats::predict(
        fit,
        newdata = tibble::tibble(x = dx[predIdx]),
        type = "response"
      ),
      silent = TRUE
    )

    if (!inherits(pred, "try-error") && all(is.finite(pred))) {
      yOut[predIdx] <- exp(pred)
      return(yOut)
    }
  }

  # Isotonic fallback: fit a nonincreasing log-density to the thinned points,
  # then interpolate it back onto the KDE grid.
  iso <- try(
    stats::isoreg(
      seq_along(fitTblThin$x),
      -fitTblThin$logDens
    ),
    silent = TRUE
  )
  if (inherits(iso, "try-error")) {
    return(NULL)
  }

  predFit <- exp(-iso$yf)
  predIso <- stats::approx(
    x = fitTblThin$x,
    y = predFit,
    xout = dx[predIdx],
    rule = 2
  )$y

  if (any(!is.finite(predIso))) {
    return(NULL)
  }

  yOut[predIdx] <- predIso
  yOut
}

#' @keywords internal
.bwNormFitDecreasingDensityIso <- function(
  x,
  y,
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
  maxPerBin = 20L,
  dx = NULL
) {
  if (!is.data.frame(fitTbl) || nrow(fitTbl) == 0L) {
    return(fitTbl)
  }

  maxPerBin <- .bwAsSafeSampleN(maxPerBin, default = 20L, lower = 1L)

  if (is.null(dx) || length(dx) < 2L || any(!is.finite(dx))) {
    breaks <- pretty(fitTbl$x, n = max(2L, ceiling(nrow(fitTbl) / maxPerBin)))
  } else {
    breaks <- sort(unique(as.numeric(dx)))
  }

  if (length(breaks) < 2L) {
    return(fitTbl)
  }

  fitTbl$bin <- cut(
    fitTbl$x,
    breaks = breaks,
    include.lowest = TRUE
  )

  fitTbl |>
    dplyr::filter(!is.na(.data$bin)) |>
    dplyr::group_by(.data$bin) |>
    dplyr::mutate(.bw_rand = stats::runif(dplyr::n())) |>
    dplyr::arrange(.data$bin, .data$.bw_rand) |>
    dplyr::slice_head(n = maxPerBin) |>
    dplyr::ungroup() |>
    dplyr::select(-.data$bin, -.data$.bw_rand)
}

#' @keywords internal

.bwNormPreferentialUpsample <- function(
  x,
  rate,
  nTarget = NULL
) {
  x <- suppressWarnings(as.numeric(x))
  rate <- suppressWarnings(as.numeric(rate))

  ok <- is.finite(x) & is.finite(rate) & rate > 0
  x <- x[ok]
  rate <- rate[ok]

  if (length(x) == 0L) {
    return(numeric(0L))
  }

  nTarget <- .bwAsSafeSampleN(
    nTarget,
    default = 0L,
    lower = 0L
  )

  if (is.null(nTarget) || nTarget <= 0L) {
    return(numeric(0L))
  }

  rate <- pmin(1, pmax(0, rate))
  if (!any(rate > 0)) {
    return(numeric(0L))
  }

  sample(
    x = x,
    size = nTarget,
    replace = TRUE,
    prob = rate
  )
}

#' @keywords internal
.bwAsSafeSampleN <- function(
  x,
  default = NULL,
  lower = 0L,
  upper = .Machine$integer.max
) {
  if (is.null(x) || length(x) == 0L) {
    return(default)
  }

  x <- suppressWarnings(as.numeric(x)[1])

  if (!is.finite(x)) {
    return(default)
  }

  x <- floor(x)
  x <- max(as.numeric(lower), min(x, as.numeric(upper)))

  as.integer(x)
}


.sdSum <- function(x1, x2) {
  .sdOne <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) < 2L) {
      return(0)
    }
    out <- stats::sd(x)
    if (!is.finite(out) || out < 0) {
      return(0)
    }
    out
  }

  sd1 <- .sdOne(x1)
  sd2 <- .sdOne(x2)
  sqrt(sd1^2 + sd2^2)
}
