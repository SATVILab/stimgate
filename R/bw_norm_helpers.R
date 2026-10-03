# Shared bandwidth helpers for standard and *Norm bandwidth methods.
# Intended method names:
#   nrd0, sj, hpi0, hpi1, hpi2, hpi3
#   nrd0Norm, sjNorm, hpi0Norm, hpi1Norm, hpi2Norm, hpi3Norm

#' Calculate bandwidth using ordinary or background-normalised methods
#' @keywords internal
.bwCalcOne <- function(
  x,
  bwMtd,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  normPeakMinRel = 0.75,
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normAdaptiveNcell = 2500L,
  bwAdaptiveCore = NULL,
  bwAdaptiveExtra = NULL,
  bwAdaptiveCrossover = NULL,
  bwAdaptiveTransitionWidth = 0,
  normMtd = "moments",
  adaptive = FALSE
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(structure(NA_real_, adaptive = adaptive))
  }

  bwMtd <- as.character(bwMtd)[1]
  bwMtdBase <- sub("Norm$", "", bwMtd)

  if (grepl("Norm$", bwMtd)) {
    return(.bwCalcOneNorm(
      x = x,
      bwMtd = bwMtdBase,
      bwAdj = bwAdj,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax,
      normPeakMinRel = normPeakMinRel,
      normExtraFrac = normExtraFrac,
      normExtraMax = normExtraMax,
      normLambda = normLambda,
      normDensityN = normDensityN,
      normExcessBwMtd = normExcessBwMtd,
      normExcessNcell = normExcessNcell,
      normAdaptiveNcell = normAdaptiveNcell,
      bwAdaptiveCore = bwAdaptiveCore,
      bwAdaptiveExtra = bwAdaptiveExtra,
      bwAdaptiveCrossover = bwAdaptiveCrossover,
      bwAdaptiveTransitionWidth = bwAdaptiveTransitionWidth,
      normMtd = normMtd,
      adaptive = adaptive
    ))
  }

  xBw <- .bwCalcOneSampleOrdinary(
    x = x,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax
  )

  bwOut <- .bwCalcOneBase(
    x = xBw,
    bwMtd = bwMtdBase
  )

  if (!is.finite(bwOut) || bwOut <= 0) {
    return(structure(NA_real_, adaptive = FALSE))
  }

  structure(as.numeric(bwOut)[1] * bwAdj, adaptive = FALSE)
}

#' @keywords internal
.bwNormManualBw <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  x <- suppressWarnings(as.numeric(x)[1])
  if (!is.finite(x) || x <= 0) {
    return(NULL)
  }

  x
}

#' @keywords internal
.bwNormBwFromCrossover <- function(
  bin,
  bwCore,
  bwExtra,
  crossover,
  transitionWidth = 0
) {
  bin <- suppressWarnings(as.numeric(bin))
  bwCore <- suppressWarnings(as.numeric(bwCore)[1])
  bwExtra <- suppressWarnings(as.numeric(bwExtra)[1])
  crossover <- suppressWarnings(as.numeric(crossover)[1])
  transitionWidth <- suppressWarnings(as.numeric(transitionWidth)[1])

  if (
    length(bin) == 0L ||
      !is.finite(bwCore) ||
      bwCore <= 0 ||
      !is.finite(bwExtra) ||
      bwExtra <= 0 ||
      !is.finite(crossover)
  ) {
    return(rep(NA_real_, length(bin)))
  }

  if (!is.finite(transitionWidth) || transitionWidth <= 0) {
    return(ifelse(bin <= crossover, bwCore, bwExtra))
  }

  wExtra <- stats::plogis((bin - crossover) / transitionWidth)
  (1 - wExtra) * bwCore + wExtra * bwExtra
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

  bwCalc <- switch(bwMtd,
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
.bwNormTooFew <- function(x) {
  length(x) < 20L || length(unique(x)) < 5L
}

#' @keywords internal
.bwCalcOneNorm <- function(
  x,
  bwMtd,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  normPeakMinRel = 0.75,
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normAdaptiveNcell = 2500L,
  bwAdaptiveCore = NULL,
  bwAdaptiveExtra = NULL,
  bwAdaptiveCrossover = NULL,
  bwAdaptiveTransitionWidth = 0,
  normMtd = c("moments", "boxcox"),
  adaptive = FALSE
) {
  normMtd <- match.arg(normMtd)

  if (isTRUE(adaptive) && identical(normMtd, "boxcox")) {
    stop("Cannot use adaptive bandwidth with boxcox normalisation method.")
  }

  .fallback_scalar <- function() {
    bwFallback <- .bwCalcOneBase(x, bwMtd)
    if (!is.finite(bwFallback) || bwFallback <= 0) {
      return(structure(NA_real_, adaptive = FALSE))
    }
    structure(as.numeric(bwFallback)[1] * bwAdj, adaptive = FALSE)
  }

  if (.bwNormTooFew(x)) {
    return(.fallback_scalar())
  }

  corePilotN <- if (
    identical(normMtd, "moments") &&
      !isTRUE(adaptive)
  ) {
    10000L
  } else {
    100000L
  }

  coreObj <- .bwNormFindBackgroundCore(
    x = x,
    peakMinRel = normPeakMinRel,
    densityN = normDensityN,
    pilotN = corePilotN
  )

  if (is.null(coreObj)) {
    return(.fallback_scalar())
  }

  xCore <- x[x <= coreObj$thresholdX]

  if (.bwNormTooFew(xCore)) {
    return(.fallback_scalar())
  }

  boxObj <- NULL
  if (identical(normMtd, "boxcox")) {
    boxObj <- .bwNormChooseBoxCox(
      xCore = xCore,
      lambda = normLambda
    )

    if (is.null(boxObj)) {
      return(.fallback_scalar())
    }
  }

  # Synthetic high-side values represent the "extra" component for the
  # normalised bandwidth calculation.
  xExtra <- .bwNormSampleExcess(
    x = x,
    coreObj = coreObj,
    normExtraFrac = normExtraFrac,
    normExtraMax = normExtraMax,
    densityN = normDensityN,
    normExcessBwMtd = normExcessBwMtd,
    normExcessNcell = normExcessNcell
  )
  xExtra <- xExtra[is.finite(xExtra)]

  if (identical(normMtd, "boxcox")) {
    zBw <- .bwBoxCoxTransform(
      x = c(x, xExtra),
      lambda = boxObj$lambda,
      winsoriseMin = boxObj$winsoriseMin
    )

    zBw <- zBw[is.finite(zBw)]
    zBw <- .bwCalcOneSampleOrdinary(
      x = zBw,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax
    )

    if (.bwNormTooFew(zBw)) {
      return(.fallback_scalar())
    }

    bwZ <- .bwCalcOneBase(
      x = zBw,
      bwMtd = bwMtd
    )

    if (!is.finite(bwZ) || bwZ <= 0) {
      return(.fallback_scalar())
    }

    zCore <- .bwBoxCoxTransform(
      x = xCore,
      lambda = boxObj$lambda,
      winsoriseMin = boxObj$winsoriseMin
    )

    scaleX <- stats::IQR(xCore, na.rm = TRUE)
    scaleZ <- stats::IQR(zCore, na.rm = TRUE)

    if (
      !is.finite(scaleX) || scaleX <= 0 || !is.finite(scaleZ) || scaleZ <= 0
    ) {
      return(.fallback_scalar())
    }

    return(structure(
      as.numeric(bwZ)[1] * scaleX / scaleZ * bwAdj,
      adaptive = FALSE
    ))
  }

  if (!isTRUE(adaptive)) {
    sdCore <- .bwRobustSd(xCore)

    sdExtra <- if (length(xExtra) >= 2L) {
      .bwRobustSd(xExtra)
    } else {
      sdCore
    }

    sdAll <- .bwRobustSd(x)

    sigmaCoreShrink <-
      4 / 5 * sdCore + 1 / 5 * sdExtra

    sigmaExtraShrink <-
      1 / 2 * sdCore + 1 / 2 * sdExtra

    nExtra <- length(xExtra)
    nCore <- max(
      20L,
      length(x) - nExtra
    )

    muCore <- mean(
      xCore,
      na.rm = TRUE
    )

    muExtra <- if (nExtra > 0L) {
      mean(
        xExtra,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }

    zBw <- .bwNormSampleNormalMixture(
      muCore = muCore,
      sdCore = sigmaCoreShrink,
      nCore = nCore,
      fallbackSdCore = sdAll,
      muExtra = muExtra,
      sdExtra = sigmaExtraShrink,
      nExtra = nExtra,
      fallbackSdExtra = sdCore,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax
    )

    zBw <- zBw[is.finite(zBw)]

    if (.bwNormTooFew(zBw)) {
      return(.fallback_scalar())
    }

    bwZ <- .bwCalcOneBase(
      x = zBw,
      bwMtd = bwMtd
    )

    if (
      !is.finite(bwZ) ||
        bwZ <= 0
    ) {
      return(.fallback_scalar())
    }

    return(
      structure(
        as.numeric(bwZ)[1] * bwAdj,
        adaptive = FALSE
      )
    )
  }

  # Adaptive normalised bandwidth: estimate separate component bandwidths on
  # fixed-size normalised core/extra components, then blend them by component
  # density over an expression grid.
  if (.bwNormTooFew(xExtra)) {
    return(.fallback_scalar())
  }

  nAdaptive <- .bwAsSafeSampleN(
    normAdaptiveNcell,
    default = 2500L,
    lower = 20L
  )

  zCore <- .bwNormSampleNormalComponent(
    mu = mean(xCore, na.rm = TRUE),
    sd = .bwRobustSd(xCore),
    n = nAdaptive,
    fallbackSd = .bwRobustSd(x)
  )

  zExtra <- .bwNormSampleNormalComponent(
    mu = mean(xExtra, na.rm = TRUE),
    sd = .bwRobustSd(xExtra),
    n = nAdaptive,
    fallbackSd = .bwRobustSd(xCore)
  )

  zCore <- zCore[is.finite(zCore)]
  zExtra <- zExtra[is.finite(zExtra)]

  if (.bwNormTooFew(zCore) || .bwNormTooFew(zExtra)) {
    return(.fallback_scalar())
  }

  bwZCore <- .bwCalcOneBase(
    x = zCore,
    bwMtd = bwMtd
  )
  bwZExtra <- .bwCalcOneBase(
    x = zExtra,
    bwMtd = bwMtd
  )

  bwManualCore <- .bwNormManualBw(bwAdaptiveCore)
  bwManualExtra <- .bwNormManualBw(bwAdaptiveExtra)

  bwAdjSafe <- suppressWarnings(as.numeric(bwAdj)[1])
  if (!is.finite(bwAdjSafe) || bwAdjSafe <= 0) {
    bwAdjSafe <- 1
  }

  if (!is.null(bwManualCore)) {
    bwZCore <- bwManualCore / bwAdjSafe
  }
  if (!is.null(bwManualExtra)) {
    bwZExtra <- bwManualExtra / bwAdjSafe
  }

  if (
    !is.finite(bwZCore) || bwZCore <= 0 || !is.finite(bwZExtra) || bwZExtra <= 0
  ) {
    return(.fallback_scalar())
  }

  bwZCore <- as.numeric(bwZCore)[1] * bwAdj
  bwZExtra <- as.numeric(bwZExtra)[1] * bwAdj

  rangeVec <- range(c(xCore, xExtra), na.rm = TRUE)
  if (!all(is.finite(rangeVec)) || diff(rangeVec) <= 0) {
    return(.fallback_scalar())
  }

  rangePad <- 0.01 * diff(rangeVec)
  rangeVec <- rangeVec + c(-rangePad, rangePad)

  binVec <- seq(
    from = rangeVec[[1]],
    to = rangeVec[[2]],
    length.out = .bwAsSafeSampleN(normDensityN, default = 512L, lower = 32L)
  )

  densZCore <- try(
    stats::density(
      x = zCore,
      bw = bwZCore,
      n = length(binVec),
      from = min(binVec, na.rm = TRUE),
      to = max(binVec, na.rm = TRUE)
    ),
    silent = TRUE
  )

  densZExtra <- try(
    stats::density(
      x = zExtra,
      bw = bwZExtra,
      n = length(binVec),
      from = min(binVec, na.rm = TRUE),
      to = max(binVec, na.rm = TRUE)
    ),
    silent = TRUE
  )

  if (inherits(densZCore, "try-error") || inherits(densZExtra, "try-error")) {
    return(.fallback_scalar())
  }

  densZCoreY <- pmax(suppressWarnings(as.numeric(densZCore$y)), 0)
  densZExtraY <- pmax(suppressWarnings(as.numeric(densZExtra$y)), 0)

  crossover <- suppressWarnings(as.numeric(bwAdaptiveCrossover)[1])
  if (is.finite(crossover)) {
    bwVec <- .bwNormBwFromCrossover(
      bin = binVec,
      bwCore = bwZCore,
      bwExtra = bwZExtra,
      crossover = crossover,
      transitionWidth = bwAdaptiveTransitionWidth
    )
  } else {
    denom <- densZCoreY + densZExtraY
    bwVec <- ifelse(
      is.finite(denom) & denom > 0,
      (densZCoreY * bwZCore + densZExtraY * bwZExtra) / denom,
      mean(c(bwZCore, bwZExtra))
    )
  }

  bwVec <- pmax(
    suppressWarnings(as.numeric(bwVec)),
    .Machine$double.eps
  )

  structure(
    list(
      bin = binVec,
      bw = bwVec,
      bwCore = bwZCore,
      bwExtra = bwZExtra,
      bwAdaptiveCoreManual = bwManualCore,
      bwAdaptiveExtraManual = bwManualExtra,
      bwAdaptiveCrossover = crossover,
      bwAdaptiveTransitionWidth = suppressWarnings(as.numeric(
        bwAdaptiveTransitionWidth
      )[1]),
      coreObj = coreObj,
      nAdaptive = nAdaptive
    ),
    adaptive = TRUE
  )
}


#' @keywords internal
.bwNormFindBackgroundCore <- function(
  x,
  peakMinRel = 0.75,
  densityN = 1024L,
  pilotN = 100000L
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (.bwNormTooFew(x)) {
    return(NULL)
  }

  pilotN <- .bwAsSafeSampleN(
    pilotN,
    default = 100000L,
    lower = 20L
  )

  xPilot <- if (
    !is.null(pilotN) &&
      length(x) > pilotN
  ) {
    sample(
      x,
      size = pilotN,
      replace = FALSE
    )
  } else {
    x
  }

  bwPilot <- try(
    suppressWarnings(
      ks::hpi(
        xPilot,
        deriv.order = 0
      )
    ),
    silent = TRUE
  )

  if (
    inherits(bwPilot, "try-error") ||
      !is.finite(bwPilot) ||
      bwPilot <= 0
  ) {
    bwPilot <- stats::IQR(
      xPilot,
      na.rm = TRUE
    ) /
      20
  }

  if (
    !is.finite(bwPilot) ||
      bwPilot <= 0
  ) {
    bwPilot <- .bwRobustSd(
      xPilot
    ) /
      5
  }

  # Use the complete sample for the actual density used to identify
  # the background core. Only bandwidth selection is performed on
  # the smaller pilot sample.
  dens <- try(
    stats::density(
      x,
      bw = bwPilot,
      n = densityN
    ),
    silent = TRUE
  )

  if (inherits(dens, "try-error")) {
    return(NULL)
  }

  dx <- suppressWarnings(
    as.numeric(dens$x)
  )

  dy <- suppressWarnings(
    as.numeric(dens$y)
  )

  dy <- pmax(
    dy,
    0
  )

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
    peakMainLeftIdx <- which.max(
      dy
    )
  }

  peakMainLeftX <- dx[
    peakMainLeftIdx
  ]

  peakHeight <- dy[
    peakMainLeftIdx
  ]

  if (
    !is.finite(peakHeight) ||
      peakHeight <= 0
  ) {
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

  dxMax <- max(
    dx,
    na.rm = TRUE
  )

  thresholdVec <- c(
    thresholdTrough,
    thresholdFlat
  )

  thresholdVec <- thresholdVec[
    is.finite(thresholdVec) &
      thresholdVec > peakMainLeftX &
      thresholdVec <= dxMax
  ]

  thresholdX <- if (length(thresholdVec) > 0L) {
    min(
      thresholdVec
    )
  } else {
    dxMax
  }

  thresholdIdx <- which(
    dx >= thresholdX
  )[1L]

  if (!is.finite(thresholdIdx)) {
    thresholdIdx <- length(dx)
    thresholdX <- dx[thresholdIdx]
  }

  list(
    thresholdX = thresholdX,
    thresholdIdx = thresholdIdx,
    peakX = peakMainLeftX,
    peakHeight = peakHeight,
    density = tibble::tibble(
      x = dx,
      y = dy
    )
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

  if (.bwNormTooFew(xCore)) {
    return(NULL)
  }

  xMin <- min(xCore, na.rm = TRUE)
  xCore <- xCore[xCore > xMin]

  if (.bwNormTooFew(xCore)) {
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

  xCore <- pmin(pmax(xCore, xCoreQuantVec[[1]]), xCoreQuantVec[[2]])
  winsoriseMin <- max(xCoreQuantVec[[1]], .Machine$double.eps)
  xCore <- pmax(xCore, winsoriseMin)

  if (.bwNormTooFew(xCore)) {
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
  densityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normScamK = 30L
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (.bwNormTooFew(x)) {
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
    nCell = normExcessNcell,
    densityN = densityN,
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

  excessCand <- is.finite(xCand) &
    is.finite(samplingRate) &
    samplingRate > 0

  if (!any(excessCand)) {
    return(numeric(0L))
  }

  xExcessCand <- xCand[excessCand]

  extraLower <- coreObj$thresholdX

  extraUpper <- stats::quantile(
    xExcessCand,
    probs = 0.95,
    na.rm = TRUE,
    names = FALSE
  )

  if (
    !is.finite(extraLower) ||
      !is.finite(extraUpper) ||
      extraUpper <= extraLower
  ) {
    return(numeric(0L))
  }

  xExtra <- stats::runif(
    n = nExtraTarget,
    min = extraLower,
    max = extraUpper
  )

  xExtra <- xExtra[is.finite(xExtra)]

  if (length(xExtra) == 0L) {
    return(numeric(0L))
  }

  # .bwRobustSd() always returns a finite positive value.
  sdCore <- .bwRobustSd(x[x <= coreObj$thresholdX])

  sdExtra <- if (length(unique(xExtra)) >= 2L) {
    stats::sd(xExtra, na.rm = TRUE)
  } else {
    0
  }

  if (!is.finite(sdExtra) || sdExtra < 0) {
    sdExtra <- 0
  }

  sdJitter <- sqrt(
    pmax(0, sdCore^2 - sdExtra^2)
  )

  if (is.finite(sdJitter) && sdJitter > 0) {
    xExtra <- xExtra +
      stats::rnorm(
        length(xExtra),
        mean = 0,
        sd = sdJitter
      )
  }

  xExtra
}


#' @keywords internal
.bwNormSampleNormalComponent <- function(
  mu,
  sd,
  n = NULL,
  fallbackSd = NULL
) {
  n <- .bwAsSafeSampleN(n, default = 0L, lower = 0L)
  if (is.null(n) || n <= 0L) {
    return(numeric(0L))
  }

  mu <- suppressWarnings(as.numeric(mu)[1])
  if (!is.finite(mu)) {
    return(numeric(0L))
  }

  sd <- suppressWarnings(as.numeric(sd)[1])
  fallbackSd <- suppressWarnings(as.numeric(fallbackSd)[1])

  if (!is.finite(sd) || sd <= 0) {
    sd <- fallbackSd
  }
  if (!is.finite(sd) || sd <= 0) {
    sd <- .Machine$double.eps
  }

  stats::rnorm(n = n, mean = mu, sd = sd)
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
.bwNormExcessDensityDecreasing <- function(
  x,
  coreObj,
  bwMtd = "hpi3",
  nCell = 10000L,
  densityN = 512L,
  scamK = 30L
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (.bwNormTooFew(x)) {
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

  # The fit only uses observations on or to the right of the main peak.
  xFit <- x[x >= peakX]

  if (length(xFit) < 6L) {
    return(NULL)
  }

  # Work out the representative log-density using the complete set of
  # observations in the core region. This is deliberately done before
  # thinning so that thinning cannot change this quantity.
  xRep <- xFit[
    xFit <= thresholdX
  ]

  .minLogDens <- function(xout) {
    densOut <- stats::approx(x = dx, y = dy, xout = xout, rule = 2)$y
    min(log(pmax(densOut, 1e2 * .Machine$double.eps)), na.rm = TRUE)
  }

  logDensRepVal <- if (length(xRep) > 0L) .minLogDens(xRep) else NA_real_

  # Preserve the previous fallback. This should rarely be needed because
  # interpolation from a valid KDE should be finite.
  if (!is.finite(logDensRepVal)) {
    logDensRepVal <- .minLogDens(xFit)
  }

  if (!is.finite(logDensRepVal)) {
    return(NULL)
  }

  # Thin the raw observations BEFORE interpolating the KDE onto them and
  # before constructing a data frame. At most maxPerBin observations from
  # each KDE-grid interval are required by the downstream SCAM fit.
  xThin <- .bwNormThinXByDensityGrid(
    x = xFit,
    maxPerBin = 20L,
    dx = dx
  )

  if (length(xThin) < 6L) {
    return(NULL)
  }

  densThin <- stats::approx(
    x = dx,
    y = dy,
    xout = xThin,
    rule = 2
  )$y

  logDensThin <- log(
    pmax(
      densThin,
      1e2 * .Machine$double.eps
    )
  )

  # Values beyond the coreset boundary may determine where high-x points
  # occur, but must not pull the decreasing background continuation upwards.
  aboveThreshold <- xThin > thresholdX

  logDensThin[aboveThreshold] <- pmin(
    logDensThin[aboveThreshold],
    logDensRepVal
  )

  fitTblThin <- tibble::tibble(
    x = xThin,
    logDens = logDensThin
  )

  k <- min(
    as.integer(scamK),
    max(
      4L,
      nrow(fitTblThin) - 1L
    )
  )

  fit <- try(
    scam::scam(
      logDens ~ s(
        x,
        bs = "mpd",
        k = k,
        m = c(2, 1)
      ),
      data = fitTblThin,
      family = stats::gaussian(),
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
        newdata = tibble::tibble(
          x = dx[predIdx]
        ),
        type = "response"
      ),
      silent = TRUE
    )

    if (
      !inherits(pred, "try-error") &&
        all(is.finite(pred))
    ) {
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

  predFit <- exp(
    -iso$yf
  )

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
.bwNormThinXByDensityGrid <- function(
  x,
  maxPerBin = 20L,
  dx
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(x)
  }

  maxPerBin <- .bwAsSafeSampleN(
    maxPerBin,
    default = 20L,
    lower = 1L
  )

  # The only caller passes a finite KDE grid with at least five points.
  breaks <- sort(unique(as.numeric(dx)))

  if (length(breaks) < 2L) {
    return(sort(x))
  }

  # The old implementation arranged by x before binning.
  x <- sort(x)

  bin <- cut(
    x,
    breaks = breaks,
    include.lowest = TRUE
  )

  keep <- !is.na(bin)
  x <- x[keep]
  bin <- bin[keep]

  if (length(x) == 0L) {
    return(x)
  }

  indByBin <- split(
    seq_along(x),
    bin,
    drop = TRUE
  )

  keepInd <- unlist(
    lapply(
      indByBin,
      function(ind) {
        # Generate random priorities for every point, as in the previous
        # grouped runif/arrange/slice implementation.
        rand <- stats::runif(
          length(ind)
        )

        ind[
          order(rand)[
            seq_len(
              min(
                length(ind),
                maxPerBin
              )
            )
          ]
        ]
      }
    ),
    use.names = FALSE
  )

  x[keepInd]
}
#' @keywords internal

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


#' @keywords internal
.bwNormSampleNormalMixture <- function(
  muCore,
  sdCore,
  nCore,
  fallbackSdCore,
  muExtra,
  sdExtra,
  nExtra,
  fallbackSdExtra,
  bwNcellMin = NULL,
  bwNcellMax = NULL
) {
  nCore <- .bwAsSafeSampleN(
    nCore,
    default = 0L,
    lower = 0L
  )
  nExtra <- .bwAsSafeSampleN(
    nExtra,
    default = 0L,
    lower = 0L
  )

  nTotal <- nCore + nExtra

  bwNcellMinSafe <- .bwAsSafeSampleN(
    bwNcellMin,
    default = NULL,
    lower = 2L
  )

  bwNcellMaxSafe <- .bwAsSafeSampleN(
    bwNcellMax,
    default = NULL,
    lower = 2L
  )

  # If the original synthetic mixture would be downsampled, determine how
  # many sampled observations come from each component before generating
  # them. Sampling without replacement from the complete mixture gives a
  # hypergeometric component count.
  canGenerateCapped <-
    !is.null(bwNcellMaxSafe) &&
      nTotal > bwNcellMaxSafe &&
      (is.null(bwNcellMinSafe) ||
        nTotal >= bwNcellMinSafe)

  if (canGenerateCapped) {
    nExtraOut <- if (nExtra > 0L) {
      as.integer(stats::rhyper(
        nn = 1L,
        m = nExtra,
        n = nCore,
        k = bwNcellMaxSafe
      ))
    } else {
      0L
    }
    nCoreOut <- bwNcellMaxSafe - nExtraOut
  } else {
    nCoreOut <- nCore
    nExtraOut <- nExtra
  }

  zBw <- c(
    .bwNormSampleNormalComponent(
      mu = muCore,
      sd = sdCore,
      n = nCoreOut,
      fallbackSd = fallbackSdCore
    ),
    .bwNormSampleNormalComponent(
      mu = muExtra,
      sd = sdExtra,
      n = nExtraOut,
      fallbackSd = fallbackSdExtra
    )
  )

  if (canGenerateCapped) {
    return(zBw)
  }

  # Keep the existing route when no downsampling is needed, or when
  # bwNcellMin would first cause upsampling with jitter.
  zBw <- zBw[is.finite(zBw)]

  .bwCalcOneSampleOrdinary(
    x = zBw,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax
  )
}
