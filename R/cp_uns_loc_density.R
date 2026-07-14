# Local-FDR density and raw probability estimation
#
# Fits stimulated and unstimulated densities, supports fixed and adaptive
# bandwidths, constructs raw response probabilities, and prepares the data used
# by the monotone smoother.

.getCpUnsLocGetProb <- function(
  exTblStimNoMin,
  exTblStimThreshold,
  exTblUnsThreshold,
  exTblUnsBias,
  bias,
  exTblUnsOrig,
  stage,
  pathProject,
  chnlSettings
) {
  ind <- .getInd(exTblStimNoMin)
  chnl <- .getCpUnsLocGetChnl(exTblStimNoMin)
  stageChnl <- file.path(stage, chnl)

  # get raw densities
  densTblRaw <- .getCpUnsLocGetDensRaw(
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold,
    stage = stage,
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )
  .intSave(ind, stageChnl, pathProject, densTblRaw)

  # get probabilities
  probTblList <- .getCpUnsLocGetProbTbl(
    densTblRaw = densTblRaw,
    stage = stage,
    cpMin = chnlSettings$cpMin + bias,
    exVecStimThreshold = .getCut(exTblStimThreshold),
    exVecUnsThreshold = .getCut(exTblUnsThreshold)
  )
  .intSave(ind, stageChnl, pathProject, probTblList)

  # get .data to smooth over
  dataMod <- .getCpUnsLocGetDataMod(
    exTblStimThreshold = exTblStimThreshold,
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsThreshold = exTblUnsThreshold,
    exTblUnsBias = exTblUnsBias,
    probTblList = probTblList,
    cpMin = chnlSettings$cpMin + bias,
    stage = stage
  )
  .intSave(ind, stageChnl, pathProject, dataMod)

  # smooth
  .getCpUnsLocGetProbSmooth(
    dataMod = dataMod,
    stage = stage,
    pathProject = pathProject,
    chnl = chnl,
    chnlSettings = chnlSettings
  )
}

# get densTblRaw
# -----------------------
#' @keywords internal
.getCpUnsLocGetDensRaw <- function(
  exTblStimThreshold,
  exTblUnsThreshold,
  stage,
  pathProject,
  chnlSettings
) {
  .debug("Calculating densities") # nolint

  densList <- .getCpUnsLocGetDensRawDensities(
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold,
    stage = stage,
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )

  # Keep the bandwidth object used here so the later antimode density can use
  # exactly half of the same fixed or adaptive bandwidth.
  densTblRaw <- .getCpUnsLocGetDensRawTabulate(
    stimX = densList$stim$x,
    stimY = densList$stim$y,
    unsX = densList$uns$x,
    unsY = densList$uns$y
  )
  attr(densTblRaw, "locDensityBw") <- densList$bw
  densTblRaw
}

#' @keywords internal
.getCpUnsLocGetDensRawDensities <- function(
  exTblStimThreshold,
  exTblUnsThreshold,
  stage,
  pathProject,
  chnlSettings
) {
  useAdaptive <- .getCpUnsLocUseAdaptiveBw(chnlSettings)

  if (isTRUE(useAdaptive)) {
    densAdaptive <- .getCpUnsLocGetDensRawDensitiesAdaptive(
      exTblStimThreshold = exTblStimThreshold,
      exTblUnsThreshold = exTblUnsThreshold,
      chnlSettings = chnlSettings
    )

    if (!is.null(densAdaptive)) {
      chnl <- .getCpUnsLocGetChnl(exTblStimThreshold)
      stageChnl <- file.path(stage, chnl)
      .intSaveNm(
        "bwCpUnsLocAdaptive",
        densAdaptive$bw,
        .getInd(exTblStimThreshold),
        stageChnl,
        pathProject
      )
      return(densAdaptive)
    }
  }

  bw <- .getCpUnsLocGetDensRawDensitiesBw(
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold,
    bw = chnlSettings$bw,
    bwMin = chnlSettings$bwMin,
    bwMax = chnlSettings$bwMax,
    bwFallback = chnlSettings$bwFallback,
    bwMtd = chnlSettings$bwMtd,
    bwAdj = chnlSettings$bwAdj,
    bwNcellMin = chnlSettings$bwNcellMin,
    bwNcellMax = chnlSettings$bwNcellMax,
    chnlSettings = chnlSettings
  )
  chnl <- .getCpUnsLocGetChnl(exTblStimThreshold)
  stageChnl <- file.path(stage, chnl)
  .intSaveNm(
    "bwCpUnsLoc",
    bw,
    .getInd(exTblStimThreshold),
    stageChnl,
    pathProject
  )
  densStim <- .getCpUnsLocGetDensRawDensitiesStim(
    exTblStimThreshold = exTblStimThreshold,
    bw = bw
  )
  densUns <- .getCpUnsLocGetDensRawDensitiesUns(
    exTblUnsThreshold = exTblUnsThreshold,
    densStim = densStim,
    bw = bw
  )
  list(stim = densStim, uns = densUns, bw = bw)
}


#' @keywords internal
.getCpUnsLocUseAdaptiveBw <- function(chnlSettings) {
  (isTRUE(chnlSettings$bwAdaptive %||% FALSE) ||
    .getCpUnsLocHasManualAdaptiveBw(chnlSettings)) &&
    is.null(chnlSettings$bw)
}

#' @keywords internal
.getCpUnsLocHasManualAdaptiveBw <- function(chnlSettings) {
  is.finite(suppressWarnings(as.numeric(chnlSettings$bwAdaptiveCore)[1])) ||
    is.finite(suppressWarnings(as.numeric(chnlSettings$bwAdaptiveExtra)[1])) ||
    is.finite(suppressWarnings(as.numeric(chnlSettings$bwAdaptiveCrossover)[1]))
}

#' @keywords internal
.getCpUnsLocManualAdaptiveBwOnGrid <- function(chnlSettings, grid) {
  bwCore <- suppressWarnings(as.numeric(chnlSettings$bwAdaptiveCore)[1])
  bwExtra <- suppressWarnings(as.numeric(chnlSettings$bwAdaptiveExtra)[1])
  crossover <- suppressWarnings(as.numeric(chnlSettings$bwAdaptiveCrossover)[1])
  transitionWidth <- suppressWarnings(
    as.numeric(chnlSettings$bwAdaptiveTransitionWidth %||% 0)[1]
  )

  if (
    !is.finite(bwCore) ||
      bwCore <= 0 ||
      !is.finite(bwExtra) ||
      bwExtra <= 0 ||
      !is.finite(crossover)
  ) {
    return(NULL)
  }

  bw <- .bwNormBwFromCrossover(
    bin = grid,
    bwCore = bwCore,
    bwExtra = bwExtra,
    crossover = crossover,
    transitionWidth = transitionWidth
  )

  bw <- .getCpUnsLocRepairBwGrid(bw)

  list(
    bin = grid,
    bw = bw,
    bwCore = bwCore,
    bwExtra = bwExtra,
    bwAdaptiveCoreManual = bwCore,
    bwAdaptiveExtraManual = bwExtra,
    bwAdaptiveCrossover = crossover,
    bwAdaptiveTransitionWidth = transitionWidth,
    adaptive = TRUE,
    manual = TRUE
  )
}

#' @keywords internal
.getCpUnsLocGetDensRawDensitiesAdaptive <- function(
  exTblStimThreshold,
  exTblUnsThreshold,
  chnlSettings
) {
  xStim <- .getCut(exTblStimThreshold)
  xUns <- .getCut(exTblUnsThreshold)
  xStim <- suppressWarnings(as.numeric(xStim))
  xUns <- suppressWarnings(as.numeric(xUns))
  xStim <- xStim[is.finite(xStim)]
  xUns <- xUns[is.finite(xUns)]

  if (
    length(xStim) < 20L ||
      length(unique(xStim)) < 5L ||
      length(xUns) < 20L ||
      length(unique(xUns)) < 5L
  ) {
    return(NULL)
  }

  densityN <- .bwAsSafeSampleN(
    chnlSettings$bwAdaptiveDensityN %||%
      chnlSettings$normDensityN %||%
      512L,
    default = 512L,
    lower = 64L
  )

  grid <- .getCpUnsLocAdaptiveGrid(
    x = c(xStim, xUns),
    n = densityN,
    padFrac = chnlSettings$bwAdaptivePadFrac %||% 0.15
  )

  if (length(grid) < 10L) {
    return(NULL)
  }

  manualBwObj <- .getCpUnsLocManualAdaptiveBwOnGrid(
    chnlSettings = chnlSettings,
    grid = grid
  )

  if (!is.null(manualBwObj)) {
    bwStimObj <- manualBwObj
    bwUnsObj <- manualBwObj
    bwStimGrid <- manualBwObj$bw
    bwUnsGrid <- manualBwObj$bw
  } else {
    bwStimObj <- .getCpUnsLocGetDensRawDensitiesBwAdaptiveOne(
      x = xStim,
      chnlSettings = chnlSettings
    )
    bwUnsObj <- .getCpUnsLocGetDensRawDensitiesBwAdaptiveOne(
      x = xUns,
      chnlSettings = chnlSettings
    )

    bwStimGrid <- .getCpUnsLocAdaptiveBwOnGrid(
      bwObj = bwStimObj,
      grid = grid,
      fallback = chnlSettings$bwFallback
    )
    bwUnsGrid <- .getCpUnsLocAdaptiveBwOnGrid(
      bwObj = bwUnsObj,
      grid = grid,
      fallback = chnlSettings$bwFallback
    )
  }

  if (
    length(bwStimGrid) != length(grid) ||
      length(bwUnsGrid) != length(grid) ||
      all(!is.finite(bwStimGrid)) ||
      all(!is.finite(bwUnsGrid))
  ) {
    return(NULL)
  }

  bwStimGrid <- .getCpUnsLocRepairBwGrid(bwStimGrid)
  bwUnsGrid <- .getCpUnsLocRepairBwGrid(bwUnsGrid)

  # Preliminary densities use each sample's own adaptive bandwidth curve. These
  # are only weighting curves for constructing the shared bandwidth curve.
  densStimPre <- .getCpUnsLocDensityAdaptiveGrid(
    x = xStim,
    grid = grid,
    bwGrid = bwStimGrid,
    normalise = TRUE,
    probGMin = attr(exTblStimThreshold, "probGMin")
  )

  densUnsPre <- .getCpUnsLocDensityAdaptiveGrid(
    x = xUns,
    grid = grid,
    bwGrid = bwUnsGrid,
    normalise = TRUE,
    probGMin = attr(exTblUnsThreshold, "probGMin")
  )

  if (is.null(densStimPre) || is.null(densUnsPre)) {
    return(NULL)
  }

  denom <- densStimPre$y + densUnsPre$y
  bwShared <- ifelse(
    is.finite(denom) & denom > 0,
    (densStimPre$y * bwStimGrid + densUnsPre$y * bwUnsGrid) / denom,
    rowMeans(cbind(bwStimGrid, bwUnsGrid), na.rm = TRUE)
  )
  bwShared <- .getCpUnsLocRepairBwGrid(bwShared)

  # Final densities use the same grid and the same location-specific bandwidth
  # vector, then are normalised separately over the full padded grid.
  densStim <- .getCpUnsLocDensityAdaptiveGrid(
    x = xStim,
    grid = grid,
    bwGrid = bwShared,
    normalise = TRUE,
    probGMin = attr(exTblStimThreshold, "probGMin")
  )

  densUns <- .getCpUnsLocDensityAdaptiveGrid(
    x = xUns,
    grid = grid,
    bwGrid = bwShared,
    normalise = TRUE,
    probGMin = attr(exTblUnsThreshold, "probGMin")
  )

  if (is.null(densStim) || is.null(densUns)) {
    return(NULL)
  }

  list(
    stim = densStim,
    uns = densUns,
    bw = list(
      adaptive = TRUE,
      grid = grid,
      stim = bwStimObj,
      uns = bwUnsObj,
      stimGrid = bwStimGrid,
      unsGrid = bwUnsGrid,
      sharedGrid = bwShared,
      densStimWeight = densStimPre$y,
      densUnsWeight = densUnsPre$y
    )
  )
}

#' @keywords internal
.getCpUnsLocGetDensRawDensitiesBwAdaptiveOne <- function(
  x,
  chnlSettings
) {
  .bwCalcOne(
    x = x,
    bwMtd = chnlSettings$bwMtd %||% "hpi1Norm",
    bwAdj = chnlSettings$bwAdj %||% 1,
    bwNcellMin = chnlSettings$bwNcellMin,
    bwNcellMax = chnlSettings$bwNcellMax,
    normPeakFrac = chnlSettings$normPeakFrac %||% 0.1,
    normPeakMinRel = chnlSettings$normPeakMinRel %||% 0.75,
    normExtraFrac = chnlSettings$normExtraFrac %||% 0.2,
    normExtraMax = chnlSettings$normExtraMax %||% Inf,
    normExtraJitterFrac = chnlSettings$normExtraJitterFrac %||% 0.25,
    normLambda = chnlSettings$normLambda %||% seq(-2, 2, length.out = 81),
    normDensityN = chnlSettings$normDensityN %||% 512L,
    normExcessBwMtd = chnlSettings$normExcessBwMtd %||% "hpi3",
    normExcessNcell = chnlSettings$normExcessNcell %||% 10000L,
    normAdaptiveNcell = chnlSettings$normAdaptiveNcell %||%
      chnlSettings$bwAdaptiveNcell %||%
      2500L,
    bwAdaptiveCore = chnlSettings$bwAdaptiveCore,
    bwAdaptiveExtra = chnlSettings$bwAdaptiveExtra,
    bwAdaptiveCrossover = chnlSettings$bwAdaptiveCrossover,
    bwAdaptiveTransitionWidth = chnlSettings$bwAdaptiveTransitionWidth %||% 0,
    normMtd = chnlSettings$normMtd %||% "moments",
    adaptive = TRUE
  )
}

#' @keywords internal
.getCpUnsLocAdaptiveGrid <- function(
  x,
  n = 512L,
  padFrac = 0.15
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(numeric(0L))
  }

  n <- .bwAsSafeSampleN(n, default = 512L, lower = 16L)
  rangeVec <- range(x, na.rm = TRUE)
  rangeWidth <- diff(rangeVec)

  if (!is.finite(rangeWidth) || rangeWidth <= 0) {
    rangeWidth <- .bwRobustSd(x)
  }
  if (!is.finite(rangeWidth) || rangeWidth <= 0) {
    rangeWidth <- 1
  }

  padFrac <- suppressWarnings(as.numeric(padFrac)[1])
  if (!is.finite(padFrac) || padFrac < 0) {
    padFrac <- 0.15
  }

  pad <- padFrac * rangeWidth
  seq(
    from = rangeVec[[1]] - pad,
    to = rangeVec[[2]] + pad,
    length.out = n
  )
}

#' @keywords internal
.getCpUnsLocAdaptiveBwOnGrid <- function(
  bwObj,
  grid,
  fallback = NULL
) {
  grid <- suppressWarnings(as.numeric(grid))
  grid <- grid[is.finite(grid)]

  if (
    is.list(bwObj) &&
      all(c("bin", "bw") %in% names(bwObj)) &&
      length(bwObj$bin) >= 2L &&
      length(bwObj$bin) == length(bwObj$bw)
  ) {
    bwGrid <- stats::approx(
      x = suppressWarnings(as.numeric(bwObj$bin)),
      y = suppressWarnings(as.numeric(bwObj$bw)),
      xout = grid,
      rule = 2
    )$y
    return(.getCpUnsLocRepairBwGrid(bwGrid))
  }

  bwScalar <- suppressWarnings(as.numeric(bwObj)[1])
  if (!is.finite(bwScalar) || bwScalar <= 0) {
    bwScalar <- suppressWarnings(as.numeric(fallback)[1])
  }
  if (!is.finite(bwScalar) || bwScalar <= 0) {
    bwScalar <- NA_real_
  }

  rep(bwScalar, length(grid))
}

#' @keywords internal
.getCpUnsLocRepairBwGrid <- function(bwGrid) {
  bwGrid <- suppressWarnings(as.numeric(bwGrid))

  if (length(bwGrid) == 0L) {
    return(bwGrid)
  }

  good <- is.finite(bwGrid) & bwGrid > 0
  if (!any(good)) {
    return(rep(NA_real_, length(bwGrid)))
  }

  bwGrid[!good] <- stats::median(bwGrid[good], na.rm = TRUE)
  pmax(bwGrid, .Machine$double.eps)
}

#' @keywords internal
.getCpUnsLocDensityAdaptiveGrid <- function(
  x,
  grid,
  bwGrid,
  normalise = TRUE,
  probGMin = NULL
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  grid <- suppressWarnings(as.numeric(grid))
  bwGrid <- suppressWarnings(as.numeric(bwGrid))

  ok <- is.finite(grid) & is.finite(bwGrid) & bwGrid > 0
  grid <- grid[ok]
  bwGrid <- bwGrid[ok]

  if (
    length(x) < 2L ||
      length(unique(x)) < 2L ||
      length(grid) < 2L ||
      length(grid) != length(bwGrid)
  ) {
    return(NULL)
  }

  y <- vapply(
    seq_along(grid),
    function(i) {
      mean(stats::dnorm(
        x = grid[[i]],
        mean = x,
        sd = bwGrid[[i]]
      ))
    },
    numeric(1)
  )

  y <- pmax(y, 0)

  y <- .getCpUnsLocDensityNormalizeAndScale(
    grid = grid,
    y = y,
    normalise = normalise,
    probGMin = probGMin
  )

  out <- list(
    x = grid,
    y = y,
    bw = bwGrid,
    n = length(x),
    call = match.call(),
    data.name = deparse(substitute(x)),
    has.na = FALSE
  )
  class(out) <- "density"
  out
}

#' @keywords internal
.getCpUnsLocTrapz <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2L || length(x) != length(y)) {
    return(NA_real_)
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

#' @keywords internal
.getCpUnsLocDensityNormalizeAndScale <- function(
  grid,
  y,
  normalise = TRUE,
  probGMin = NULL
) {
  grid <- suppressWarnings(as.numeric(grid))
  y <- suppressWarnings(as.numeric(y))

  if (length(grid) != length(y)) {
    return(y)
  }

  y <- pmax(y, 0)

  if (isTRUE(normalise)) {
    area <- .getCpUnsLocTrapz(grid, y)
    if (is.finite(area) && area > 0) {
      y <- y / area
    }
  }

  probGMin <- suppressWarnings(as.numeric(probGMin)[1])
  if (is.finite(probGMin)) {
    probGMin <- max(0, min(1, probGMin))
    y <- y * probGMin
  }

  y
}


#' @keywords internal
.getCpUnsLocGetDensRawDensitiesBw <- function(
  exTblStimThreshold,
  exTblUnsThreshold,
  bw,
  bwMin,
  bwMax,
  bwFallback,
  bwMtd,
  bwAdj,
  bwNcellMin,
  bwNcellMax,
  chnlSettings = NULL
) {
  if (!is.null(bw)) {
    return(bw)
  }
  bwStim <- .getCpUnsLocGetDensRawDensitiesBwInit(
    .data = .getCut(exTblStimThreshold),
    bwMin = bwMin,
    bwMax = bwMax,
    bwFallback = bwFallback,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax,
    chnlSettings = chnlSettings
  )
  bwUns <- .getCpUnsLocGetDensRawDensitiesBwInit(
    .data = .getCut(exTblUnsThreshold),
    bwMin = bwMin,
    bwMax = bwMax,
    bwFallback = bwFallback,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax,
    chnlSettings = chnlSettings
  )
  min(bwUns, bwStim)
}


#' @keywords internal
.getCpUnsLocGetDensRawDensitiesStim <- function(
  exTblStimThreshold,
  bw
) {
  densObj <- density(.getCut(exTblStimThreshold), bw = bw)
  if (is.null(attr(exTblStimThreshold, "probGMin"))) {
    return(densObj)
  }
  densObj$y <- densObj$y * attr(exTblStimThreshold, "probGMin")
  densObj
}

#' @keywords internal
.getCpUnsLocGetDensRawDensitiesUns <- function(
  exTblUnsThreshold,
  densStim,
  bw
) {
  densObj <- density(
    .getCut(exTblUnsThreshold),
    from = min(densStim$x),
    to = max(densStim$x),
    bw = bw
  )
  if (is.null(attr(exTblUnsThreshold, "probGMin"))) {
    return(densObj)
  }
  densObj$y <- densObj$y * attr(exTblUnsThreshold, "probGMin")
  densObj
}

#' @keywords internal
.getCpUnsLocGetDensRawTabulate <- function(
  stimX,
  stimY,
  unsX,
  unsY
) {
  densTblRawStim <- tibble::tibble(xStim = stimX, yStim = stimY)
  densTblRawWide <- .getCpUnsLocGetDensRawTabulateUnsInterp(
    .data = densTblRawStim,
    unsX = unsX,
    unsY = unsY
  )
  .getCpUnsLocGetDensRawTabulateFormat(densTblRawWide)
}

#' @keywords internal
.getCpUnsLocGetDensRawTabulateUnsInterp <- function(
  .data,
  unsX,
  unsY
) {
  yUns <- stats::approx(
    x = suppressWarnings(as.numeric(unsX)),
    y = suppressWarnings(as.numeric(unsY)),
    xout = suppressWarnings(as.numeric(.data$xStim)),
    rule = 2
  )$y

  .data |>
    dplyr::mutate(yUns = .env$yUns)
}


#' @keywords internal
.getCpUnsLocGetDensRawTabulateFormat <- function(.data) {
  .data |>
    tidyr::pivot_longer(yStim:yUns, names_to = "stim", values_to = "dens") |> # nolint
    dplyr::mutate(stim = ifelse(stim == "yStim", "yes", "no")) # nolint
}

# get probabilities
# -------------------
#' @keywords internal
.getCpUnsLocGetProbTbl <- function(
  densTblRaw,
  stage,
  cpMin,
  exVecStimThreshold,
  exVecUnsThreshold
) {
  .debug("Normalising probabilities") # nolint

  probTbl <- .getCpUnsLocGetProbTblInit(densTblRaw, cpMin)

  densTblStim <- densTblRaw |>
    dplyr::filter(.data$stim == "yes") |>
    dplyr::arrange(.data$xStim)
  densTblUns <- densTblRaw |>
    dplyr::filter(.data$stim == "no") |>
    dplyr::arrange(.data$xStim)
  peakStimIdx <- .getPeakMainLeftIdx(densTblStim$dens)
  peakUnsIdx <- .getPeakMainLeftIdx(densTblUns$dens)
  peakStimX <- densTblStim$xStim[peakStimIdx]
  peakUnsX <- densTblUns$xStim[peakUnsIdx]

  probTblPos <- .getCpUnsLocProbTblFilter(
    densTbl = densTblRaw,
    probTbl = probTbl,
    exVecStim = exVecStimThreshold,
    exVecUns = exVecUnsThreshold,
    stage = stage,
    peakStimX = peakStimX,
    peakUnsX = peakUnsX
  )

  list(
    all = probTbl,
    pos = probTblPos,
    densityBw = attr(densTblRaw, "locDensityBw"),
    stimDensity = densTblStim |>
      dplyr::transmute(
        x = suppressWarnings(as.numeric(.data$xStim)),
        y = suppressWarnings(as.numeric(.data$dens))
      ),
    stimPeakX = suppressWarnings(as.numeric(peakStimX)[1L])
  )
}

#' @keywords internal
.getCpUnsLocGetProbTblInit <- function(densTblRaw, cpMin) {
  densTblRaw |>
    tidyr::pivot_wider(
      id_cols = xStim,
      names_from = stim,
      values_from = dens # nolint
    ) |>
    dplyr::mutate(
      probStim = 1 - no / yes, # nolint
      probStim = ifelse(yes == 0 & no == 0, 0, probStim), # nolint
      probStimNorm = pmin(1, probStim), # nolint
      probStimNorm = pmax(0, probStimNorm) # nolint
    ) |>
    dplyr::filter(xStim > cpMin)
}

#' @keywords internal
.getCpUnsLocGetDensRawDensitiesBwInit <- function(
  .data,
  bwMin,
  bwMax,
  bwFallback,
  bwMtd,
  bwAdj,
  bwNcellMin,
  bwNcellMax,
  chnlSettings = NULL
) {
  chnlSettings <- chnlSettings %||% list()
  .data <- suppressWarnings(as.numeric(.data))
  .data <- .data[is.finite(.data)]

  if (length(.data) < 2L || length(unique(.data)) < 2L) {
    return(bwFallback)
  }

  bwCalc <- .bwCalcOne(
    x = .data,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax,
    normPeakFrac = chnlSettings$normPeakFrac %||% 0.1,
    normPeakMinRel = chnlSettings$normPeakMinRel %||% 0.75,
    normExtraFrac = chnlSettings$normExtraFrac %||% 0.2,
    normExtraMax = chnlSettings$normExtraMax %||% Inf,
    normExtraJitterFrac = chnlSettings$normExtraJitterFrac %||% 0.25,
    normLambda = chnlSettings$normLambda %||% seq(-2, 2, length.out = 81),
    normDensityN = chnlSettings$normDensityN %||% 512L,
    normExcessBwMtd = chnlSettings$normExcessBwMtd %||% "hpi3",
    normExcessNcell = chnlSettings$normExcessNcell %||% 10000L,
    normAdaptiveNcell = chnlSettings$normAdaptiveNcell %||%
      chnlSettings$bwAdaptiveNcell %||%
      2500L,
    bwAdaptiveCore = chnlSettings$bwAdaptiveCore,
    bwAdaptiveExtra = chnlSettings$bwAdaptiveExtra,
    bwAdaptiveCrossover = chnlSettings$bwAdaptiveCrossover,
    bwAdaptiveTransitionWidth = chnlSettings$bwAdaptiveTransitionWidth %||% 0,
    normMtd = chnlSettings$normMtd %||% "moments",
    adaptive = FALSE
  )

  if (!is.finite(bwCalc) || bwCalc <= 0) {
    return(bwFallback)
  }

  max(bwMin, min(as.numeric(bwCalc)[1], bwMax))
}


#' @keywords internal
.getCpUnsLocProbTblFilter <- function(
  densTbl,
  probTbl,
  exVecStim,
  exVecUns,
  stage,
  peakStimX = NULL,
  peakUnsX = NULL
) {
  .debug("Filtering before smoothing") # nolint
  densTblStim <- densTbl |>
    dplyr::filter(stim == "yes")
  densTblUns <- densTbl |>
    dplyr::filter(stim == "no")

  if (is.null(peakStimX) || !is.finite(peakStimX)) {
    peakStimIdx <- .getPeakMainLeftIdx(densTblStim$dens)
    peakStimX <- densTblStim$xStim[peakStimIdx]
  }
  if (is.null(peakUnsX) || !is.finite(peakUnsX)) {
    peakUnsIdx <- .getPeakMainLeftIdx(densTblUns$dens)
    peakUnsX <- densTblUns$xStim[peakUnsIdx]
  }
  peakX <- max(peakStimX, peakUnsX)

  windowWidthStim <- 1 /
    3 *
    abs(diff(quantile(exVecStim[exVecStim < peakStimX], c(0.05, 1))))
  windowWidthUns <- 1 /
    3 *
    abs(diff(quantile(exVecUns[exVecUns < peakUnsX], c(0.05, 1))))
  windowWidth <- max(windowWidthStim, windowWidthUns, na.rm = TRUE)

  probTbl <- probTbl |>
    dplyr::filter(xStim > peakX + windowWidth) # nolint

  if (nrow(probTbl) <= 5L) {
    return(probTbl)
  }

  # Find the threshold index
  valid_idx <- probTbl |>
    dplyr::arrange(xStim) |>
    dplyr::mutate(
      ge_0025 = probStimNorm >= 0.025,
      ge_0075 = probStimNorm >= 0.075,
      n_remaining = rev(seq_len(dplyr::n())),

      # Calculate tail proportions
      prop_0025 = rev(cumsum(rev(ge_0025))) / n_remaining,
      prop_0075 = rev(cumsum(rev(ge_0075))) / n_remaining,

      # Do both conditions match?
      both_met = prop_0025 >= 0.90 & prop_0075 >= 0.25
    ) |>
    # Get the row index of the first time this becomes true
    dplyr::pull(both_met) |>
    which() |>
    which.min()

  # Slice the table from that first valid point all the way to the end
  if (length(valid_idx) > 0) {
    probTbl |>
      dplyr::arrange(xStim) |>
      dplyr::slice(valid_idx:dplyr::n())
  } else {
    probTbl |> dplyr::slice(0) # Empty if nothing matches
  }
}

#' @keywords internal
.getCpUnsLocGetMinProbX <- function(probTblPos) {
  min(probTblPos$xStim)
}

#' @keywords internal
.getCpUnsLocCheckResponse <- function(probTblPos, exTblStimOrig) {
  nrow(probTblPos) == 0 ||
    max(.getCut(exTblStimOrig)) < .getCpUnsLocGetMinProbX(probTblPos)
}

#' @keywords internal
.getCpUnsLocGetDataMod <- function(
  exTblStimThreshold,
  exTblStimNoMin,
  exTblUnsThreshold,
  exTblUnsBias,
  probTblList,
  cpMin,
  stage
) {
  if (.getCpUnsLocCheckResponse(probTblList$pos, exTblStimNoMin)) {
    return(.getCpUnsLocConditionCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "No responding cells"
    ))
  }

  minProbXPos <- min(probTblList$pos$xStim, na.rm = TRUE)

  margin <- .getCpUnsLocGetDataModMargin(
    exTblStimNoMin = exTblStimNoMin,
    exTblUnsNoMin = exTblUnsThreshold
  )

  minMod <- minProbXPos - margin

  # Keep the lower values in dataMod to anchor/clamp the spline at the periphery
  dataMod <- exTblStimThreshold[
    .getCut(exTblStimThreshold) >= minMod,
    ,
    drop = FALSE
  ]

  # Interpolate using the 'all' table so baseline cells get real probabilities
  probVec <- try(
    approx(
      x = probTblList$all$xStim,
      y = probTblList$all$probStimNorm,
      xout = .getCut(dataMod),
      method = "linear",
      rule = 2
    )$y,
    silent = TRUE
  )

  if (inherits(probVec, "try-error")) {
    return(.getCpUnsLocConditionCheckOut(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      stage = stage,
      msg = "No responding cells"
    ))
  }

  dataMod <- dataMod |>
    dplyr::mutate(probSmooth = probVec)

  # get the bins that we'll need for thinning
  binVec <- .getCpUnsLocGetDataModBinVec(
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold
  )
  attr(dataMod, "binVec") <- binVec

  # Attach this for filtering when it comes to estimating
  # the final response proportion
  attr(dataMod, "minProbXPos") <- minProbXPos

  # Retain the exact bandwidth object used for the stimulated and unstimulated
  # density estimates. The later antimode density uses half this bandwidth.
  attr(dataMod, "locDensityBw") <- probTblList$densityBw

  # Retain the original stimulated density and its already identified
  # left-main peak for the post-smoothing marginal-density safeguard.
  attr(dataMod, "locStimDensity") <- probTblList$stimDensity
  attr(dataMod, "locStimPeakX") <- probTblList$stimPeakX

  .thinDataMod(dataMod, maxCellsPerBin = 20)
}

#' @keywords internal
.getCpUnsLocGetDataModMargin <- function(
  exTblStimNoMin,
  exTblUnsNoMin
) {
  spanStim <- diff(quantile(
    .getCut(exTblStimNoMin),
    probs = c(0.05, 0.95),
    na.rm = TRUE
  ))
  spanUns <- diff(quantile(
    .getCut(exTblUnsNoMin),
    probs = c(0.05, 0.95),
    na.rm = TRUE
  ))

  max(spanStim, spanUns) * 0.05
}

.getCpUnsLocGetDataModBinVec <- function(
  exTblStimThreshold,
  exTblUnsThreshold
) {
  stimVals <- .getCut(exTblStimThreshold)
  unsVals <- .getCut(exTblUnsThreshold)

  rng <- range(stimVals, unsVals, na.rm = TRUE)

  seq.int(from = rng[1], to = rng[2], length.out = 512L)
}

# smooth
# ---------------------
#' @keywords internal

.thinDataMod <- function(dataMod, maxCellsPerBin = 20) {
  binVec <- attr(dataMod, "binVec")
  if (is.null(binVec)) {
    return(dataMod)
  }

  breaks <- c(-Inf, binVec[-length(binVec)] + diff(binVec) / 2, Inf)

  x_vals <- .getCut(dataMod)
  bin_indices <- findInterval(x_vals, breaks)

  # Sample down if a bin has too many cells
  idxMod <- unlist(lapply(split(seq_along(x_vals), bin_indices), function(idx) {
    if (length(idx) > maxCellsPerBin) {
      sample(idx, size = maxCellsPerBin, replace = FALSE)
    } else {
      idx
    }
  }))

  attr(dataMod, "idxMod") <- sort(unname(idxMod))

  dataMod
}
