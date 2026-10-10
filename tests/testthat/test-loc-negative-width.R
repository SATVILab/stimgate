# These boundaries intentionally replace the old 5th-percentile span and
# width / 3 offset: distant negative modes must not set the response-search start.
.filterNegWidth <- function(densityStim, densityUns = densityStim,
                            peakStim = 0, peakUns = 0, bw = 0.1,
                            shiftedPeak = NULL, x = seq(0, 5, by = 0.01)) {
  stimgate:::.getCpUnsLocProbTblFilter(
    probTbl = data.frame(xStim = x, probStimNorm = 1),
    exVecStim = densityStim$x, exVecUns = densityUns$x, stage = "test",
    peakStimX = peakStim, peakUnsX = peakUns, shiftedPeak = shiftedPeak,
    densityStim = densityStim, densityUns = densityUns, densityBw = bw
  )
}

test_that("a normal negative width is its density half-height distance", {
  # Deterministic normal quantiles avoid sample-noise tolerance and RNG state.
  ex <- stats::qnorm(seq(0.0001, 0.9999, length.out = 10001))
  bw <- 0.1
  dens <- stats::density(ex, bw = bw, n = 4096)
  peak <- dens$x[which.max(dens$y)]
  width <- stimgate:::.getCpUnsLocNegWidth(dens, peak, ex, bw)
  expect_identical(width$source, "half_height")
  expect_equal(width$width, sqrt(2 * log(2)), tolerance = 0.02)
  expect_true(is.na(width$dipX))
  filtered <- .filterNegWidth(dens, peakStim = peak, peakUns = peak, bw = bw)
  start <- peak + 0.5 * width$width
  expect_equal(filtered$windowWidthInfo$searchStartX, start)
  expect_true(all(filtered$probTbl$xStim > start))
  expect_lte(min(filtered$probTbl$xStim) - start, 0.01)
})

test_that("a deep second negative mode does not inflate the main width", {
  x <- seq(-8, 3, length.out = 11001)
  dens <- list(x = x, y = 0.7 * stats::dnorm(x, 0, 0.4) +
    0.3 * stats::dnorm(x, -5, 0.4))
  ex <- c(stats::qnorm(seq(0.001, 0.999, length.out = 7000), 0, 0.4),
    stats::qnorm(seq(0.001, 0.999, length.out = 3000), -5, 0.4))
  width <- stimgate:::.getCpUnsLocNegWidth(dens, 0, ex, 0.1)
  oldWidth <- abs(diff(stats::quantile(ex[ex < 0], c(0.05, 1))))
  expect_identical(width$source, "half_height")
  expect_equal(width$width, 0.4 * sqrt(2 * log(2)), tolerance = 0.002)
  expect_gt(width$halfHeightX, width$dipX)
  expect_lt(width$width, oldWidth / 5)
})

test_that("clear dips win but shallow or too-close dips are ignored", {
  # The first half-height crossing is -3; the trough at -1 is closer.
  dens <- list(x = -4:1, y = c(0.2, 0.5, 0.9, 0.6, 1, 0.8))
  width <- stimgate:::.getCpUnsLocNegWidth(dens, 0, dens$x, 1)
  expect_identical(width$source, "dip")
  expect_equal(width$width, 1)
  expect_equal(width$halfHeightX, -3)
  expect_equal(width$dipX, -1)

  boundary <- dens
  boundary$y[4] <- 0.75
  expect_identical(
    stimgate:::.getCpUnsLocNegWidth(boundary, 0, dens$x, 1)$source, "dip"
  )

  shallow <- dens
  shallow$y[4] <- 0.8
  ignored <- stimgate:::.getCpUnsLocNegWidth(shallow, 0, dens$x, 1)
  expect_identical(ignored$source, "half_height")
  expect_true(is.na(ignored$dipX))
  close <- stimgate:::.getCpUnsLocNegWidth(dens, 0, dens$x, 1.01)
  expect_identical(close$source, "half_height")
  expect_true(is.na(close$dipX))
})

test_that("crossings interpolate and the closest qualifying dip is used", {
  dens <- list(x = -6:1, y = c(0.2, 0.4, 0.9, 0.65, 0.9, 0.7, 1, 0.8))
  width <- stimgate:::.getCpUnsLocNegWidth(dens, 0, dens$x, 0.5)
  expect_identical(width$source, "dip")
  expect_equal(width$dipX, -1)
  expect_equal(width$halfHeightX, -4.8)
  # A flat trough counts at its right edge, not a shoulder on a descending run.
  flat <- list(x = -5:1, y = c(0.2, 0.4, 0.9, 0.7, 0.7, 1, 0.8))
  expect_equal(stimgate:::.getCpUnsLocNegWidth(flat, 0, flat$x, 0.5)$dipX, -1)
})

test_that("missing candidates record the tube minimum fallback", {
  dens <- list(x = -2:1, y = c(0.8, 0.9, 1, 0.9))
  width <- stimgate:::.getCpUnsLocNegWidth(dens, 0, c(-7, -1, 0, NA), 0.2)
  expect_identical(width$source, "fallback")
  expect_equal(width$width, 7)
  expect_equal(width$boundaryX, -7)
  expect_true(is.na(width$halfHeightX))
  expect_true(is.na(width$dipX))
})

test_that("the response-search offset floors at one bandwidth", {
  dens <- list(x = c(-0.1, -0.02, 0, 0.02), y = c(0, 0.5, 1, 0.5))
  filtered <- .filterNegWidth(dens, bw = 0.2, x = c(0.01, 0.1, 0.2, 0.21, 1))
  expect_equal(filtered$windowWidth, 0.02)
  expect_equal(filtered$windowWidthInfo$searchStartX, 0.2)
  expect_equal(filtered$probTbl$xStim, c(0.21, 1))
})

test_that("tube widths combine by max and shifted peaks use only the uns width", {
  x <- seq(-3, 5, by = 0.01)
  uns <- list(x = x, y = stats::dnorm(x, 0, 0.3))
  stim <- list(x = x, y = stats::dnorm(x, 2, 0.8))
  ordinary <- .filterNegWidth(stim, uns, peakStim = 2)
  shifted <- .filterNegWidth(stim, uns, peakStim = 2,
    shiftedPeak = list(mult = 2, bw = 0.1, bwSource = "shared"))
  expect_equal(ordinary$windowWidth, ordinary$windowWidthInfo$stim$width)
  expect_true(shifted$shiftedPeakRef$applied)
  expect_equal(shifted$peakX, 0)
  expect_equal(shifted$windowWidth, 0.3 * sqrt(2 * log(2)), tolerance = 0.001)
  expect_equal(shifted$windowWidth, shifted$windowWidthInfo$uns$width)
  expect_equal(shifted$shiftedPeakRef$windowWidthUns, shifted$windowWidth)
  expect_equal(shifted$windowWidthInfo$searchStartX, 0.5 * shifted$windowWidth)
  expect_equal(shifted$shiftedPeakRef$windowWidthUnsInfo, shifted$windowWidthInfo$uns)
  # A widened density bandwidth floors the start without changing the trigger.
  widened <- .filterNegWidth(stim, uns, peakStim = 2, bw = 0.6,
    shiftedPeak = list(mult = 2, bw = 0.1, bwSource = "shared"))
  expect_true(widened$shiftedPeakRef$applied)
  expect_equal(widened$windowWidthInfo$searchStartX, 0.6)
})

test_that("adaptive widths and offsets use the shared bandwidth at each peak", {
  dens <- list(x = -4:1, y = c(0.2, 0.5, 0.9, 0.6, 1, 0.8))
  adaptive <- list(adaptive = TRUE, grid = c(-1, 1), sharedGrid = c(2, 4))
  filtered <- .filterNegWidth(dens, bw = adaptive, x = c(1, 2, 3, 4))
  expect_equal(filtered$windowWidthInfo$stim$bw, 3)
  expect_true(is.na(filtered$windowWidthInfo$stim$dipX))
  expect_equal(filtered$windowWidth, 3)
  expect_equal(filtered$windowWidthInfo$searchStartX, 3)
  expect_equal(filtered$probTbl$xStim, 4)
})

test_that("width diagnostics survive derivative row subsetting", {
  data <- data.frame(val = 1:3)
  attr(data, "locWindowWidth") <- 1
  info <- list(stim = list(source = "fallback", boundaryX = -2))
  attr(data, "locWindowWidthInfo") <- info
  subset <- stimgate:::.getCpUnsLocSubsetRows(data, c(TRUE, FALSE, TRUE))
  expect_equal(attr(subset, "locWindowWidth"), 1)
  expect_equal(attr(subset, "locWindowWidthInfo"), info)
})

test_that("cytokine-positive refinement uses the same negative-width definition", {
  x <- stats::qnorm(seq(0.0001, 0.9999, length.out = 2001))
  ref <- stimgate:::.getCytPosMarginalReference(data.frame(cyt = x), "cyt")
  expect_identical(ref$reason, "marginal_reference_available")
  expect_identical(ref$windowWidthInfo$source, "half_height")
  expect_equal(ref$windowWidth, sqrt(2 * log(2)), tolerance = 0.06)
  expect_equal(ref$lowerX, ref$peakX + max(0.5 * ref$windowWidth, ref$densityBw))
})

test_that("shape-tailgate margins use the recorded negative width", {
  testthat::local_mocked_bindings(
    .getCpUnsLocMarginalDensityLowerBound = function(...) list(lowerBoundX = 2),
    .getCpUnsLocAntimodeDensity = function(...) NULL,
    .package = "stimgate"
  )
  info <- list(stim = list(source = "half_height", width = 0.8))
  ex <- data.frame(val = 1:10)
  attr(ex, "chnlCut") <- "val"
  shape <- stimgate:::.getCpUnsLocGetShapeThreshold(
    ex, list(peakX = 0, windowWidth = 0.8, windowWidthInfo = info,
      stimDensity = data.frame(x = 0:2, y = c(1, 0.5, 0))),
    list(locShapeTailgateMarginFrac = 0.25)
  )
  expect_equal(shape$info$adjustedTailgateX, 2.2)
  expect_equal(shape$info$windowWidthInfo, info)
})

test_that("shape refits retain the ordinary width diagnostics", {
  testthat::local_mocked_bindings(
    .getCpUnsLocGetProbSmooth = function(dataMod, ...) dataMod,
    .package = "stimgate"
  )
  ex <- data.frame(val = stats::qnorm(seq(0.001, 0.999, length.out = 301)))
  attr(ex, "chnlCut") <- "val"
  info <- list(stim = list(source = "dip", width = 0.6), searchStartX = 0.3)
  fit <- stimgate:::.getCpUnsLocGetProbFit(
    exTblStimNoMin = ex, exTblStimThreshold = ex,
    exTblUnsThreshold = ex, exTblUnsBias = ex, bias = 0,
    exTblUnsOrig = ex, stage = "test", pathProject = tempdir(),
    chnlSettings = list(bw = 0.2, cpMin = -Inf),
    applyPreliminaryFilter = FALSE, peakX = 0, windowWidth = 0.6,
    windowWidthInfo = info
  )
  expect_equal(fit$probTblList$windowWidth, 0.6)
  expect_equal(fit$probTblList$windowWidthInfo, info)
  expect_equal(attr(fit$dataMod, "locWindowWidth"), 0.6)
  expect_equal(attr(fit$dataMod, "locWindowWidthInfo"), info)
})
