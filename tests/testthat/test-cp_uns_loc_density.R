# Focused tests for the ordinary fixed-bandwidth local-FDR density / raw
# probability stage.

.makeCpUnsLocDensityFixture <- function(seed = 42L) {
  set.seed(seed)
  exStim <- c(
    stats::rnorm(350, mean = 0, sd = 1),
    stats::rnorm(170, mean = 4.2, sd = 0.6)
  )
  exUns <- stats::rnorm(500, mean = 0, sd = 1)

  exTblStim <- data.frame(val = exStim)
  attr(exTblStim, "chnlCut") <- "val"

  exTblUns <- data.frame(val = exUns)
  attr(exTblUns, "chnlCut") <- "val"

  list(
    stim = exTblStim,
    uns = exTblUns,
    stimVec = exStim,
    unsVec = exUns
  )
}

.makeCpUnsLocDensityChnlSettings <- function(bw = 0.55) {
  list(
    bw = bw,
    bwMin = 0.05,
    bwMax = 2,
    bwFallback = 0.2,
    bwMtd = "nrd0",
    bwAdj = 1,
    bwNcellMin = 20L,
    bwNcellMax = 5000L,
    cpMin = -Inf,
    locEnforceShapeThreshold = FALSE
  )
}

# Basic density alignment and bandwidth propagation are critical for the current
# fixed-bandwidth path, because downstream filtering expects a single comparison
# grid and a consistent bandwidth object paired with the raw densities.
test_that("ordinary fixed-bandwidth densities stay finite, aligned and retain the selected bandwidth", {
  withr::local_preserve_seed()
  fixture <- .makeCpUnsLocDensityFixture(seed = 101L)
  chnlSettings <- .makeCpUnsLocDensityChnlSettings(bw = 0.55)

  densRaw <- stimgate:::.getCpUnsLocGetDensRaw(
    exTblStimThreshold = fixture$stim,
    exTblUnsThreshold = fixture$uns,
    stage = "test",
    pathProject = tempdir(),
    chnlSettings = chnlSettings
  )

  expect_true(is.data.frame(densRaw))
  expect_true(all(is.finite(densRaw$dens)))
  expect_true(all(densRaw$dens >= 0))
  expect_equal(
    densRaw |> dplyr::filter(.data$stim == "yes") |> dplyr::pull(.data$xStim),
    densRaw |> dplyr::filter(.data$stim == "no") |> dplyr::pull(.data$xStim),
    tolerance = 1e-8
  )
  expect_equal(as.numeric(attr(densRaw, "locDensityBw")), chnlSettings$bw)
})

# The raw response probabilities from the fixed-bandwidth route must stay inside
# [0, 1] and increase in the clearly separated positive region, as a guardrail
# against negative or mis-scaled densities being fed into the smoother.
test_that("ordinary fixed-bandwidth raw probabilities remain valid and rise in the response region", {
  withr::local_preserve_seed()
  fixture <- .makeCpUnsLocDensityFixture(seed = 202L)
  chnlSettings <- .makeCpUnsLocDensityChnlSettings(bw = 0.6)

  densRaw <- stimgate:::.getCpUnsLocGetDensRaw(
    exTblStimThreshold = fixture$stim,
    exTblUnsThreshold = fixture$uns,
    stage = "test",
    pathProject = tempdir(),
    chnlSettings = chnlSettings
  )

  probTbl <- stimgate:::.getCpUnsLocGetProbTbl(
    densTblRaw = densRaw,
    stage = "test",
    cpMin = -Inf,
    exVecStimThreshold = fixture$stimVec,
    exVecUnsThreshold = fixture$unsVec
  )

  allProb <- probTbl[["all"]]
  expect_true(all(is.finite(allProb$probStimNorm)))
  expect_true(all(allProb$probStimNorm >= 0 & allProb$probStimNorm <= 1))

  negRegion <- allProb$probStimNorm[
    allProb$xStim > -1 & allProb$xStim < 1
  ]
  responseRegion <- allProb$probStimNorm[
    allProb$xStim > 2.5 & allProb$xStim < 6
  ]

  expect_true(length(responseRegion) > 0)
  expect_true(length(negRegion) > 0)
  expect_gt(mean(responseRegion, na.rm = TRUE), mean(negRegion, na.rm = TRUE))
})

# The ordinary fit should package the observed density metadata and the smoothable
# probability inputs that the later local-FDR stages expect.
test_that("ordinary fixed-bandwidth fit preserves the density metadata and smoothing inputs", {
  withr::local_preserve_seed()
  fixture <- .makeCpUnsLocDensityFixture(seed = 303L)
  chnlSettings <- .makeCpUnsLocDensityChnlSettings(bw = 0.5)

  fit <- stimgate:::.getCpUnsLocGetProbFit(
    exTblStimNoMin = fixture$stim,
    exTblStimThreshold = fixture$stim,
    exTblUnsThreshold = fixture$uns,
    exTblUnsBias = fixture$uns,
    bias = 0,
    exTblUnsOrig = fixture$uns,
    stage = "test",
    pathProject = tempdir(),
    chnlSettings = chnlSettings,
    applyPreliminaryFilter = TRUE
  )

  dataMod <- fit[["dataMod"]]
  probTblList <- fit[["probTblList"]]

  expect_true(is.data.frame(dataMod))
  expect_true(all(c("val", "probSmooth", "pred") %in% names(dataMod)))
  expect_true(all(is.finite(dataMod$probSmooth)))
  expect_true(all(dataMod$probSmooth >= 0 & dataMod$probSmooth <= 1))
  expect_true(all(is.finite(dataMod$pred)))
  expect_equal(
    as.numeric(attributes(dataMod)$locDensityBw),
    as.numeric(probTblList[["densityBw"]]),
    tolerance = 1e-8
  )
  expect_true(all(c("all", "pos", "densityBw", "peakX", "windowWidth") %in% names(probTblList)))
  expect_true(is.finite(probTblList[["peakX"]]))
  expect_true(is.finite(probTblList[["windowWidth"]]) && probTblList[["windowWidth"]] > 0)
  expect_equal(
    attr(dataMod, "locWindowWidthInfo"), probTblList$windowWidthInfo
  )
})

# A narrow negative peak with a small bandwidth: when the stimulated negatives
# are shifted above every unstimulated cell, the unstimulated peak lies below
# the stimulated density's range. It must still be found, rather than a small
# bump among the positives, or the response region starts inside the positives.
test_that("the unstimulated peak is found when it lies below every stimulated cell", {
  withr::local_preserve_seed()
  set.seed(404L)
  exStim <- c(
    stats::rnorm(800, mean = 1.05, sd = 0.02),
    stats::rnorm(200, mean = 2, sd = 0.3)
  )
  exUnsNeg <- stats::rnorm(960, mean = 0.95, sd = 0.02)
  exUns <- c(exUnsNeg, stats::rnorm(40, mean = 2, sd = 0.3))
  exStim <- exStim[exStim > max(exUnsNeg)]
  exTblStim <- data.frame(val = exStim)
  attr(exTblStim, "chnlCut") <- "val"
  exTblUns <- data.frame(val = exUns)
  attr(exTblUns, "chnlCut") <- "val"
  chnlSettings <- .makeCpUnsLocDensityChnlSettings(bw = 0.01)

  densRaw <- stimgate:::.getCpUnsLocGetDensRaw(
    exTblStimThreshold = exTblStim,
    exTblUnsThreshold = exTblUns,
    stage = "test",
    pathProject = tempdir(),
    chnlSettings = chnlSettings
  )
  expect_gt(min(densRaw$xStim), stats::median(exUns))

  peakX <- stimgate:::.getCpUnsLocGetPeakX(densRaw, exStim, exUns)
  expect_equal(peakX[["uns"]], 0.95, tolerance = 0.03)
  expect_equal(peakX[["stim"]], 1.05, tolerance = 0.03)

  probTbl <- stimgate:::.getCpUnsLocGetProbTbl(
    densTblRaw = densRaw,
    stage = "test",
    cpMin = -Inf,
    exVecStimThreshold = exStim,
    exVecUnsThreshold = exUns
  )
  expect_lt(min(probTbl$pos$xStim), 1.5)
  # The unstimulated width must use its full peak-search density, too.
  expect_identical(probTbl$windowWidthInfo$uns$source, "half_height")
  expect_equal(probTbl$windowWidthInfo$uns$peakX, peakX[["uns"]])
  expect_lt(probTbl$windowWidthInfo$uns$width, 0.06)
})

# When the unstimulated cells lie within the stimulated range, the joint range
# is the existing grid, so the peaks come from the raw densities unchanged.
test_that("peaks come from the raw densities when the unstimulated cells lie inside", {
  withr::local_preserve_seed()
  fixture <- .makeCpUnsLocDensityFixture(seed = 505L)
  exUns <- fixture$unsVec
  exUns <- exUns[exUns >= min(fixture$stimVec) & exUns <= max(fixture$stimVec)]
  exTblUns <- data.frame(val = exUns)
  attr(exTblUns, "chnlCut") <- "val"

  densRaw <- stimgate:::.getCpUnsLocGetDensRaw(
    exTblStimThreshold = fixture$stim,
    exTblUnsThreshold = exTblUns,
    stage = "test",
    pathProject = tempdir(),
    chnlSettings = .makeCpUnsLocDensityChnlSettings(bw = 0.5)
  )
  rawPeak <- function(level) {
    d <- densRaw[densRaw$stim == level, ]
    d <- d[order(d$xStim), ]
    d$xStim[stimgate:::.getPeakMainLeftIdx(d$dens)]
  }

  expect_identical(
    stimgate:::.getCpUnsLocGetPeakX(densRaw, fixture$stimVec, exUns),
    c(stim = rawPeak("yes"), uns = rawPeak("no"))
  )
})
