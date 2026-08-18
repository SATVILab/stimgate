test_that("analysis helpers do not expose or forward calcSinglePosGates", {
  skip_if_not_installed("stimgate")

  root_dir <- testthat::test_path("../..")
  script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
  script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
  script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

  skip_if_not(file.exists(script_bw), "sim-bandwidth.R not found")

  source(script_misc, local = TRUE)
  source(script_bw, local = TRUE)
  source(script_comp, local = TRUE)

  expect_false("calcSinglePosGates" %in% names(formals(.simBandwidthBsFreq)))
  expect_false("calcSinglePosGates" %in% names(formals(.simCompareStimgateRows)))
  expect_false("calcSinglePosGates" %in% names(formals(.simCompareFreqBs)))

  # Test that .simBandwidthBsFreq calls gateStim without unused-argument error
  res <- .simBandwidthBsFreq(
    nSample = 1L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0,
    bwMtd = "hpi1",
    nCellStim = 50L,
    probResponse = 0.1,
    meanPos = 5,
    transformation = "gaussian",
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.1,
    ncellUnsRelativeToStim = 1
  )

  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
})

test_that("direct bandwidth simulation helpers agree numerically with package implementation", {
  skip_if_not_installed("stimgate")

  root_dir <- testthat::test_path("../..")
  script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
  script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
  script_adaptive <- file.path(root_dir, "scripts", "r", "sim-bw-adaptive.R")

  skip_if_not(file.exists(script_bw), "sim-bandwidth.R not found")
  skip_if_not(file.exists(script_adaptive), "sim-bw-adaptive.R not found")

  source(script_misc, local = TRUE)
  source(script_bw, local = TRUE)
  source(script_adaptive, local = TRUE)

  set.seed(42)
  x_sample <- rnorm(200, mean = 5, sd = 1.5)

  # Check standard normalized bandwidth via .simBandwidthBwOne vs .bwCalcOne
  bw_direct_std <- .simBandwidthBwOne(
    x = x_sample,
    bwMtd = "hpi1Norm",
    bwMin = 1e-10,
    bwMax = 1e10,
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L,
    bwFallback = NULL,
    normPeakFrac = 0.15,
    normExtraFrac = 0.25
  )

  bw_pkg_std <- stimgate:::.bwCalcOne(
    x = x_sample,
    bwMtd = "hpi1Norm",
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L,
    normPeakFrac = 0.15,
    normExtraFrac = 0.25,
    adaptive = FALSE
  )

  expect_equal(as.numeric(bw_direct_std), as.numeric(bw_pkg_std))

  # Check adaptive normalized bandwidth via .bwDiagBwOne vs .bwCalcOne
  set.seed(42)
  bw_diag_adapt <- .bwDiagBwOne(
    x = x_sample,
    bw_mtd = "hpi1Norm",
    bw_adaptive = TRUE,
    bw_adj = 1,
    bw_ncell_min = 10L,
    bw_ncell_max = 1000L,
    norm_peak_frac = 0.15,
    norm_extra_frac = 0.25
  )

  set.seed(42)
  bw_pkg_adapt <- stimgate:::.bwCalcOne(
    x = x_sample,
    bwMtd = "hpi1Norm",
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L,
    normPeakFrac = 0.15,
    normExtraFrac = 0.25,
    adaptive = TRUE
  )

  expect_true(is.list(bw_diag_adapt))
  expect_true(is.list(bw_pkg_adapt))
  expect_equal(bw_diag_adapt$bw, bw_pkg_adapt$bw)
  expect_equal(bw_diag_adapt$bwCore, bw_pkg_adapt$bwCore)
  expect_equal(bw_diag_adapt$bwExtra, bw_pkg_adapt$bwExtra)
})
