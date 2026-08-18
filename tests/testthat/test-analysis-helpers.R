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
