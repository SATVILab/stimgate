root_dir <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = FALSE)
if (!file.exists(file.path(root_dir, "scripts", "r", "sim-bandwidth.R"))) {
  root_dir <- normalizePath(
    file.path(testthat::test_path(), "../../.."),
    mustWork = FALSE
  )
}

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

test_that("analysis wrapper functions do not expose removed calcSinglePosGates argument", {
  skip_if_not(file.exists(script_bw), "sim-bandwidth.R not found")
  skip_if_not(file.exists(script_comp), "sim-compare-freq_bs.R not found")

  env <- new.env(parent = baseenv())
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  expect_false("calcSinglePosGates" %in% names(formals(env$.simBandwidthBsFreq)))
  expect_false("calcSinglePosGates" %in% names(formals(env$.simCompareStimgateRows)))
  expect_false("calcSinglePosGates" %in% names(formals(env$.simCompareFreqBs)))
})

test_that(".simBandwidthBsFreq forwards to gateStim without unknown-argument error", {
  skip_if_not_installed("stimgate")
  skip_if_not(file.exists(script_bw), "sim-bandwidth.R not found")

  env <- new.env(parent = baseenv())
  source(script_misc, local = env)
  source(script_bw, local = env)

  res <- env$.simBandwidthBsFreq(
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
