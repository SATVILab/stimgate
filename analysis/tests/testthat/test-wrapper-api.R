root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_cyt <- file.path(root_dir, "scripts", "r", "functionsForBenchmarking-Cyt.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_bw_io <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")
script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

.load_analysis_env <- function() {
  for (f in c(script_misc, script_cyt, script_bw, script_bw_io, script_bw_plot, script_comp)) {
    if (!file.exists(f)) stop("Expected analysis helper not found: ", f)
  }
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_cyt, local = env)
  source(script_bw, local = env)
  source(script_bw_io, local = env)
  source(script_bw_plot, local = env)
  source(script_comp, local = env)
  env
}

test_that("analysis wrapper functions do not expose removed calcSinglePosGates argument", {
  env <- .load_analysis_env()

  expect_false("calcSinglePosGates" %in% names(formals(env$.simBandwidthBsFreq)))
  expect_false("calcSinglePosGates" %in% names(formals(env$.simCompareStimgateRows)))
  expect_false("calcSinglePosGates" %in% names(formals(env$.simCompareFreqBs)))
})

test_that(".simBandwidthBsFreq exposes adaptive bandwidth settings directly", {
  env <- .load_analysis_env()

  adaptive_nm <- c(
    "bwAdaptive",
    "bwAdaptiveCore",
    "bwAdaptiveExtra",
    "bwAdaptiveCrossover",
    "bwAdaptiveTransitionWidth"
  )

  expect_true(all(adaptive_nm %in% names(formals(env$.simBandwidthBsFreq))))
})

test_that(".simBandwidthBsFreq forwards to gateStim without unknown-argument error", {
  env <- .load_analysis_env()

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

test_that(".simCompareFreqBs forwards to gateStim via .simCompareStimgateRows without unknown-argument error", {
  env <- .load_analysis_env()

  res <- env$.simCompareFreqBs(
    nSample = 1L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0,
    bw = 0.1,
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

  stimgate_res <- res[res[["approach"]] == "stimgate", , drop = FALSE]
  expect_true(nrow(stimgate_res) > 0)
  expect_false(any(stimgate_res[["method"]] == "stimgate_error"))
  expect_true(all(is.na(stimgate_res[["error"]])))
})
