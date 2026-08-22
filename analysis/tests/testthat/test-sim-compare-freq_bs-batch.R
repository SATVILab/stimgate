root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

test_that("analysis/8-sim-compare-freq_bs-batch.qmd does not source functionsForBenchmarking-Cyt.R", {
  qmd_path <- file.path(root_dir, "analysis", "8-sim-compare-freq_bs-batch.qmd")
  expect_true(file.exists(qmd_path))

  lines <- readLines(qmd_path, warn = FALSE)
  expect_false(
    any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
    info = "analysis/8-sim-compare-freq_bs-batch.qmd should not source functionsForBenchmarking-Cyt.R"
  )
})

test_that(".simCompareFreqBs forwards stimMeanShift and stimSdMultiplier to simcyto::simCytExperiment", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  orig_simcyto_experiment <- simcyto::simCytExperiment
  captured_shift <- NULL
  captured_sd_mult <- NULL

  testthat::with_mocked_bindings(
    simCytExperiment = function(..., stimMeanShift = 0, stimSdMultiplier = 1) {
      captured_shift <<- stimMeanShift
      captured_sd_mult <<- stimSdMultiplier
      orig_simcyto_experiment(..., stimMeanShift = stimMeanShift, stimSdMultiplier = stimSdMultiplier)
    },
    .package = "simcyto",
    {
      set.seed(42)
      res <- env$.simCompareFreqBs(
        nSample = 1L,
        nMarker = 1L,
        nCondition = 2L,
        nCluster = 2L,
        nIter = 1L,
        biasUns = 0,
        bw = 0.1,
        bwMtd = "hpi1",
        nCellStim = 200L,
        probResponse = 0.1,
        meanPos = 5,
        transformation = "gaussian",
        samplePerturbationSd = 0,
        conditionPerturbationSd = 0,
        clusterPerturbationSd = 0,
        backgroundRelativeToResponse = 0.1,
        ncellUnsRelativeToStim = 1,
        tailgateAutoTol = TRUE,
        stimMeanShift = 0.05,
        stimSdMultiplier = 1.05
      )

      expect_equal(captured_shift, 0.05)
      expect_equal(captured_sd_mult, 1.05)
      expect_s3_class(res, "data.frame")
      expect_true("stimMeanShift" %in% names(res))
      expect_true("stimSdMultiplier" %in% names(res))
      expect_equal(res$stimMeanShift[[1]], 0.05)
      expect_equal(res$stimSdMultiplier[[1]], 1.05)
    }
  )
})

test_that(".simCompareFreqBs with zero mismatch reproduces clean baseline", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  set.seed(123)
  res_clean <- env$.simCompareFreqBs(
    nSample = 2L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0,
    bw = 0.1,
    bwMtd = "hpi1",
    nCellStim = 300L,
    probResponse = 0.05,
    meanPos = 5,
    transformation = "gaussian",
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.1,
    ncellUnsRelativeToStim = 1,
    tailgateAutoTol = TRUE
  )

  set.seed(123)
  res_zero_mismatch <- env$.simCompareFreqBs(
    nSample = 2L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0,
    bw = 0.1,
    bwMtd = "hpi1",
    nCellStim = 300L,
    probResponse = 0.05,
    meanPos = 5,
    transformation = "gaussian",
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.1,
    ncellUnsRelativeToStim = 1,
    tailgateAutoTol = TRUE,
    stimMeanShift = 0,
    stimSdMultiplier = 1
  )

  # Disregarding the added stimMeanShift and stimSdMultiplier columns
  common_cols <- intersect(names(res_clean), names(res_zero_mismatch))
  expect_equal(res_clean[common_cols], res_zero_mismatch[common_cols])
})

test_that(".simCompareSummariseFreqBs correctly handles mismatch scenarios and summary metrics", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  grid <- data.frame(
    transformation = c("gaussian", "gaussian"),
    mean_pos = c(5, 5),
    prob_response = c(0.1, 0.1),
    n_cell = c(200, 200),
    bias_uns = c(0.15, 0.15),
    bw = c(0.1, 0.1),
    sample_perturbation_sd = c(0, 0),
    condition_perturbation_sd = c(0, 0),
    cluster_perturbation_sd = c(0, 0),
    background_relative_to_response = c(0.1, 0.1),
    n_cell_uns_relative_to_stim = c(1, 1),
    stim_mean_shift = c(0, 0.05),
    stim_sd_multiplier = c(1, 1),
    mismatch_type = c("mean_shift", "mean_shift"),
    mismatch_val = c(0, 0.05),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  raw_res <- env$.simCompareFreqBsGrid(
    sim_grid = grid,
    nSample = 1,
    nIter = 1,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    probExact = TRUE,
    tailgateAutoTol = TRUE
  )

  expect_s3_class(raw_res, "data.frame")
  expect_true(nrow(raw_res) > 0)

  summ <- env$.simCompareSummariseFreqBs(raw_res)
  expect_s3_class(summ, "data.frame")
  expect_true(all(c("med_abs_rel_error", "q90_abs_rel_error", "q95_abs_rel_error", "max_abs_rel_error") %in% names(summ)))
  expect_true("mismatch_val" %in% names(summ))
  expect_true(nrow(summ) > 0)
})
