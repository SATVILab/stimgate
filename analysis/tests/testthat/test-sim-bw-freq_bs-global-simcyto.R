root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_bw_io <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")

test_that("global bandwidth simulation helpers source cleanly without legacy functionsForBenchmarking-Cyt.R", {
  for (f in c(script_misc, script_bw, script_bw_io, script_bw_plot)) {
    if (!file.exists(f)) stop("Expected analysis helper not found: ", f)
  }

  env <- new.env(parent = getNamespace("stimgate"))
  expect_no_error(source(script_misc, local = env))
  expect_no_error(source(script_bw, local = env))
  expect_no_error(source(script_bw_io, local = env))
  expect_no_error(source(script_bw_plot, local = env))

  # Ensure no local unexported simCytExperiment is introduced into the helper environment
  expect_false(exists("simCytExperiment", envir = env, inherits = FALSE))
})

test_that("analysis/2-sim-bw-freq_bs-global.qmd does not source functionsForBenchmarking-Cyt.R", {
  qmd_path <- file.path(root_dir, "analysis", "2-sim-bw-freq_bs-global.qmd")
  expect_true(file.exists(qmd_path))

  lines <- readLines(qmd_path, warn = FALSE)
  expect_false(
    any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
    info = "analysis/2-sim-bw-freq_bs-global.qmd should not source functionsForBenchmarking-Cyt.R"
  )
})

test_that(".simBandwidthBsFreq calls simcyto::simCytExperiment and produces valid Gaussian bandwidth results", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_bw_io, local = env)
  source(script_bw_plot, local = env)

  orig_simcyto_experiment <- simcyto::simCytExperiment
  called_simcyto <- FALSE

  testthat::with_mocked_bindings(
    simCytExperiment = function(...) {
      called_simcyto <<- TRUE
      orig_simcyto_experiment(...)
    },
    .package = "simcyto",
    {
      set.seed(42)
      res_gauss <- env$.simBandwidthBsFreq(
        nSample = 2L,
        nMarker = 1L,
        nCondition = 2L,
        nCluster = 2L,
        nIter = 1L,
        biasUns = 0.05,
        bw = 0.1,
        bwMin = "none",
        bwMax = "none",
        bwFallback = 0.1,
        nCellStim = 200L,
        probResponse = 0.05,
        probExact = TRUE,
        meanPos = 8,
        transformation = "gaussian",
        samplePerturbationSd = 0,
        conditionPerturbationSd = 0,
        clusterPerturbationSd = 0,
        backgroundRelativeToResponse = 0.2,
        ncellUnsRelativeToStim = 1,
        covEvMin = 1.5,
        covEvMax = 1.5,
        tolClust = NULL,
        locEnforceShapeThreshold = FALSE,
        calcCytPosGates = FALSE
      )

      expect_true(called_simcyto)
      expect_s3_class(res_gauss, "tbl_df")
      expect_true(nrow(res_gauss) > 0)
      expect_equal(unique(res_gauss$nCellStim), 200L)
      expect_equal(unique(res_gauss$nCellUns), 200L)

      # Check truth proportions
      expect_true(all(is.finite(res_gauss$propStimTruth)))
      expect_true(all(is.finite(res_gauss$propUnsTruth)))
      expect_true(all(is.finite(res_gauss$propRespTruth)))
      expect_equal(
        res_gauss$propRespTruth,
        res_gauss$propStimTruth - res_gauss$propUnsTruth,
        tolerance = 1e-10
      )

      # Check methods and estimated proportions
      expect_true(all(c("propRespSmooth", "propRespPred", "loc_condition", "loc_sample") %in% res_gauss$method))
      expect_true(all(is.finite(res_gauss$propRespEst)))
    }
  )
})

test_that(".simBandwidthBsFreq works with gamma and skew transformations from simcyto", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_bw_io, local = env)
  source(script_bw_plot, local = env)

  for (tr in c("gamma", "skew")) {
    set.seed(123)
    res <- env$.simBandwidthBsFreq(
      nSample = 2L,
      nMarker = 1L,
      nCondition = 2L,
      nCluster = 2L,
      nIter = 1L,
      biasUns = 0.05,
      bw = if (identical(tr, "gamma")) 0.02 else 0.25,
      bwMin = "none",
      bwMax = "none",
      bwFallback = if (identical(tr, "gamma")) 0.02 else 0.25,
      nCellStim = 200L,
      probResponse = 0.05,
      probExact = TRUE,
      meanPos = if (identical(tr, "gamma")) 4 else 6,
      transformation = tr,
      samplePerturbationSd = 0,
      conditionPerturbationSd = 0,
      clusterPerturbationSd = 0,
      backgroundRelativeToResponse = 0.2,
      ncellUnsRelativeToStim = 1,
      covEvMin = 1.5,
      covEvMax = 1.5,
      tolClust = NULL,
      locEnforceShapeThreshold = FALSE,
      calcCytPosGates = FALSE
    )

    expect_s3_class(res, "tbl_df")
    expect_true(nrow(res) > 0)
    expect_equal(unique(res$nCellStim), 200L)
    expect_equal(unique(res$nCellUns), 200L)
    expect_true(all(is.finite(res$propRespEst)))
  }
})

test_that(".simBandwidthBsFreq correctly preserves perturbations and cell count ratios", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_bw_io, local = env)
  source(script_bw_plot, local = env)

  set.seed(456)
  res_ratio <- env$.simBandwidthBsFreq(
    nSample = 2L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0.05,
    bw = 0.25,
    bwMin = "none",
    bwMax = "none",
    bwFallback = 0.25,
    nCellStim = 300L,
    probResponse = 0.05,
    probExact = TRUE,
    meanPos = 8,
    transformation = "gaussian",
    samplePerturbationSd = 0.5,
    conditionPerturbationSd = 0.5,
    clusterPerturbationSd = 0.2,
    backgroundRelativeToResponse = 0.2,
    ncellUnsRelativeToStim = 0.5,
    covEvMin = 1.5,
    covEvMax = 1.5,
    tolClust = NULL,
    locEnforceShapeThreshold = FALSE,
    calcCytPosGates = FALSE
  )

  expect_s3_class(res_ratio, "tbl_df")
  expect_equal(unique(res_ratio$nCellStim), 300L)
  expect_equal(unique(res_ratio$nCellUns), 150L)
  expect_true(all(is.finite(res_ratio$propRespEst)))
})
