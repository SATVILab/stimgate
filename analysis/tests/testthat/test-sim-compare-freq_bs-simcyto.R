root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

test_that("comparison simulation helpers source cleanly without legacy functionsForBenchmarking-Cyt.R", {
  for (f in c(script_misc, script_bw, script_comp)) {
    if (!file.exists(f)) stop("Expected analysis helper not found: ", f)
  }

  env <- new.env(parent = getNamespace("stimgate"))
  expect_no_error(source(script_misc, local = env))
  expect_no_error(source(script_bw, local = env))
  expect_no_error(source(script_comp, local = env))

  # Ensure no local unexported simCytExperiment is introduced into the helper environment
  expect_false(exists("simCytExperiment", envir = env, inherits = FALSE))
})

test_that("analysis/7-sim-compare-freq_bs.qmd does not source functionsForBenchmarking-Cyt.R", {
  qmd_path <- file.path(root_dir, "analysis", "7-sim-compare-freq_bs.qmd")
  expect_true(file.exists(qmd_path))

  lines <- readLines(qmd_path, warn = FALSE)
  expect_false(
    any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
    info = "analysis/7-sim-compare-freq_bs.qmd should not source functionsForBenchmarking-Cyt.R"
  )
})

test_that(".simCompareFreqBs calls simcyto::simCytExperiment and preserves comparison structure", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

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
      res_gauss <- env$.simCompareFreqBs(
        nSample = 2L,
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
        tailgateAutoTol = TRUE
      )

      expect_true(called_simcyto)
      expect_s3_class(res_gauss, "data.frame")
      expect_true(nrow(res_gauss) > 0)

      # Check approaches present
      expect_true(all(c("stimgate", "fbeta", "tailgate") %in% res_gauss$approach))

      # Check StimGate results
      stimgate_res <- res_gauss[res_gauss$approach == "stimgate", , drop = FALSE]
      expect_true(nrow(stimgate_res) > 0)
      expect_false(any(stimgate_res$method == "stimgate_error"))
      expect_true(all(is.na(stimgate_res$error)))
      expect_true(all(is.finite(stimgate_res$threshold)))
      expect_true(all(is.finite(stimgate_res$propRespEst)))

      # Check truth proportions
      expect_true(all(is.finite(res_gauss$propStimTruth)))
      expect_true(all(is.finite(res_gauss$propUnsTruth)))
      expect_true(all(is.finite(res_gauss$propRespTruth)))
      expect_equal(
        res_gauss$propRespTruth,
        res_gauss$propStimTruth - res_gauss$propUnsTruth,
        tolerance = 1e-10
      )
    }
  )
})

test_that(".simCompareFreqBs works with gamma and skew transformations and rare responses", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  scenarios <- list(
    list(trans = "gamma", meanPos = 4, probResp = 0.002),
    list(trans = "skew", meanPos = 6, probResp = 0.05)
  )

  for (sc in scenarios) {
    set.seed(99)
    res <- env$.simCompareFreqBs(
      nSample = 1L,
      nMarker = 1L,
      nCondition = 2L,
      nCluster = 2L,
      nIter = 1L,
      biasUns = 0,
      bw = 0.1,
      bwMtd = "hpi1",
      nCellStim = 300L,
      probResponse = sc$probResp,
      meanPos = sc$meanPos,
      transformation = sc$trans,
      samplePerturbationSd = 0,
      conditionPerturbationSd = 0,
      clusterPerturbationSd = 0,
      backgroundRelativeToResponse = 0.1,
      ncellUnsRelativeToStim = 1,
      tailgateAutoTol = TRUE
    )

    expect_s3_class(res, "data.frame")
    expect_true(nrow(res) > 0)

    stimgate_res <- res[res$approach == "stimgate", , drop = FALSE]
    expect_true(nrow(stimgate_res) > 0)
    expect_true(all(is.finite(stimgate_res$threshold)))
    expect_true(all(is.finite(stimgate_res$propRespEst)))
  }
})
