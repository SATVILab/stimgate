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

test_that("changing biasUns modifies StimGate while leaving F-beta and tailgate inputs and estimates untouched (#307)", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  set.seed(42)
  trans <- simcyto::simCytTransformGaussian()
  exp_data <- simcyto::simCytExperiment(
    nSample = 2,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    nCellByCondition = c(500, 500),
    transformationFunc = trans,
    mixtureType = "gaussianOnly",
    meanExprMat = matrix(c(0, 5), ncol = 1),
    clusterLabelVec = c("gn", "gp"),
    probVecUns = c(0.99, 0.01),
    probExact = TRUE,
    probResponseVecByStimCondition = list(c(-0.05, 0.05)),
    covEvMin = 1.5,
    covEvMax = 1.5
  )

  # 1. Direct alternative rows test: biasUns must NOT affect F-beta or Tailgate
  alt_bias0 <- env$.simCompareAlternativeRows(
    flowFrameList = exp_data$flowFrameList,
    labelsList = exp_data$labelsList,
    nSample = 2,
    nCondition = 2,
    chnl = "F1",
    biasUns = 0,
    tailgateX = "stim"
  )

  alt_bias1 <- env$.simCompareAlternativeRows(
    flowFrameList = exp_data$flowFrameList,
    labelsList = exp_data$labelsList,
    nSample = 2,
    nCondition = 2,
    chnl = "F1",
    biasUns = 0.5,
    tailgateX = "stim"
  )

  # F-beta and Tailgate thresholds, proportions, and estimates MUST be identical
  expect_equal(alt_bias0$threshold, alt_bias1$threshold)
  expect_equal(alt_bias0$propStim, alt_bias1$propStim)
  expect_equal(alt_bias0$propUns, alt_bias1$propUns)
  expect_equal(alt_bias0$propRespEst, alt_bias1$propRespEst)

  # 2. Full .simCompareFreqBs test: biasUns changes StimGate, but not competitors
  set.seed(123)
  res_bias0 <- env$.simCompareFreqBs(
    nSample = 2,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    nIter = 1,
    biasUns = 0,
    bw = 0.1,
    bwMtd = "hpi1",
    nCellStim = 500,
    probResponse = 0.05,
    meanPos = 5,
    transformation = "gaussian",
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.1,
    ncellUnsRelativeToStim = 1,
    tailgateX = "stim"
  )

  set.seed(123)
  res_bias_pos <- env$.simCompareFreqBs(
    nSample = 2,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    nIter = 1,
    biasUns = 0.2,
    bw = 0.1,
    bwMtd = "hpi1",
    nCellStim = 500,
    probResponse = 0.05,
    meanPos = 5,
    transformation = "gaussian",
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.1,
    ncellUnsRelativeToStim = 1,
    tailgateX = "stim"
  )

  # StimGate threshold shifts with biasUns
  stim0 <- res_bias0[res_bias0$method == "stimgate", ]
  stim_pos <- res_bias_pos[res_bias_pos$method == "stimgate", ]
  expect_true(all(stim_pos$threshold >= stim0$threshold))

  # Competitors remain identical
  comp0 <- res_bias0[res_bias0$method %in% c("fbeta", "tailgate"), ]
  comp_pos <- res_bias_pos[res_bias_pos$method %in% c("fbeta", "tailgate"), ]
  expect_equal(comp0$threshold, comp_pos$threshold)
  expect_equal(comp0$propRespEst, comp_pos$propRespEst)
})

test_that("primary StimGate comparator scores full cluster-refined procedure (#308)", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  set.seed(42)
  res <- env$.simCompareFreqBs(
    nSample = 3,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    nIter = 1,
    biasUns = 0,
    bw = 0.1,
    bwMtd = "hpi1",
    nCellStim = 500,
    probResponse = 0.05,
    meanPos = 5,
    transformation = "gaussian",
    samplePerturbationSd = 0.5,
    conditionPerturbationSd = 0.5,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.1,
    ncellUnsRelativeToStim = 1,
    tolClust = 0.01,
    includeLocCondition = TRUE,
    tailgateX = "stim"
  )

  # Primary stimgate rows must exist and be named 'stimgate'
  stimgate_primary <- res[res$method == "stimgate", ]
  expect_equal(nrow(stimgate_primary), 3L)
  expect_true(all(stimgate_primary$approach == "stimgate"))
  expect_true(all(stimgate_primary$detailLevel == "cluster_final"))
  expect_true(all(is.finite(stimgate_primary$threshold)))
  expect_true(all(is.finite(stimgate_primary$propRespEst)))

  # Pre-clustering diagnostic rows stimgate_loc_sample are also preserved
  stimgate_loc_sample <- res[res$method == "stimgate_loc_sample", ]
  expect_equal(nrow(stimgate_loc_sample), 3L)
  expect_true(all(stimgate_loc_sample$detailLevel == "sample"))

  # Summaries default to keeping c("stimgate", "fbeta", "tailgate")
  summ <- env$.simCompareSummariseFreqBs(res)
  expect_setequal(unique(summ$method), c("stimgate", "fbeta", "tailgate"))
})

test_that("tailgate derives threshold from stimulated sample and applies to both conditions (#308)", {
  testthat::skip_if_not_installed("cytoUtils")
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  set.seed(42)
  trans <- simcyto::simCytTransformGaussian()
  exp_data <- simcyto::simCytExperiment(
    nSample = 1,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    nCellByCondition = c(500, 500),
    transformationFunc = trans,
    mixtureType = "gaussianOnly",
    meanExprMat = matrix(c(0, 5), ncol = 1),
    clusterLabelVec = c("gn", "gp"),
    probVecUns = c(0.99, 0.01),
    probExact = TRUE,
    probResponseVecByStimCondition = list(c(-0.05, 0.05)),
    covEvMin = 1.5,
    covEvMax = 1.5
  )

  alt_stim <- env$.simCompareAlternativeRows(
    flowFrameList = exp_data$flowFrameList,
    labelsList = exp_data$labelsList,
    nSample = 1,
    nCondition = 2,
    chnl = "F1",
    tailgateX = "stim"
  )

  tailgate_row <- alt_stim[alt_stim$method == "tailgate", ]
  expect_equal(nrow(tailgate_row), 1L)

  xStim <- as.numeric(flowCore::exprs(exp_data$flowFrameList[[2]])[, "F1"])
  xUns <- as.numeric(flowCore::exprs(exp_data$flowFrameList[[1]])[, "F1"])

  # Direct tailgate call on stimulated cells must match threshold
  direct_tg <- env$.simCompareTailgateThreshold(x = xStim, autoTol = TRUE)
  expect_equal(tailgate_row$threshold, direct_tg$threshold)

  # Background-subtracted propRespEst must equal propStim - propUns on raw data
  expected_prop_stim <- sum(xStim > direct_tg$threshold) / length(xStim)
  expected_prop_uns <- sum(xUns > direct_tg$threshold) / length(xUns)
  expect_equal(tailgate_row$propStim, expected_prop_stim)
  expect_equal(tailgate_row$propUns, expected_prop_uns)
  expect_equal(tailgate_row$propRespEst, expected_prop_stim - expected_prop_uns)
})

test_that("F-beta threshold metric is preserved as diagnostic without altering threshold (#308)", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  set.seed(42)
  trans <- simcyto::simCytTransformGaussian()
  exp_data <- simcyto::simCytExperiment(
    nSample = 1,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    nCellByCondition = c(300, 300),
    transformationFunc = trans,
    mixtureType = "gaussianOnly",
    meanExprMat = matrix(c(0, 5), ncol = 1),
    clusterLabelVec = c("gn", "gp"),
    probVecUns = c(0.99, 0.01),
    probExact = TRUE,
    probResponseVecByStimCondition = list(c(-0.05, 0.05)),
    covEvMin = 1.5,
    covEvMax = 1.5
  )

  alt <- env$.simCompareAlternativeRows(
    flowFrameList = exp_data$flowFrameList,
    labelsList = exp_data$labelsList,
    nSample = 1,
    nCondition = 2,
    chnl = "F1"
  )

  fbeta_row <- alt[alt$method == "fbeta", ]
  expect_equal(nrow(fbeta_row), 1L)
  expect_true(is.numeric(fbeta_row$thresholdMetric))
  expect_true(is.finite(fbeta_row$threshold))
})
