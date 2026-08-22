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

test_that(".simBandwidthBsFreq fixed-seed parity checks match simcyto for gamma and gaussian scenarios", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_bw_io, local = env)
  source(script_bw_plot, local = env)

  run_case <- function(
      seed,
      transformation,
      mean_pos,
      bw,
      bias_uns,
      expected_abs_err) {
    n_sample <- 2L
    n_condition <- 2L
    n_cell_stim <- 240L
    n_cell_uns <- 120L
    prob_response <- 0.05
    background_relative_to_response <- 0.2
    prob_response_uns <- prob_response * background_relative_to_response

    captured_sim <- NULL
    orig_simcyto_experiment <- simcyto::simCytExperiment
    set.seed(seed)
    res <- testthat::with_mocked_bindings(
      simCytExperiment = function(...) {
        out <- orig_simcyto_experiment(...)
        captured_sim <<- out
        out
      },
      .package = "simcyto",
      env$.simBandwidthBsFreq(
        nSample = n_sample,
        nMarker = 1L,
        nCondition = n_condition,
        nCluster = 2L,
        nIter = 1L,
        biasUns = bias_uns,
        bw = bw,
        bwMin = "none",
        bwMax = "none",
        bwFallback = bw,
        nCellStim = n_cell_stim,
        probResponse = prob_response,
        probExact = TRUE,
        meanPos = mean_pos,
        transformation = transformation,
        samplePerturbationSd = 0.2,
        conditionPerturbationSd = 0.3,
        clusterPerturbationSd = 0.1,
        backgroundRelativeToResponse = background_relative_to_response,
        ncellUnsRelativeToStim = 0.5,
        covEvMin = 1.5,
        covEvMax = 1.5,
        tolClust = NULL,
        locEnforceShapeThreshold = FALSE,
        calcCytPosGates = FALSE
      )
    )

    truth_from_helper <- res |>
      dplyr::distinct(
        .data$sample,
        .data$ind,
        .data$propStimTruth,
        .data$propUnsTruth,
        .data$propRespTruth,
        .data$nCellStim,
        .data$nCellUns
      ) |>
      dplyr::arrange(.data$sample, .data$ind)

    set.seed(seed)
    sim <- simcyto::simCytExperiment(
      nSample = n_sample,
      nMarker = 1L,
      nCondition = n_condition,
      nCluster = 2L,
      nCellByCondition = c(n_cell_uns, n_cell_stim),
      transformationFunc = env$.simMiscGetTrans(transformation),
      mixtureType = "gaussianOnly",
      meanExprMat = matrix(c(0, mean_pos), byrow = TRUE, ncol = 1),
      clusterLabelVec = c("gn", "gp"),
      probVecUns = c(1 - prob_response_uns, prob_response_uns),
      probExact = TRUE,
      probResponseVecByStimCondition = list(c(-prob_response, prob_response)),
      conditionPerturbationSd = 0.3,
      clusterPerturbationSd = 0.1,
      samplePerturbationSd = 0.2,
      covEvMin = 1.5,
      covEvMax = 1.5
    )

    truth_from_simcyto <- purrr::map_df(seq_len(n_sample), function(sample_curr) {
      ind_uns <- (sample_curr - 1L) * n_condition + 1L
      ind_stim <- ind_uns + 1L
      labels_uns <- sim$labelsList[[ind_uns]]
      labels_stim <- sim$labelsList[[ind_stim]]
      prop_uns_truth <- sum(labels_uns %in% "gp") / length(labels_uns)
      prop_stim_truth <- sum(labels_stim %in% "gp") / length(labels_stim)

      tibble::tibble(
        sample = as.character(sample_curr),
        ind = as.character(ind_stim),
        propStimTruth = prop_stim_truth,
        propUnsTruth = prop_uns_truth,
        propRespTruth = prop_stim_truth - prop_uns_truth,
        nCellStim = n_cell_stim,
        nCellUns = n_cell_uns
      )
    }) |>
      dplyr::arrange(.data$sample, .data$ind)

    expect_equal(truth_from_helper, truth_from_simcyto, tolerance = 1e-12)

    expect_type(captured_sim, "list")
    expect_true(all(c("flowFrameList", "labelsList") %in% names(captured_sim)))

    expr_means_helper <- vapply(
      captured_sim$flowFrameList,
      function(ff) mean(flowCore::exprs(ff)[, 1]),
      numeric(1)
    )
    expr_sds_helper <- vapply(
      captured_sim$flowFrameList,
      function(ff) stats::sd(flowCore::exprs(ff)[, 1]),
      numeric(1)
    )
    expr_means_direct <- vapply(
      sim$flowFrameList,
      function(ff) mean(flowCore::exprs(ff)[, 1]),
      numeric(1)
    )
    expr_sds_direct <- vapply(
      sim$flowFrameList,
      function(ff) stats::sd(flowCore::exprs(ff)[, 1]),
      numeric(1)
    )
    expect_equal(unname(expr_means_helper), unname(expr_means_direct), tolerance = 1e-12)
    expect_equal(unname(expr_sds_helper), unname(expr_sds_direct), tolerance = 1e-12)
    expect_true(expr_means_helper[[2]] > expr_means_helper[[1]])
    expect_true(expr_means_helper[[4]] > expr_means_helper[[3]])

    abs_err <- res |>
      dplyr::filter(.data$method %in% c("loc_condition", "loc_sample")) |>
      dplyr::arrange(.data$sample, .data$ind, .data$method) |>
      dplyr::transmute(abs_err = abs(.data$propRespEst - .data$propRespTruth)) |>
      dplyr::pull(.data$abs_err)
    expect_equal(abs_err, expected_abs_err, tolerance = 1e-8)
  }

  run_case(
    seed = 2026L,
    transformation = "gamma",
    mean_pos = 4,
    bw = 0.02,
    bias_uns = 0.0025,
    expected_abs_err = c(0.05, 0.05, 0.0041666667, 0.0041666667)
  )

  run_case(
    seed = 2028L,
    transformation = "gaussian",
    mean_pos = 8,
    bw = 0.25,
    bias_uns = 0.05,
    expected_abs_err = c(0.0041666667, 0.0041666667, 0, 0)
  )
})
