root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_bw_io <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")

test_that("Analysis 2a executes its original bandwidth and fixed-bias design", {
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "sim-misc.R", "sim-bandwidth.R",
    "sim-bandwidth-analysis-io.R", "sim-bandwidth-analysis-run.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  lines <- readLines(file.path(
    root_dir, "analysis", "2a-sim-bw-freq_bs-global.qmd"
  ))
  chunk <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1L]
    lines[(start + 1L):(end - 1L)]
  }
  env$analysis_quick <- FALSE
  env$analysis_dev <- FALSE
  env$analysis_semantics_version <- "test"
  env$simulation_seed <- 12345L
  env$sim_grid_shuffle_seed <- 8L
  env$sim_grid_chunk_index <- 1L
  env$sim_grid_n_chunks <- 1L
  eval(parse(text = chunk("actual-settings")), envir = env)
  invisible(utils::capture.output(
    eval(parse(text = chunk("bw-manual-grid")), envir = env)
  ))
  eval(parse(text = chunk("bw-manual-settings")), envir = env)

  expect_equal(env$scenario_settings$nSample, 200)
  expect_null(env$scenario_settings$tolClust)
  expect_false(env$scenario_settings$locEnforceShapeThreshold)
  expect_false(env$scenario_settings$calcCytPosGates)
  expect_equal(sort(unique(env$sim_grid_all$n_cell)), c(1e3, 5e3, 2e4, 1e5))
  expected_pairs <- tidyr::expand_grid(
    n_cell = c(1e3, 5e3, 2e4, 1e5),
    prob_response = c(1 / 2e5, 1 / 5e4, 1 / 1e4, 1 / 2e3, 1 / 5e2, 1 / 5)
  ) |>
    dplyr::filter(!(.data$n_cell * .data$prob_response < 5 &
      .data$prob_response < 0.04)) |>
    dplyr::arrange(.data$n_cell, .data$prob_response)
  actual_pairs <- env$sim_grid_all |>
    dplyr::distinct(.data$n_cell, .data$prob_response) |>
    dplyr::arrange(.data$n_cell, .data$prob_response)
  expect_equal(actual_pairs, expected_pairs)
  expect_equal(
    sort(unique(env$sim_grid_all$condition_perturbation_sd)), c(0, 0.5)
  )
  expect_true(all(env$sim_grid_all$sample_perturbation_sd == 0))
  expect_true(all(env$sim_grid_all$cluster_perturbation_sd == 0))
  expect_true(all(env$sim_grid_all$background_relative_to_response == 0.2))
  expect_true(all(env$sim_grid_all$n_cell_uns_relative_to_stim == 1))
  for (transformation in c("gaussian", "skew", "gamma")) {
    rows <- env$sim_grid_all |>
      dplyr::filter(.data$transformation == .env$transformation)
    expect_equal(
      sort(unique(rows$mean_pos)),
      switch(transformation, gaussian = c(4.5, 8), skew = c(6, 8.5), gamma = c(4, 7))
    )
    expected_bw <- if (transformation == "gamma") {
      c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015,
        0.0175, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15)
    } else {
      c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1, 1.25, 1.5)
    }
    expect_equal(sort(unique(rows$bw)), expected_bw)
    bias_rules <- rows |>
      dplyr::distinct(.data$bias_uns_setting, .data$bias_uns) |>
      dplyr::arrange(.data$bias_uns)
    expect_identical(bias_rules$bias_uns_setting, c("none", "low", "high"))
    expect_equal(
      bias_rules$bias_uns,
      if (transformation == "gamma") c(0, 0.0025, 0.01) else c(0, 0.05, 0.25)
    )
  }
})

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

test_that("analysis/2a-sim-bw-freq_bs-global.qmd does not source functionsForBenchmarking-Cyt.R", {
  qmd_path <- file.path(root_dir, "analysis", "2a-sim-bw-freq_bs-global.qmd")
  expect_true(file.exists(qmd_path))

  lines <- readLines(qmd_path, warn = FALSE)
  expect_false(
    any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
    info = "analysis/2a-sim-bw-freq_bs-global.qmd should not source functionsForBenchmarking-Cyt.R"
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
