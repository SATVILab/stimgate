root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_bw_io <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")

test_that("adaptive bandwidth simulation helpers source cleanly without legacy functionsForBenchmarking-Cyt.R", {
  for (f in c(script_misc, script_bw, script_bw_io, script_bw_plot)) {
    if (!file.exists(f)) stop("Expected analysis helper not found: ", f)
  }

  env <- new.env(parent = getNamespace("stimgate"))
  expect_no_error(source(script_misc, local = env))
  expect_no_error(source(script_bw, local = env))
  expect_no_error(source(script_bw_io, local = env))
  expect_no_error(source(script_bw_plot, local = env))

  expect_false(exists("simCytExperiment", envir = env, inherits = FALSE))
})

test_that("analysis/6-sim-bw-freq_bs-adaptive.qmd does not source functionsForBenchmarking-Cyt.R", {
  qmd_path <- file.path(root_dir, "analysis", "6-sim-bw-freq_bs-adaptive.qmd")
  expect_true(file.exists(qmd_path))

  lines <- readLines(qmd_path, warn = FALSE)
  expect_false(
    any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
    info = "analysis/6-sim-bw-freq_bs-adaptive.qmd should not source functionsForBenchmarking-Cyt.R"
  )
})

test_that(".simBandwidthBsFreq adaptive fixed-seed parity checks match simcyto for gamma and skew scenarios", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_bw_io, local = env)
  source(script_bw_plot, local = env)

  run_case <- function(
    seed,
    transformation,
    mean_pos,
    bias_uns,
    bw_core,
    bw_extra,
    bw_fallback,
    bw_crossover,
    bw_transition_width,
    expected_abs_err
  ) {
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
        bw = NULL,
        bwAdaptive = TRUE,
        bwAdaptiveCore = bw_core,
        bwAdaptiveExtra = bw_extra,
        bwAdaptiveCrossover = bw_crossover,
        bwAdaptiveTransitionWidth = bw_transition_width,
        bwFallback = bw_fallback,
        bwMin = "none",
        bwMax = "none",
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

    expect_equal(unique(res$bwAdaptive), TRUE)
    expect_equal(unique(res$bwAdaptiveCore), bw_core)
    expect_equal(unique(res$bwAdaptiveExtra), bw_extra)
    expect_equal(unique(res$bwAdaptiveTransitionWidth), bw_transition_width)
    expect_equal(unique(res$nCellStim), n_cell_stim)
    expect_equal(unique(res$nCellUns), n_cell_uns)

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

    expr_helper <- lapply(captured_sim$flowFrameList, function(ff) flowCore::exprs(ff)[, 1])
    expr_direct <- lapply(sim$flowFrameList, function(ff) flowCore::exprs(ff)[, 1])
    expect_equal(expr_helper, expr_direct, tolerance = 1e-12)

    abs_err <- res |>
      dplyr::filter(.data$method %in% c("loc_condition", "loc_sample")) |>
      dplyr::arrange(.data$sample, .data$ind, .data$method) |>
      dplyr::transmute(abs_err = abs(.data$propRespEst - .data$propRespTruth)) |>
      dplyr::pull(.data$abs_err)
    expect_equal(abs_err, expected_abs_err, tolerance = 1e-8)
  }

  run_case(
    seed = 2926L,
    transformation = "gamma",
    mean_pos = 4,
    bias_uns = 0.0025,
    bw_core = 0.02,
    bw_extra = 0.03,
    bw_fallback = 0.01,
    bw_crossover = NA_real_,
    bw_transition_width = 0,
    expected_abs_err = c(0, 0, 0.0416666667, 0.0416666667)
  )

  run_case(
    seed = 2928L,
    transformation = "skew",
    mean_pos = 6,
    bias_uns = 0.05,
    bw_core = 0.25,
    bw_extra = 0.5,
    bw_fallback = 0.5,
    bw_crossover = 5.5,
    bw_transition_width = 0.25,
    expected_abs_err = c(0.0041666667, 0.0041666667, 0.0083333333, 0.0083333333)
  )
})
