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
  withr::local_preserve_seed()
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


test_that("analysis 6 uses shared transactional runners and full-grid reruns", {
  content <- paste(readLines(file.path(
    root_dir, "analysis", "6-sim-bw-freq_bs-adaptive.qmd"
  )), collapse = "\n")
  for (contract in c(
    'analysis_semantics_version <- "adaptive-bw-freq-v3"',
    "sim_grid_full <- sim_grid",
    "sim_grid_spec = analysis_grid_spec",
    "scenario_settings = scenario_settings",
    "nSample = n_sample_sim", "nIter = n_iter_sim",
    ".simBandwidthRunRow(", ".simBandwidthRunGrid(",
    ".simBandwidthFinishChunk(", "retry_errors = TRUE",
    "validate_fn = .simBandwidthFreqBsAdaptiveValidate",
    "required_params = analysis_required_params",
    "run_ctx <- .analysis_results_context(",
    "Skipping plots during a multi-chunk simulation render."
  )) {
    expect_true(grepl(contract, content, fixed = TRUE), info = contract)
  }
  expect_false(grepl("promote_analysis6_if_ready", content, fixed = TRUE))
  expect_false(grepl("purrr::flatten", content, fixed = TRUE))
  expect_false(grepl("dens_tbl", content, fixed = TRUE))
  expect_false(grepl("rug_tbl", content, fixed = TRUE))
  expect_false(grepl(".write_rds_atomic", content, fixed = TRUE))
  expect_equal(length(gregexpr("#| eval: false", content,
                               fixed = TRUE)[[1]]), 1L)
  expect_lt(regexpr("sim_grid_full <-", content, fixed = TRUE)[[1]],
            regexpr("if (analysis_quick", content, fixed = TRUE)[[1]])
})

.load_adaptive_run_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "sim-misc.R", "sim-bandwidth.R",
                 "sim-bandwidth-analysis-io.R",
                 "sim-bandwidth-analysis-run.R")) {
    source(file.path(root_dir, "scripts", "r", file), local = env)
  }
  env
}

.adaptive_grid_row <- function() {
  tibble::tibble(
    sim_id = 7L, sim_seed = 124L, bias_uns = 0.05,
    bw_core = 0.15, bw_extra = 0.25, bw_crossover = NA_real_,
    bw_transition_width = 0, bw_fallback = 0.5, n_cell = 1000,
    prob_response = 0.02, mean_pos = 6, transformation = "skew",
    sample_perturbation_sd = 0, condition_perturbation_sd = 0,
    cluster_perturbation_sd = 0, background_relative_to_response = 0.2,
    n_cell_uns_relative_to_stim = 1
  )
}

test_that("adaptive scenario forwards settings and reruns identical random draws", {
  withr::local_preserve_seed()
  env <- .load_adaptive_run_env()
  captured <- NULL
  env$.simBandwidthBsFreq <- function(...) {
    captured <<- list(...)
    tibble::tibble(iter = 1L, ind = 2L, sample = "1",
                   method = "loc_sample", threshold = stats::runif(1),
                   propRespTruth = 0.02, propRespEst = 0.025)
  }
  row <- .adaptive_grid_row()
  settings <- list(nSample = 3L, nIter = 2L, bw = NULL, bwAdaptive = TRUE)
  withr::local_rng_version("4.4.0")
  set.seed(99L, kind = "L'Ecuyer-CMRG")
  before <- .Random.seed
  first <- env$.simBandwidthRunRow(
    row, env$.simBandwidthFreqBsAdaptiveScenario, settings
  )
  expect_identical(.Random.seed, before)
  expect_identical(captured$nSample, 3L)
  expect_identical(captured$nIter, 2L)
  expect_null(captured$bwAdaptiveCrossover)
  expect_identical(captured$bwAdaptiveCore, row$bw_core)
  expect_identical(captured$bwAdaptiveExtra, row$bw_extra)
  set.seed(55L, kind = "Mersenne-Twister")
  second <- env$.simBandwidthRunRow(
    row, env$.simBandwidthFreqBsAdaptiveScenario, settings
  )
  expect_identical(first, second)
  row$bw_crossover <- 5.5
  env$.simBandwidthRunRow(row, env$.simBandwidthFreqBsAdaptiveScenario,
                         settings)
  expect_identical(captured$bwAdaptiveCrossover, 5.5)
  error <- env$.simBandwidthErrorRow(row, "failed")
  expect_type(dplyr::bind_rows(first, error)$ind, "integer")
})

test_that("adaptive collation scores final sample frequencies and keeps summaries", {
  env <- .load_adaptive_run_env()
  tbl <- tibble::tibble(
    sim_id = c(7L, 7L, 7L), sim_seed = 124L,
    iter = 1L, ind = c(2L, 4L, 2L), sample = c("1", "2", "1"),
    method = c("loc_sample", "loc_sample", "loc_condition"),
    threshold = c(1, 3, 99), propRespTruth = 0.02,
    propRespEst = c(0.01, 0.03, 1), propBsEst = 0.9
  )
  expect_identical(env$.simBandwidthFreqBsAdaptiveValidate(tbl), character())
  result <- env$.simBandwidthFreqBsAdaptiveCollate(tbl,
                                                 c("sim_id", "sim_seed"))
  expect_named(result, c("bw_tbl_results_raw", "bw_tbl_results_summary",
                         "summary_tbl"))
  expect_equal(nrow(result$bw_tbl_results_raw), 2L)
  expect_equal(result$summary_tbl$threshold_median, 2)
  expect_equal(result$summary_tbl$threshold_iqr_length, 1)
  expect_equal(result$bw_tbl_results_summary$propRespEst_median_diff, 0,
               tolerance = 1e-12)
  duplicate <- dplyr::bind_rows(tbl, tbl[1, ])
  expect_match(env$.simBandwidthFreqBsAdaptiveValidate(duplicate), "Duplicate")
  missing <- tbl
  missing$threshold <- NA_real_
  expect_match(env$.simBandwidthFreqBsAdaptiveValidate(missing),
               "no finite final")
  expect_match(env$.simBandwidthFreqBsAdaptiveValidate(tbl["sim_id"]),
               "Missing final")
})

test_that("adaptive failed rows retry and promoted reads enforce grid settings", {
  env <- .load_adaptive_run_env()
  project <- withr::local_tempdir()
  withr::local_dir(project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")
  row <- .adaptive_grid_row()
  required <- list(
    analysis_semantics_version = "adaptive-bw-freq-v3",
    sim_grid_spec = row[, setdiff(names(row), "sim_seed")],
    scenario_settings = list(nSample = 5L, nIter = 5L)
  )
  ctx <- env$.analysis_run_context(
    c("sim", "bw", "freq_bs", "adaptive"), run_id = "adaptive-retry",
    path_root = project, params = required
  )
  fail <- TRUE
  env$.simBandwidthBsFreq <- function(...) {
    if (fail) stop("temporary failure")
    tibble::tibble(iter = 1L, ind = 2L, sample = "1",
                   method = "loc_sample", threshold = stats::runif(1),
                   propRespTruth = 0.02, propRespEst = 0.025)
  }
  run <- function() env$.simBandwidthRunRowResumable(
    row, env$.simBandwidthFreqBsAdaptiveScenario,
    required$scenario_settings, ctx, total_sims = 1L
  )
  expect_match(run()$error_message, "temporary failure")
  expect_error(env$.simBandwidthFinishChunk(
    ctx, row, row,
    validate_fn = env$.simBandwidthFreqBsAdaptiveValidate,
    collate_fn = function(tbl) {
      env$.simBandwidthFreqBsAdaptiveCollate(tbl, names(row))
    }
  ), "simulation errors")
  fail <- FALSE
  retried <- run()
  expect_true(all(is.na(retried$error_message)))
  expect_false(file.exists(file.path(ctx$chunk_jobs_dir, "error-7")))
  expect_identical(retried, env$.simBandwidthRunRow(
    row, env$.simBandwidthFreqBsAdaptiveScenario, required$scenario_settings
  ))
  expect_true(env$.simBandwidthFinishChunk(
    ctx, row, row,
    validate_fn = env$.simBandwidthFreqBsAdaptiveValidate,
    collate_fn = function(tbl) {
      env$.simBandwidthFreqBsAdaptiveCollate(tbl, names(row))
    }
  ))
  read_ctx <- env$.analysis_results_context(ctx$analysis_key,
                                           path_root = project)
  summary_path <- env$.analysis_current_file(
    read_ctx, c("collated", "summary_tbl.rds"), required_params = required
  )
  expect_equal(readRDS(summary_path)$threshold_median, retried$threshold)
  expect_true(read_ctx$read_only)
  changed <- required
  changed$sim_grid_spec$bw_core <- 0.5
  expect_error(env$.analysis_current_file(
    read_ctx, c("collated", "summary_tbl.rds"), required_params = changed
  ), "sim_grid_spec")
  changed <- required
  changed$scenario_settings$nSample <- 6L
  expect_error(env$.analysis_current_file(
    read_ctx, c("collated", "summary_tbl.rds"), required_params = changed
  ), "scenario_settings")
})
