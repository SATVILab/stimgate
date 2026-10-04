.bias_uns_test_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "sim-misc.R", "sim-bandwidth.R",
    "sim-bandwidth-analysis-io.R", "sim-bandwidth-analysis-run.R"
  )) {
    source(file.path(root, "scripts", "r", fn), local = env)
  }
  env
}

test_that("Analysis 2b executes the agreed grid with shared biological seeds", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  lines <- readLines(file.path(root, "analysis", "2b-sim-bias_uns-freq_bs.qmd"))
  chunk <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1L]
    lines[(start + 1L):(end - 1L)]
  }
  run_grid <- function(quick, dev) {
    env <- .bias_uns_test_env()
    # Mirrors QMD set-up: dev takes precedence over quick.
    env$analysis_quick <- quick && !dev
    env$analysis_dev <- dev
    env$simulation_seed <- 12345L
    env$sim_grid_shuffle_seed <- 8L
    env$sim_grid_chunk_index <- 1L
    env$sim_grid_n_chunks <- 1L
    env$analysis_semantics_version <- "test"
    eval(parse(text = chunk("bias-uns-settings")), envir = env)
    invisible(utils::capture.output(
      eval(parse(text = chunk("bias-uns-grid")), envir = env)
    ))
    eval(parse(text = chunk("bias-uns-run-settings")), envir = env)
    env
  }
  full <- run_grid(FALSE, FALSE)
  expect_equal(nrow(full$sim_grid_all), 6720L)
  expect_equal(sort(unique(full$sim_grid_all$n_cell)), c(1e4, 5e4))
  expect_equal(sort(unique(full$sim_grid_all$prob_response)), c(0.002, 0.05))
  expect_equal(full$scenario_settings$nSample, 25)
  expect_equal(full$scenario_settings$biasUnsWidthHeightFrac, 0.15)
  expect_equal(full$scenario_settings$covEvMin, 1.5)
  expect_equal(full$scenario_settings$covEvMax, 1.5)
  expect_null(full$scenario_settings$tolClust)
  expect_false(full$scenario_settings$calcCytPosGates)
  expect_false(full$scenario_settings$locEnforceShapeThreshold)
  expect_identical(full$analysis_required_params$simulation_seed, 12345L)
  expect_identical(
    full$analysis_required_params$sim_grid_spec,
    dplyr::select(full$sim_grid_all, -sim_seed)
  )
  expect_equal(
    full$bias_uns_settings_tbl$bias_uns_multiplier[
      full$bias_uns_settings_tbl$bias_uns_basis == "bandwidth"
    ],
    c(0, 0.1, 0.25, 0.33, 0.5, 0.75, 1, 1.25, 1.5, 2)
  )
  expect_equal(
    full$bias_uns_settings_tbl$bias_uns_multiplier[
      full$bias_uns_settings_tbl$bias_uns_basis == "negative_width"
    ],
    c(0.1, 0.25, 0.5, 0.75)
  )
  for (transformation in c("gaussian", "skew", "gamma")) {
    rows <- full$sim_grid_all |>
      dplyr::filter(.data$transformation == .env$transformation)
    expect_equal(
      sort(unique(rows$mean_pos)),
      switch(transformation, gaussian = c(4.5, 8), skew = c(6, 8.5), gamma = c(4, 7))
    )
    expect_true(all(rows$background_relative_to_response == 0.2))
    expect_true(all(rows$n_cell_uns_relative_to_stim == 1))
    expect_equal(
      sort(unique(rows$bw)),
      if (transformation == "gamma") c(0.001, 0.0025, 0.005, 0.01) else {
        c(0.05, 0.1, 0.25, 0.5)
      }
    )
    shifted <- rows[rows$mismatch_type == "mean_shift", ]
    expect_equal(
      sort(unique(shifted$stim_mean_shift)),
      if (transformation == "gamma") c(0, 0.005, 0.01, 0.05) else {
        c(0, 0.025, 0.05, 0.2)
      }
    )
    expect_true(all(shifted$stim_mean_shift_clusters == "gn"))
    expect_true(all(shifted$stim_sd_multiplier == 1))
    inflated <- rows[rows$mismatch_type == "sd_inflation", ]
    expect_true(all(inflated$stim_mean_shift == 0))
    expect_true(all(inflated$stim_sd_multiplier == 1.10))
    expect_true(all(inflated$stim_sd_multiplier_clusters == "gn"))
  }
  grouped <- full$sim_grid_all |>
    dplyr::group_by(.data$base_scenario_id) |>
    dplyr::summarise(
      n_seed = dplyr::n_distinct(.data$sim_seed),
      n_rule = dplyr::n_distinct(.data$bias_uns_rule),
      n_mismatch = dplyr::n_distinct(.data$mismatch_label),
      n_bw = dplyr::n_distinct(.data$bw), .groups = "drop"
    )
  expect_true(all(grouped$n_seed == 1L))
  expect_true(all(grouped$n_rule == 14L))
  expect_true(all(grouped$n_mismatch == 5L))
  expect_true(all(grouped$n_bw == 4L))
  quick <- run_grid(TRUE, FALSE)
  expect_equal(nrow(quick$sim_grid_all), 120L)
  expect_setequal(quick$sim_grid_all$n_cell, c(1e4, 5e4))
  expect_setequal(quick$sim_grid_all$mismatch_label, c("mean shift 0", "SD inflation 10%"))
  expect_setequal(quick$sim_grid_all$bias_uns_basis, c("bandwidth", "negative_width"))
  for (transformation in unique(quick$sim_grid_all$transformation)) {
    expect_equal(dplyr::n_distinct(quick$sim_grid_all$bw[
      quick$sim_grid_all$transformation == transformation
    ]), 2L)
  }
  expect_identical(quick$scenario_settings$nSample, 1L)
  expect_identical(run_grid(TRUE, TRUE)$sim_grid_all, run_grid(FALSE, TRUE)$sim_grid_all)
  for (mode in list(c(TRUE, FALSE), c(FALSE, TRUE), c(TRUE, TRUE))) {
    reduced <- run_grid(mode[[1L]], mode[[2L]])
    expect_gt(nrow(reduced$sim_grid_all), 0L)
    expect_identical(reduced$sim_grid_full, full$sim_grid_full)
    expected <- full$sim_grid_full |>
      dplyr::filter(.data$sim_id %in% reduced$sim_grid_all$sim_id)
    expect_identical(reduced$sim_grid_all, expected)
  }
})

test_that("Analysis 2b runtime guards simulations and reads canonical results", {
  env <- .bias_uns_test_env()
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  lines <- readLines(file.path(root, "analysis", "2b-sim-bias_uns-freq_bs.qmd"))
  chunk <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1L]
    lines[(start + 1L):(end - 1L)]
  }
  unexpected <- function(...) stop("A disabled simulation or plot executed.")
  env$run_simulations <- FALSE
  env$run_plots <- FALSE
  env$.analysis_run_context <- unexpected
  env$.simBandwidthRunGrid <- unexpected
  env$.simBandwidthFinishChunk <- unexpected
  env$.analysis_cache_dir <- unexpected
  env$ggplot <- unexpected
  for (label in c(
    "bias-uns-parallel", "bias-uns-collate",
    "fig-relative-error", "bias-uns-signed-error", "fig-signed-error",
    "fig-estimated-frequency"
  )) {
    expect_no_error(eval(parse(text = chunk(label)), envir = env))
  }
  expect_false(exists("run_ctx", envir = env, inherits = FALSE))

  env$run_plots <- TRUE
  env$analysis_key <- c("sim", "bias_uns", "freq_bs")
  env$root_dir <- root
  env$analysis_required_params <- list(
    analysis_semantics_version = "bias-uns-freq-v1", simulation_seed = 12345L
  )
  env$analysis_qmd <- "analysis/2b-sim-bias_uns-freq_bs.qmd"
  env$.analysis_results_context <- function(analysis_key, path_root, qmd_path) {
    expect_identical(analysis_key, env$analysis_key)
    expect_identical(path_root, root)
    list(read_only = TRUE)
  }
  reads <- list()
  env$.analysis_read_current <- function(run_ctx, relative_path, required_params) {
    expect_true(run_ctx$read_only)
    expect_identical(required_params, env$analysis_required_params)
    reads[[length(reads) + 1L]] <<- relative_path
    relative_path[[2L]]
  }
  env$readRDS <- identity
  eval(parse(text = chunk("bias-uns-parallel")), envir = env)
  eval(parse(text = chunk("bias-uns-collate")), envir = env)
  expect_identical(reads, list(
    c("collated", "bias_uns_results_raw.rds"),
    c("collated", "bias_uns_results_summary.rds")
  ))

  rerun <- chunk("rerun-one-simulation")
  expect_true("#| eval: false" %in% rerun)
  expect_true(any(grepl(".simBandwidthRunRow(", rerun, fixed = TRUE)))
  expect_true(any(grepl(
    "scenario_fn = .simBandwidthBiasUnsScenario", rerun, fixed = TRUE
  )))
})

test_that("negative shoulder width matches a normal KDE and stops at an antimode", {
  env <- .bias_uns_test_env()
  x <- stats::qnorm(seq(0.0001, 0.9999, length.out = 10000L))
  width <- env$.simBandwidthNegativeShoulderWidth(x, bw = 0.25)
  expected <- sqrt(1 + 0.25^2) * sqrt(-2 * log(0.15))
  expect_equal(width, expected, tolerance = 0.04)

  mixture <- c(
    stats::qnorm(seq(0.0001, 0.9999, length.out = 6000L), mean = -3),
    stats::qnorm(seq(0.0001, 0.9999, length.out = 4000L), mean = 0)
  )
  width_15 <- env$.simBandwidthNegativeShoulderWidth(mixture, bw = 0.25)
  width_01 <- env$.simBandwidthNegativeShoulderWidth(
    mixture, bw = 0.25, heightFrac = 0.01
  )
  expect_true(is.finite(width_15))
  expect_gt(width_15, 1)
  expect_lt(width_15, 2.5)
  expect_identical(width_15, width_01)
  expect_true(is.na(env$.simBandwidthNegativeShoulderWidth(
    rep(1, 10), bw = 0.1
  )))
})

test_that("Analysis 2b retains invalid final frequency diagnostics", {
  env <- .bias_uns_test_env()
  tbl <- tibble::tibble(
    sim_id = 1L, sim_seed = 12L, iter = 1L,
    sample = c("1", "2", "3"), ind = c("2", "4", "6"),
    method = "loc_sample", propRespTruth = 0.1,
    propRespEst = c(0.1, 0, NA_real_), threshold = c(3, Inf, NA_real_),
    nCellStim = 1000L, nCellUns = 1000L,
    nPosStim = c(120L, 0L, NA_integer_),
    nPosUns = c(20L, 0L, NA_integer_),
    propStim = c(0.12, 0, NA_real_), propUns = c(0.02, 0, NA_real_),
    thresholdOrigin = c("calculated", "fallback", "failed"),
    gateReturnPoint = NA_character_, locGenerated = c(TRUE, FALSE, FALSE),
    locGeneratedDirect = c(TRUE, FALSE, FALSE),
    locSource = c("direct", "fallback", "fallback"),
    locReason = c("selected", "no_response", "unavailable"),
    biasUns = 0.025, biasUnsNegativeWidth = NA_real_
  )
  res <- env$.simBandwidthBiasUnsCollate(tbl, c("sim_id", "sim_seed"))
  expect_identical(res$bias_uns_results_raw$ind, tbl$ind)
  expect_identical(res$bias_uns_results_raw$threshold, tbl$threshold)
  expect_equal(res$bias_uns_results_summary$propRespEst_mean, 0.1)
  expect_equal(res$bias_uns_results_summary$n_sample, 3L)
  expect_equal(res$bias_uns_results_summary$n_valid, 1L)
  expect_equal(res$bias_uns_results_summary$n_failed, 2L)
  expect_equal(res$bias_uns_results_summary$failure_fraction, 2 / 3)
  expect_true(all(is.na(res$bias_uns_results_raw$rel_error[2:3])))

  invalid <- tbl
  invalid$threshold <- NA_real_
  all_invalid <- env$.simBandwidthBiasUnsCollate(
    invalid, c("sim_id", "sim_seed")
  )
  expect_equal(all_invalid$bias_uns_results_summary$n_valid, 0L)
  expect_true(is.na(all_invalid$bias_uns_results_summary$propRespEst_mean))
  expect_true(is.na(all_invalid$bias_uns_results_summary$q90_abs_rel_error))

  missing_scenario <- tbl[1L, ]
  missing_scenario$sim_id <- 2L
  missing_scenario$method <- "propRespPred"
  expect_error(
    env$.simBandwidthBiasUnsCollate(
      dplyr::bind_rows(tbl, missing_scenario), c("sim_id", "sim_seed")
    ),
    "loc_sample"
  )
  expect_error(
    env$.simBandwidthBiasUnsCollate(
      tbl, c("sim_id", "sim_seed"), n_sample_expected = 25L
    ),
    "sample"
  )
  expect_error(
    env$.simBandwidthBiasUnsCollate(
      dplyr::bind_rows(tbl, tbl[3L, ]), c("sim_id", "sim_seed")
    ),
    "duplicate result keys"
  )
})

test_that("Analysis 2b runs width-based bias with selective batch mismatch", {
  env <- .bias_uns_test_env()
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  source(file.path(root, "scripts", "r", "sim-compare-freq_bs.R"), local = env)
  simulate <- env$.simCompareSimCytExperiment
  captured <- list()
  env$.simCompareSimCytExperiment <- function(...) {
    out <- simulate(...)
    captured[[length(captured) + 1L]] <<- out
    out
  }
  withr::local_envvar(c(STIMGATE_INTERMEDIATE = NA_character_))
  row <- tibble::tibble(
    sim_id = 1L,
    sim_seed = 2028L,
    bias_uns_basis = "negative_width",
    bias_uns_multiplier = 0.5,
    bw = 0.25,
    n_cell = 500L,
    prob_response = 0.05,
    mean_pos = 8,
    transformation = "gaussian",
    background_relative_to_response = 0.2,
    n_cell_uns_relative_to_stim = 1,
    stim_mean_shift = 0.05,
    stim_sd_multiplier = 1,
    stim_mean_shift_clusters = "gn",
    stim_sd_multiplier_clusters = NA_character_
  )
  settings <- list(
    nSample = 2L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUnsWidthHeightFrac = 0.15,
    bwMin = "none",
    bwMax = "none",
    probExact = TRUE,
    covEvMin = 1.5,
    covEvMax = 1.5,
    tolClust = NULL,
    locEnforceShapeThreshold = FALSE,
    calcCytPosGates = FALSE
  )
  res <- env$.simBandwidthRunRow(
    row, scenario_fn = env$.simBandwidthBiasUnsScenario, settings = settings
  )
  final <- res |>
    dplyr::filter(.data$method == "loc_sample")
  expect_equal(nrow(final), 2L)
  expect_true(all(is.finite(final$biasUnsNegativeWidth)))
  expect_true(all(final$biasUnsNegativeWidth > 0))
  expect_equal(final$biasUns, 0.5 * final$biasUnsNegativeWidth)

  # Estimator settings and deterministic mismatch must preserve the random draws.
  baseline_row <- row
  baseline_row$stim_mean_shift <- 0
  baseline_row$bias_uns_multiplier <- 0.25
  baseline_row$bw <- 0.5
  baseline <- env$.simBandwidthRunRow(
    baseline_row,
    scenario_fn = env$.simBandwidthBiasUnsScenario,
    settings = settings
  )
  expect_equal(nrow(baseline[baseline$method == "loc_sample", ]), 2L)
  expect_equal(length(captured), 2L)
  expect_identical(captured[[1L]]$labelsList, captured[[2L]]$labelsList)
  for (ind in seq_len(4L)) {
    shifted <- flowCore::exprs(captured[[1L]]$flowFrameList[[ind]])
    original <- flowCore::exprs(captured[[2L]]$flowFrameList[[ind]])
    expected_shift <- if (ind %% 2L == 0L) {
      as.numeric(captured[[1L]]$labelsList[[ind]] == "gn") * 0.05
    } else {
      rep(0, nrow(original))
    }
    expect_equal(shifted[, 1L] - original[, 1L], expected_shift)
  }
})
