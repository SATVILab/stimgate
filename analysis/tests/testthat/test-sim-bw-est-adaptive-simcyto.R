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

test_that("analysis/5-sim-bw-est-adaptive.qmd does not source functionsForBenchmarking-Cyt.R", {
  qmd_path <- file.path(root_dir, "analysis", "5-sim-bw-est-adaptive.qmd")
  expect_true(file.exists(qmd_path))

  lines <- readLines(qmd_path, warn = FALSE)
  expect_false(
    any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
    info = "analysis/5-sim-bw-est-adaptive.qmd should not source functionsForBenchmarking-Cyt.R"
  )
})

test_that(".simBandwidthEstBwDirectAdaptive preserves simcyto simulation boundary and adaptive outputs", {
  # Fixed expected values need R's default generators; an earlier test
  # can leave the parallel (L'Ecuyer) generator selected.
  withr::local_seed(1L, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  withr::local_preserve_seed()
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)

  run_case <- function(
      seed,
      transformation,
      mean_pos,
      bias_uns,
      expected_bw_means) {
    n_sample <- 2L
    n_condition <- 2L
    n_cell_stim <- 240L
    n_cell_uns <- 240L
    prob_response <- 0.05
    background_relative_to_response <- 0.2
    prob_response_uns <- prob_response * background_relative_to_response

    captured_args <- NULL
    captured_sim <- NULL
    orig_simcyto_experiment <- simcyto::simCytExperiment

    set.seed(seed)
    res <- testthat::with_mocked_bindings(
      simCytExperiment = function(...) {
        captured_args <<- list(...)
        out <- orig_simcyto_experiment(...)
        captured_sim <<- out
        out
      },
      .package = "simcyto",
      env$.simBandwidthEstBwDirectAdaptive(
        nSample = n_sample,
        nMarker = 1L,
        nCondition = n_condition,
        nCluster = 2L,
        nIter = 1L,
        biasUns = bias_uns,
        bwMtd = "hpi1Norm",
        bwFallback = 0.234,
        bwMin = -Inf,
        bwMax = Inf,
        bwNcellMax = 500L,
        nCellStim = n_cell_stim,
        probResponse = prob_response,
        probExact = TRUE,
        meanPos = mean_pos,
        transformation = transformation,
        backgroundRelativeToResponse = background_relative_to_response,
        ncellUnsRelativeToStim = 1,
        covEvMin = 1.5,
        covEvMax = 1.5,
        summarise = FALSE
      )
    )

    expect_s3_class(res, "tbl_df")
    expect_equal(nrow(res), n_sample)

    expect_type(captured_args, "list")
    expect_equal(captured_args$nCellByCondition, c(n_cell_uns, n_cell_stim))
    expect_equal(captured_args$meanExprMat, matrix(c(0, mean_pos), byrow = TRUE, ncol = 1))
    expect_equal(captured_args$clusterLabelVec, c("gn", "gp"))
    expect_equal(captured_args$probVecUns, c(1 - prob_response_uns, prob_response_uns))
    expect_equal(captured_args$probResponseVecByStimCondition, list(c(-prob_response, prob_response)))
    expect_equal(captured_args$samplePerturbationSd, 0)
    expect_equal(captured_args$conditionPerturbationSd, 0)
    expect_equal(captured_args$clusterPerturbationSd, 0)
    expect_equal(captured_args$covEvMin, 1.5)
    expect_equal(captured_args$covEvMax, 1.5)

    expected_trans <- env$.simBandwidthGetTrans(transformation)
    expect_equal(captured_args$transformationFunc(c(-1, 0, 1)), expected_trans(c(-1, 0, 1)))

    set.seed(seed)
    # The helper draws each replicate's seed before simulating.
    set.seed(sample.int(.Machine$integer.max, 1L, replace = TRUE))
    sim_direct <- do.call(simcyto::simCytExperiment, captured_args)

    truth_from_labels <- function(sim_obj) {
      purrr::map_df(seq_len(n_sample), function(sample_curr) {
        ind_uns <- (sample_curr - 1L) * n_condition + 1L
        ind_stim <- ind_uns + 1L
        labels_uns <- sim_obj$labelsList[[ind_uns]]
        labels_stim <- sim_obj$labelsList[[ind_stim]]
        prop_uns_truth <- sum(labels_uns %in% "gp") / length(labels_uns)
        prop_stim_truth <- sum(labels_stim %in% "gp") / length(labels_stim)

        tibble::tibble(
          sample = as.character(sample_curr),
          ind = as.character(ind_stim),
          propStimTruth = prop_stim_truth,
          propUnsTruth = prop_uns_truth,
          propRespTruth = prop_stim_truth - prop_uns_truth
        )
      }) |>
        dplyr::arrange(.data$sample, .data$ind)
    }

    expect_equal(truth_from_labels(captured_sim), truth_from_labels(sim_direct), tolerance = 1e-12)

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
      sim_direct$flowFrameList,
      function(ff) mean(flowCore::exprs(ff)[, 1]),
      numeric(1)
    )
    expr_sds_direct <- vapply(
      sim_direct$flowFrameList,
      function(ff) stats::sd(flowCore::exprs(ff)[, 1]),
      numeric(1)
    )
    expect_equal(unname(expr_means_helper), unname(expr_means_direct), tolerance = 1e-12)
    expect_equal(unname(expr_sds_helper), unname(expr_sds_direct), tolerance = 1e-12)

    expect_equal(res$n_cell_stim, rep(n_cell_stim, n_sample))
    expect_equal(res$n_cell_uns, rep(n_cell_uns, n_sample))
    expect_true(all(res$n_uns_bw_core <= n_cell_uns))
    expect_true(all(res$n_stim_bw_core <= n_cell_stim))

    bw_cols <- c("bw_uns_core", "bw_stim_core", "bw_uns_extra", "bw_stim_extra")
    expect_true(all(vapply(res[bw_cols], function(x) all(is.finite(x)), logical(1))))
    expect_true(all(vapply(res[bw_cols], function(x) all(x > 0), logical(1))))

    bw_means <- colMeans(res[bw_cols], na.rm = TRUE)
    expect_equal(unname(bw_means), expected_bw_means, tolerance = 1e-7)
  }

  run_case(
    seed = 2026L,
    transformation = "gaussian",
    mean_pos = 8,
    bias_uns = 0.05,
    expected_bw_means = c(
      0.3370803958,
      0.3284741467,
      0.3999293279,
      0.4726434710
    )
  )

  run_case(
    seed = 2027L,
    transformation = "gamma",
    mean_pos = 4,
    bias_uns = 0.0025,
    expected_bw_means = c(
      0.0049646913,
      0.0053038529,
      0.0048982453,
      0.0210311014
    )
  )
})


test_that("adaptive normalised bandwidths are controlled by normAdaptiveNcell, not bwNcellMax", {
  withr::local_preserve_seed()
  set.seed(519L)
  x <- c(
    stats::rnorm(800L, mean = 0, sd = 1),
    stats::rnorm(200L, mean = 5, sd = 1)
  )

  calc_adaptive <- function(bw_ncell_max) {
    set.seed(520L)
    stimgate:::.bwCalcOne(
      x = x,
      bwMtd = "hpi1Norm",
      bwNcellMax = bw_ncell_max,
      normAdaptiveNcell = 200L,
      adaptive = TRUE
    )
  }

  bw_small_cap <- calc_adaptive(100L)
  bw_large_cap <- calc_adaptive(100000L)

  expect_true(is.list(bw_small_cap))
  expect_true(isTRUE(attr(bw_small_cap, "adaptive")))
  expect_equal(bw_small_cap$bwCore, bw_large_cap$bwCore, tolerance = 1e-12)
  expect_equal(bw_small_cap$bwExtra, bw_large_cap$bwExtra, tolerance = 1e-12)
  expect_equal(bw_small_cap$bw, bw_large_cap$bw, tolerance = 1e-12)
})

.load_adaptive_run_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R",
    "sim-misc.R",
    "sim-bandwidth.R",
    "sim-bandwidth-analysis-io.R",
    "sim-bandwidth-analysis-run.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  env
}

.adaptive_fake_result <- function(sim_id, n_rows, bw = 0.5) {
  tibble::tibble(
    sim_id = sim_id,
    transformation = "gaussian",
    bw_mtd = "hpi1Norm",
    n_cell = 1000,
    bw_uns_core = bw,
    bw_stim_core = c(bw, rep(NA_real_, n_rows - 1L)),
    bw_uns_extra = bw,
    bw_stim_extra = bw
  )
}

test_that("analysis 5 scenario forwards grid row and fixed settings", {
  env <- .load_adaptive_run_env()
  captured <- NULL
  env$.simBandwidthEstBwDirectAdaptive <- function(...) {
    captured <<- list(...)
    tibble::tibble(iter = 1L, bw_uns_core = 1)
  }
  row <- tibble::tibble(
    transformation = "gamma", prob_response = 0.2, n_cell = 5000,
    mean_pos = 4, bias_uns = 0.00125, bw_mtd = "hpi2Norm"
  )
  env$.simBandwidthEstAdaptiveScenario(
    row, list(nSample = 2L, normAdaptiveNcell = 2500L)
  )
  expect_identical(captured$nSample, 2L)
  expect_identical(captured$normAdaptiveNcell, 2500L)
  expect_identical(captured$bwMtd, "hpi2Norm")
  expect_identical(captured$transformation, "gamma")
  expect_identical(captured$nCellStim, 5000)
  expect_identical(captured$biasUns, 0.00125)
})

test_that("analysis 5 validation checks sample-row counts per sim_id", {
  env <- .load_adaptive_run_env()
  tbl <- dplyr::bind_rows(
    .adaptive_fake_result(1L, 2L), .adaptive_fake_result(2L, 2L)
  )
  expect_identical(env$.simBandwidthEstAdaptiveValidate(tbl, 2L), character())
  bw_cols <- c("bw_uns_core", "bw_stim_core", "bw_uns_extra", "bw_stim_extra")
  tbl[bw_cols] <- NA_real_
  expect_identical(env$.simBandwidthEstAdaptiveValidate(tbl, 2L), character())
  zero <- env$.simBandwidthEstAdaptiveCollate(tbl, c("sim_id", "bw_mtd"))$bw_tbl_results
  expect_true(all(zero$n_est == 0L))
  expect_true(all(zero$prop_est == 0))
  expect_true(all(is.na(zero$mean_bw)))
  expect_match(
    env$.simBandwidthEstAdaptiveValidate(tbl, 3L),
    "sim_id: 1, 2"
  )
})

test_that("analysis 5 collation summarises finite estimates per component", {
  env <- .load_adaptive_run_env()
  tbl <- .adaptive_fake_result(1L, 2L)
  res <- env$.simBandwidthEstAdaptiveCollate(tbl, c("sim_id", "bw_mtd"))
  expect_named(res, c("bw_list_raw", "bw_tbl_results"))
  expect_identical(res$bw_list_raw, tbl)
  out <- res$bw_tbl_results
  expect_identical(nrow(out), 4L)
  expect_false("bw_mtd_norm" %in% names(out))
  expect_identical(unique(out$bw_mtd_base), "hpi1")
  stim_core <- out[out$bw_component == "core" & out$bw_condition == "stim", ]
  expect_identical(stim_core$n_total, 2L)
  expect_identical(stim_core$n_est, 1L)
  expect_equal(stim_core$prop_est, 0.5)
})

test_that("analysis 5 dev filter selects values that exist in the grid", {
  lines <- readLines(file.path(
    root_dir, "analysis", "5-sim-bw-est-adaptive.qmd"
  ))
  chunk <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1]
    lines[(start + 1L):(end - 1L)]
  }
  run_grid <- function(dev) {
    env <- .load_adaptive_run_env()
    env$analysis_dev <- dev
    env$simulation_seed <- 12345L
    env$sim_grid_shuffle_seed <- 8L
    env$sim_grid_chunk_index <- 1L
    env$sim_grid_n_chunks <- 1L
    eval(parse(text = chunk("actual-settings")), envir = env)
    invisible(utils::capture.output(
      eval(parse(text = chunk("bw-estimate-grid")), envir = env)
    ))
    env
  }
  full <- run_grid(FALSE)
  dev <- run_grid(TRUE)
  expect_identical(
    full$n_cell_stim_vec, c(1e3, 5e3, 2e4, 1e5)
  )
  expect_identical(nrow(dev$sim_grid_all), 6L)
  expect_setequal(dev$sim_grid_all$n_cell, c(5e3, 2e4))
  expected <- full$sim_grid_full |>
    dplyr::filter(.data$sim_id %in% dev$sim_grid_all$sim_id)
  expect_identical(dev$sim_grid_all, expected)
})


test_that("analysis 5 plot filenames work with the controlled scenario labels", {
  lines <- readLines(file.path(root_dir, "analysis", "5-sim-bw-est-adaptive.qmd"))
  start <- grep("^          path_plot <- file.path", lines)
  end <- grep("^          dir.create", lines)
  expect_length(start, 1L)
  expect_length(end, 1L)
  env <- new.env(parent = baseenv())
  env$dir_results <- "plots"
  env$mean_pos_setting_curr <- factor("high", levels = c("low", "high"))
  env$bias_uns_setting_curr <- factor("low", levels = c("none", "low", "high"))
  env$bw_component_curr <- "core"
  env$bw_condition_curr <- "stim"
  expect_no_error(eval(parse(text = lines[start:(end - 1L)]), envir = env))
  expect_identical(
    env$path_plot,
    file.path("plots", "bw_estimate_meanpos_setting_high_bias_uns_low_component_core_condition_stim.png")
  )
})
