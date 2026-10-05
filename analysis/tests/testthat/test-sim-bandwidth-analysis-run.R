root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

.load_bw_run_env <- function() {
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

# Cheap scenario: two random draws shifted by the row's `a`; fails for the
# row with `a == settings$fail_a`.
.fake_scenario <- function(row, settings) {
  if (identical(row$a[[1]], settings$fail_a)) {
    stop("boom")
  }
  tibble::tibble(iter = 1:2, x = stats::rnorm(2) + row$a[[1]])
}

.fake_grid <- function() {
  tibble::tibble(a = c(1, 2, 3)) |>
    dplyr::mutate(
      sim_id = dplyr::row_number(),
      sim_seed = as.integer(100L + sim_id - 1L)
    )
}

.local_run_ctx <- function(env, run_id, chunk_index = 1L, n_chunks = 1L,
                           frame = parent.frame()) {
  tmp_project <- withr::local_tempdir(.local_envir = frame)
  withr::local_dir(tmp_project, .local_envir = frame)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")
  env$.analysis_run_context(
    analysis_key = c("sim", "test", "run"),
    run_id = run_id,
    params = list(analysis_semantics_version = "test-v1"),
    sim_grid_chunk_index = chunk_index,
    sim_grid_n_chunks = n_chunks
  )
}

test_that("analysis 2a scenario rerun is identical whatever the prior RNG", {
  withr::local_preserve_seed()
  originalRngKind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(originalRngKind)))
  env <- .load_bw_run_env()
  settings <- list(
    nSample = 2L, nMarker = 1L, nCondition = 2L, nCluster = 2L, nIter = 1L,
    bwMin = "none", bwMax = "none", probExact = TRUE,
    covEvMin = 1.5, covEvMax = 1.5, tolClust = NULL,
    locEnforceShapeThreshold = FALSE, calcCytPosGates = FALSE
  )
  row <- tibble::tibble(
    transformation = "gaussian", prob_response = 0.05, n_cell = 200,
    mean_pos_setting = "high", mean_pos = 8, bw = 0.25,
    bias_uns_setting = "low", bias_uns = 0.05,
    sample_perturbation_sd = 0, condition_perturbation_sd = 0,
    cluster_perturbation_sd = 0, background_relative_to_response = 0.2,
    n_cell_uns_relative_to_stim = 1, sim_id = 7L, sim_seed = 12351L
  )
  withr::local_seed(123L)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(99)
  seed_before <- .Random.seed
  res_lecuyer <- env$.simBandwidthRunRow(
    row, env$.simBandwidthFreqBsGlobalScenario, settings
  )
  expect_identical(RNGkind()[[1]], "L'Ecuyer-CMRG")
  expect_identical(.Random.seed, seed_before)

  RNGkind("default", "default", "default")
  set.seed(1)
  stats::runif(5)
  res_default <- env$.simBandwidthRunRow(
    row, env$.simBandwidthFreqBsGlobalScenario, settings
  )

  expect_identical(res_lecuyer, res_default)
  expect_true(all(is.na(res_default$error_message)))
  expect_identical(names(res_default)[seq_along(row)], names(row))
})

test_that("grid and interactive reruns agree regardless of chunking", {
  env <- .load_bw_run_env()
  grid <- .fake_grid()
  settings <- list(fail_a = NA)

  run_chunks <- function(n_chunks, order) {
    out <- list()
    for (k in seq_len(n_chunks)) {
      chunk_grid <- grid[order, ][seq_len(nrow(grid)) %% n_chunks == k - 1L, ]
      run_ctx <- .local_run_ctx(env, paste0("chunks-", n_chunks), k, n_chunks)
      out <- c(out, env$.simBandwidthRunGrid(
        chunk_grid,
        scenario_fn = .fake_scenario,
        settings = settings,
        run_ctx = run_ctx,
        sim_grid_chunk_index = k,
        sim_grid_n_chunks = n_chunks,
        workers = 1L
      ))
    }
    res <- purrr::list_rbind(out)
    res[order(res$sim_id, res$iter), ]
  }

  res_one <- run_chunks(1L, 1:3)
  res_two <- run_chunks(2L, c(3L, 1L, 2L))
  expect_identical(res_one, res_two)

  rerun <- env$.simBandwidthRunRow(grid[2, ], .fake_scenario, settings)
  expect_identical(rerun, res_one[res_one$sim_id == 2L, ])
})

test_that("error rows bind with success rows and resume retries them", {
  env <- .load_bw_run_env()
  grid <- .fake_grid()
  run_ctx <- .local_run_ctx(env, "retry")
  run_grid <- function(settings, retry_errors = TRUE) {
    env$.simBandwidthRunGrid(
      grid,
      scenario_fn = .fake_scenario,
      settings = settings,
      run_ctx = run_ctx,
      retry_errors = retry_errors,
      workers = 1L
    )
  }

  failed <- purrr::list_rbind(run_grid(list(fail_a = 2)))
  expect_identical(failed$error_message[failed$sim_id == 2L], "boom")
  expect_true(is.na(failed$x[failed$sim_id == 2L]))
  expect_type(failed$x, "double")
  expect_identical(
    env$.simBandwidthChunkMarkerCounts(run_ctx),
    list(completed = 2L, failed = 1L)
  )

  kept <- purrr::list_rbind(run_grid(list(fail_a = NA), retry_errors = FALSE))
  expect_identical(kept$error_message[kept$sim_id == 2L], "boom")

  retried <- purrr::list_rbind(run_grid(list(fail_a = NA)))
  expect_true(all(is.na(retried$error_message)))
  expect_identical(
    env$.simBandwidthChunkMarkerCounts(run_ctx),
    list(completed = 3L, failed = 0L)
  )
})

test_that("finishing a chunk validates, promotes and refuses errors", {
  env <- .load_bw_run_env()
  grid <- .fake_grid()
  collate_fn <- function(tbl) list(summary = dplyr::count(tbl, .data$sim_id))

  run_ctx <- .local_run_ctx(env, "promote-error")
  env$.simBandwidthRunGrid(
    grid,
    scenario_fn = .fake_scenario,
    settings = list(fail_a = 3),
    run_ctx = run_ctx,
    workers = 1L
  )
  expect_error(
    env$.simBandwidthFinishChunk(
      run_ctx, grid, grid,
      collate_fn = collate_fn, label = "test analysis"
    ),
    "simulation errors for sim_id: 3"
  )
  expect_false(file.exists(file.path(run_ctx$current_dir, "COMPLETE")))

  env$.simBandwidthRunGrid(
    grid,
    scenario_fn = .fake_scenario,
    settings = list(fail_a = NA),
    run_ctx = run_ctx,
    workers = 1L
  )
  expect_true(env$.simBandwidthFinishChunk(
    run_ctx, grid, grid,
    collate_fn = collate_fn, label = "test analysis"
  ))
  path <- env$.analysis_current_file(
    env$.analysis_results_context(run_ctx$analysis_key),
    c("collated", "summary.rds"),
    required_params = list(analysis_semantics_version = "test-v1")
  )
  expect_identical(readRDS(path)$sim_id, 1:3)
  expect_length(
    env$.simBandwidthFindSimOutput(run_ctx$current_dir, 2L),
    1L
  )
})

test_that("promotion waits for every chunk and checks the full grid", {
  env <- .load_bw_run_env()
  grid <- .fake_grid()
  collate_fn <- function(tbl) list(summary = tbl)

  run_ctx <- .local_run_ctx(env, "two-chunks", 1L, 2L)
  chunk_grid <- grid[1:2, ]
  env$.simBandwidthRunGrid(
    chunk_grid,
    scenario_fn = .fake_scenario,
    settings = list(fail_a = NA),
    run_ctx = run_ctx,
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 2L,
    workers = 1L
  )
  expect_false(env$.simBandwidthFinishChunk(
    run_ctx, chunk_grid, grid,
    collate_fn = collate_fn, label = "test analysis"
  ))

  validation <- env$.simBandwidthValidateOutputs(
    env$.simBandwidthReadOutputs(
      env$.find_bw_list_output_files(run_ctx$staging_run_dir)
    ),
    expected_grid = grid
  )
  expect_false(validation$ids_ok)
  expect_false(validation$validation_ok)
})


test_that("row RNG and its absence are restored even on scenario errors", {
  withr::local_preserve_seed()
  originalRngKind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(originalRngKind)))
  env <- .load_bw_run_env()
  withr::local_seed(10L)
  RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rejection")
  set.seed(18L)
  before_kind <- RNGkind()
  before_seed <- .Random.seed
  expect_error(
    env$.simBandwidthRunRow(.fake_grid()[1, ], .fake_scenario,
      list(fail_a = 1)
    ),
    "boom"
  )
  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before_seed)

  rm(".Random.seed", envir = .GlobalEnv)
  env$.simBandwidthRunRow(.fake_grid()[1, ], .fake_scenario, list(fail_a = NA))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(RNGkind(), before_kind)
  expect_error(
    env$.simBandwidthRunRow(.fake_grid()[1, ], function(row, settings) {
      tibble::tibble()
    }),
    "at least one result row"
  )
})

test_that("stale error markers and unreadable outputs are retried", {
  env <- .load_bw_run_env()
  grid <- .fake_grid()[1, ]
  ctx <- .local_run_ctx(env, "stale-error")
  run <- function(retry_errors = TRUE) {
    env$.simBandwidthRunRowResumable(
      grid, .fake_scenario, list(fail_a = NA), ctx, 1L,
      retry_errors = retry_errors
    )
  }
  reference <- run()
  output <- env$.path_sim_output(1L, ctx$chunk_output_dir, 1L, 1L)
  stale <- reference
  stale$x <- -999
  saveRDS(stale, output)
  file.create(file.path(ctx$chunk_jobs_dir, "error-1"))
  expect_identical(run(), reference)
  expect_false(file.exists(file.path(ctx$chunk_jobs_dir, "error-1")))
  writeLines("unreadable RDS", output)
  expect_identical(run(), reference)
})

test_that("full-grid validation catches wrong IDs, seeds and corrupt files", {
  env <- .load_bw_run_env()
  grid <- .fake_grid()
  tbl <- purrr::list_rbind(lapply(seq_len(nrow(grid)), function(i) {
    env$.simBandwidthRunRow(grid[i, ], .fake_scenario, list(fail_a = NA))
  }))
  expect_true(env$.simBandwidthValidateOutputs(tbl, grid)$validation_ok)
  wrong_seed <- tbl
  wrong_seed$sim_seed[1] <- 0L
  expect_match(
    env$.simBandwidthValidateOutputs(wrong_seed, grid)$problems,
    "sim_seed values"
  )
  extra <- tbl
  extra$sim_id[1] <- 99L
  expect_false(env$.simBandwidthValidateOutputs(extra, grid)$ids_ok)
  expect_identical(
    env$.simBandwidthValidateOutputs(tbl, grid,
      validate_fn = function(tbl) "analysis-specific failure"
    )$problems,
    "analysis-specific failure"
  )
  path <- withr::local_tempfile()
  writeLines("broken", path)
  expect_error(env$.simBandwidthReadOutputs(path), "Could not read saved")
})

test_that("all chunks collate together and empty chunks can complete", {
  env <- .load_bw_run_env()
  project <- withr::local_tempdir()
  withr::local_dir(project)
  grid <- .fake_grid()
  contexts <- lapply(1:2, function(k) {
    env$.analysis_run_context(
      c("sim", "test", "shared"), run_id = "shared",
      path_root = project, sim_grid_chunk_index = k, sim_grid_n_chunks = 2L
    )
  })
  collate <- function(tbl) list(summary = tbl)
  env$.simBandwidthRunGrid(
    grid, scenario_fn = .fake_scenario, settings = list(fail_a = NA),
    run_ctx = contexts[[1]], sim_grid_n_chunks = 2L, workers = 1L
  )
  expect_false(env$.simBandwidthFinishChunk(
    contexts[[1]], grid, grid, collate_fn = collate
  ))
  expect_true(env$.simBandwidthFinishChunk(
    contexts[[2]], grid[0, ], grid, collate_fn = collate
  ))
  expect_identical(
    list.files(contexts[[1]]$staging_collated_dir), "summary.rds"
  )
  expect_identical(
    sort(unique(readRDS(file.path(
      contexts[[1]]$current_dir, "collated", "summary.rds"
    ))$sim_id)),
    grid$sim_id
  )
})

test_that("parallel workers use the same row stream as an interactive rerun", {
  env <- .load_bw_run_env()
  grid <- .fake_grid()
  ctx <- .local_run_ctx(env, "parallel")
  actual <- env$.simBandwidthRunGrid(
    grid[3:1, ], scenario_fn = .fake_scenario, settings = list(fail_a = NA),
    run_ctx = ctx, workers = 2L
  )
  expected <- lapply(3:1, function(i) {
    env$.simBandwidthRunRow(grid[i, ], .fake_scenario, list(fail_a = NA))
  })
  expect_identical(actual, expected)
})

test_that("analysis 2 collates final sample estimates", {
  env <- .load_bw_run_env()
  grid <- .fake_grid()[1, ]
  tbl <- dplyr::bind_cols(
    grid[rep(1L, 3L), ],
    tibble::tibble(
      iter = 1L, sample = c("1", "1", "2"), ind = c("2", "2", "4"),
      method = c("propRespPred", "loc_sample", "loc_sample"),
      threshold = c(NA_real_, 3, 5), propRespTruth = 0.1,
      propRespEst = c(0.99, 0.15, 0.25), propBsEst = 0.95
    )
  )
  res <- env$.simBandwidthFreqBsGlobalCollate(tbl, names(grid))
  expect_identical(res$bw_tbl_results_raw$propRespEst, c(0.15, 0.25))
  expect_equal(res$bw_tbl_results_summary$propRespEst_median, 0.2)
  expect_false("propBsEst" %in% names(res$bw_tbl_results_raw))
  expect_error(
    env$.simBandwidthFreqBsGlobalCollate(
      dplyr::bind_rows(tbl, tbl), names(grid)
    ),
    "duplicate result keys"
  )
})

test_that("analysis 2a dev and quick filters preserve full-grid IDs and seeds", {
  lines <- readLines(file.path(
    root_dir, "analysis", "2a-sim-bw-freq_bs-global.qmd"
  ))
  chunk <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1]
    lines[(start + 1L):(end - 1L)]
  }
  run_grid <- function(quick, dev) {
    env <- .load_bw_run_env()
    # Mirrors QMD set-up: dev takes precedence over quick.
    env$analysis_quick <- quick && !dev
    env$analysis_dev <- dev
    env$simulation_seed <- 12345L
    env$sim_grid_shuffle_seed <- 8L
    env$sim_grid_chunk_index <- 1L
    env$sim_grid_n_chunks <- 1L
    eval(parse(text = chunk("actual-settings")), envir = env)
    invisible(utils::capture.output(
      eval(parse(text = chunk("bw-manual-grid")), envir = env)
    ))
    env
  }
  full <- run_grid(FALSE, FALSE)
  quick <- run_grid(TRUE, FALSE)
  dev <- run_grid(FALSE, TRUE)
  both <- run_grid(TRUE, TRUE)
  for (env in list(quick, dev, both)) {
    expect_identical(env$sim_grid_full, full$sim_grid_full)
    expect_gt(nrow(env$sim_grid_all), 0L)
    expected <- full$sim_grid_full |>
      dplyr::filter(.data$sim_id %in% env$sim_grid_all$sim_id)
    expect_identical(env$sim_grid_all, expected)
  }
  expect_equal(nrow(quick$sim_grid_all), 216L)
  expect_setequal(quick$sim_grid_all$n_cell, c(1e3, 5e3))
  # Only matched conditions (no condition perturbation) are simulated.
  expect_setequal(full$sim_grid_full$condition_perturbation_sd, 0)
  expect_identical(both$sim_grid_all, dev$sim_grid_all)
  expect_true(0.02 %in% quick$sim_grid_all$prob_response)
  expect_identical(nrow(both$sim_grid_all), 3L)
})

test_that("shared run helpers cannot write to a read-only results context", {
  env <- .load_bw_run_env()
  project <- withr::local_tempdir()
  withr::local_dir(project)
  ctx <- list(
    read_only = TRUE,
    sim_root = file.path(project, "cache", "sim", "absent")
  )
  grid <- .fake_grid()
  expect_false(env$.simBandwidthFinishChunk(
    ctx, grid, grid, collate_fn = function(tbl) list(summary = tbl)
  ))
  expect_false(env$.simBandwidthPromoteIfReady(
    ctx, grid, collate_fn = function(tbl) list(summary = tbl)
  ))
  expect_error(
    env$.simBandwidthRunRowResumable(
      grid[1, ], .fake_scenario, list(fail_a = NA), ctx, 1L
    ),
    "read-only results context"
  )
  expect_false(dir.exists(ctx$sim_root))
})


test_that("analysis 2b bias rules forward the intended realised settings", {
  env <- .load_bw_run_env()
  captured <- NULL
  env$.simBandwidthBsFreq <- function(...) {
    captured <<- list(...)
    tibble::tibble(method = "loc_sample")
  }

  row <- tibble::tibble(
    bias_uns_basis = "bandwidth",
    bias_uns_multiplier = 0.25,
    bw = 0.1,
    n_cell = 1e4,
    prob_response = 0.002,
    mean_pos = 8,
    transformation = "gaussian",
    background_relative_to_response = 0.2,
    n_cell_uns_relative_to_stim = 1,
    stim_mean_shift = 0.05,
    stim_sd_multiplier = 1,
    stim_mean_shift_clusters = "gn",
    stim_sd_multiplier_clusters = NA_character_
  )
  env$.simBandwidthBiasUnsScenario(row, list(nSample = 25))
  expect_equal(captured$biasUns, 0.025)
  expect_null(captured$biasUnsWidthMultiplier)
  expect_equal(captured$stimMeanShift, 0.05)
  expect_identical(captured$stimMeanShiftClusters, "gn")

  row$bias_uns_basis <- "negative_width"
  row$bias_uns_multiplier <- 0.5
  env$.simBandwidthBiasUnsScenario(row, list(nSample = 25))
  expect_equal(captured$biasUns, 0)
  expect_equal(captured$biasUnsWidthMultiplier, 0.5)
})

test_that("analysis 2b declares the agreed grid and common-random-number seeds", {
  content <- paste(readLines(file.path(
    root_dir, "analysis", "2b-sim-bias_uns-freq_bs.qmd"
  ), warn = FALSE), collapse = "\n")
  has <- function(x) grepl(x, content, fixed = TRUE)

  expect_true(has("n_cell_vec <- c(1e4, 5e4)"))
  expect_true(has("prob_response_vec <- c(0.002, 0.05)"))
  expect_true(has(
    "n_sample_sim <- if (analysis_quick) 1L else c(final = 25L, draft = 10L)[[sim_size]]"
  ))
  expect_true(has("nSample = n_sample_sim"))
  expect_true(has("if (nrow(sim_grid) != 6720L)"))
  expect_true(has("biasUnsWidthHeightFrac = 0.15"))
  expect_true(has("stim_mean_shift_clusters = \"gn\""))
  expect_true(has("stim_sd_multiplier = 1.10"))
  expect_true(has("sim_seed = as.integer(simulation_seed + .data$base_scenario_id - 1L)"))
  expect_true(has("scenario_fn = .simBandwidthBiasUnsScenario"))
  expect_true(has(".simBandwidthBiasUnsCollate("))
})


test_that("negative shoulder width extends as the density-height cutoff falls", {
  env <- .load_bw_run_env()
  x <- stats::qnorm(seq(0.001, 0.999, length.out = 2000))
  width_50 <- env$.simBandwidthNegativeShoulderWidth(
    x,
    bw = 0.25,
    heightFrac = 0.50
  )
  width_15 <- env$.simBandwidthNegativeShoulderWidth(
    x,
    bw = 0.25,
    heightFrac = 0.15
  )

  expect_true(is.finite(width_50))
  expect_true(is.finite(width_15))
  expect_gt(width_15, width_50)
})
