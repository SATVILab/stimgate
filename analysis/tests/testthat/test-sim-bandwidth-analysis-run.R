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

test_that("analysis 2 scenario rerun is identical whatever the prior RNG", {
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
  old_kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(old_kind)))

  RNGkind("L'Ecuyer-CMRG")
  set.seed(99)
  res_lecuyer <- env$.simBandwidthRunRow(
    row, env$.simBandwidthFreqBsGlobalScenario, settings
  )
  expect_identical(RNGkind()[[1]], "L'Ecuyer-CMRG")

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

test_that("grid run outputs match an interactive rerun regardless of chunking", {
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
