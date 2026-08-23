root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
script_runtime <- file.path(root_dir, "scripts", "r", "analysis-runtime.R")

.load_runtime_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_runtime, local = env)
  env
}

test_that("QMD param lookup follows param > default precedence and env override precedence", {
  env <- .load_runtime_env()
  env$params <- list(sim_grid_chunk_index = 7L)

  expect_identical(env$.get_qmd_param("sim_grid_chunk_index", 1L), 7L)
  expect_identical(env$.get_qmd_param("missing_param", "fallback"), "fallback")

  old_chunk_index <- Sys.getenv("SIM_GRID_CHUNK_INDEX", unset = NA_character_)
  old_n_chunks <- Sys.getenv("SIM_GRID_N_CHUNKS", unset = NA_character_)
  on.exit({
    if (is.na(old_chunk_index)) {
      Sys.unsetenv("SIM_GRID_CHUNK_INDEX")
    } else {
      Sys.setenv(SIM_GRID_CHUNK_INDEX = old_chunk_index)
    }
    if (is.na(old_n_chunks)) {
      Sys.unsetenv("SIM_GRID_N_CHUNKS")
    } else {
      Sys.setenv(SIM_GRID_N_CHUNKS = old_n_chunks)
    }
  }, add = TRUE)
  Sys.setenv(SIM_GRID_CHUNK_INDEX = "9", SIM_GRID_N_CHUNKS = "4")

  expect_identical(
    env$.get_qmd_param_env("sim_grid_chunk_index", "SIM_GRID_CHUNK_INDEX", 1L),
    "9"
  )
  expect_identical(
    env$.get_qmd_param_env("sim_grid_chunk_index", "SIM_GRID_MISSING", 1L),
    7L
  )
})

test_that("boolean flags parse consistently with standard QMD semantics", {
  env <- .load_runtime_env()

  expect_true(env$.as_flag(TRUE))
  expect_true(env$.as_flag("TRUE"))
  expect_true(env$.as_flag("yes"))
  expect_true(env$.as_flag("1"))
  expect_false(env$.as_flag("false"))
  expect_false(env$.as_flag(0))
})

test_that("sim grid chunk validation rejects invalid settings and formats labels", {
  env <- .load_runtime_env()

  expect_error(
    env$.validate_sim_grid_chunk_settings(0L, 3L),
    "Invalid sim grid chunk settings"
  )
  expect_error(
    env$.validate_sim_grid_chunk_settings(4L, 3L),
    "Invalid sim grid chunk settings"
  )
  expect_error(
    env$.validate_sim_grid_chunk_settings(NA_integer_, 3L),
    "Invalid sim grid chunk settings"
  )

  out <- env$.validate_sim_grid_chunk_settings(2L, 5L)
  expect_identical(out$sim_grid_chunk_index, 2L)
  expect_identical(out$sim_grid_n_chunks, 5L)
  expect_identical(env$.sim_chunk_label(2L, 5L), "002-of-005")
})

test_that("atomic RDS writes are readable and preserve object contents", {
  env <- .load_runtime_env()

  path <- tempfile(file.path(tempdir(), "analysis-runtime-"), fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)

  obj <- list(
    x = c(1, 2, 3),
    y = data.frame(a = 1:2, b = c("p", "q"))
  )

  expect_identical(env$.write_rds_atomic(obj, path), path)
  expect_equal(readRDS(path), obj)
})

test_that("run contexts are isolated by logical run ID", {
  env <- .load_runtime_env()

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  ctx_a <- env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-test"),
    run_id = "run-a",
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 1L
  )
  ctx_b <- env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-test"),
    run_id = "run-b",
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 1L
  )

  expect_false(identical(ctx_a$staging_run_dir, ctx_b$staging_run_dir))
  expect_false(identical(ctx_a$progress_run_dir, ctx_b$progress_run_dir))
  expect_true(file.exists(ctx_a$manifest_path))
  expect_true(file.exists(ctx_b$manifest_path))
  expect_true(file.exists(ctx_a$status_path))
  expect_true(file.exists(ctx_b$status_path))
})

test_that("promotion only occurs after all expected chunks complete and validate", {
  env <- .load_runtime_env()

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  ctx_1 <- env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-promotion"),
    run_id = "shared-run",
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 2L
  )
  saveRDS(tibble::tibble(x = 1L), file.path(ctx_1$staging_collated_dir, "chunk1.rds"))

  env$.analysis_mark_chunk(
    run_ctx = ctx_1,
    total_sims = 1L,
    completed_sims = 1L,
    failed_sims = 0L,
    collate_ok = TRUE,
    validation_ok = TRUE
  )

  expect_false(env$.analysis_can_promote(ctx_1))
  expect_false(isTRUE(env$.analysis_promote_run(ctx_1)))
  expect_false(dir.exists(ctx_1$current_dir))

  ctx_2 <- env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-promotion"),
    run_id = "shared-run",
    sim_grid_chunk_index = 2L,
    sim_grid_n_chunks = 2L
  )
  saveRDS(tibble::tibble(x = 2L), file.path(ctx_2$staging_collated_dir, "chunk2.rds"))

  env$.analysis_mark_chunk(
    run_ctx = ctx_2,
    total_sims = 1L,
    completed_sims = 1L,
    failed_sims = 0L,
    collate_ok = TRUE,
    validation_ok = TRUE
  )

  expect_true(env$.analysis_can_promote(ctx_2))
  expect_true(isTRUE(env$.analysis_promote_run(ctx_2)))
  expect_true(dir.exists(ctx_2$current_dir))
  expect_true(file.exists(file.path(ctx_2$staging_run_dir, "COMPLETE")))

  status <- env$.analysis_read_status(ctx_2)
  expect_identical(status$status, "completed")
  expect_true(isTRUE(status$promotion_done))
})

test_that("failed or invalid runs are never promoted", {
  env <- .load_runtime_env()

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  ctx <- env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-failure"),
    run_id = "failed-run",
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 1L
  )

  env$.analysis_mark_chunk(
    run_ctx = ctx,
    total_sims = 1L,
    completed_sims = 0L,
    failed_sims = 1L,
    collate_ok = FALSE,
    validation_ok = FALSE,
    error_message = "sim failure"
  )

  expect_false(env$.analysis_can_promote(ctx))
  expect_false(isTRUE(env$.analysis_promote_run(ctx)))
  expect_false(dir.exists(ctx$current_dir))
})

test_that("failed simulations block promotion even when collation and validation succeed", {
  env <- .load_runtime_env()

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  ctx <- env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-failed-sim-blocks-promotion"),
    run_id = "failed-sim-run",
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 1L
  )

  env$.analysis_mark_chunk(
    run_ctx = ctx,
    total_sims = 2L,
    completed_sims = 1L,
    failed_sims = 1L,
    collate_ok = TRUE,
    validation_ok = TRUE
  )

  expect_false(env$.analysis_can_promote(ctx))
  expect_false(isTRUE(env$.analysis_promote_run(ctx)))

  status <- env$.analysis_read_status(ctx)
  expect_identical(status$n_completed, 1L)
  expect_identical(status$n_failed, 1L)
  expect_identical(status$n_outstanding, 0L)
})

test_that("explicit run ID reuse rejects incompatible manifest settings", {
  env <- .load_runtime_env()

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-manifest"),
    run_id = "shared-manifest-run",
    params = list(sim_grid_n_chunks = 2L, sim_grid_shuffle_seed = 7L),
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 2L
  )

  expect_error(
    env$.analysis_run_context(
      analysis_key = c("sim", "analysis-runtime-manifest"),
      run_id = "shared-manifest-run",
      params = list(sim_grid_n_chunks = 2L, sim_grid_shuffle_seed = 99L),
      sim_grid_chunk_index = 2L,
      sim_grid_n_chunks = 2L
    ),
    "incompatible with existing manifest"
  )
})

test_that("concurrent chunk updates preserve per-chunk status and aggregate counts", {
  env <- .load_runtime_env()

  tmp_project <- tempfile("analysis-runtime-concurrent-")
  dir.create(tmp_project, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  setwd(tmp_project)
  on.exit(setwd(old_wd), add = TRUE)
  on.exit(unlink(tmp_project, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(c("directories:", "  docs:", "    path: docs"), file.path(tmp_project, "_projr.yml"))

  run_id <- "concurrent-run"
  analysis_key <- c("sim", "analysis-runtime-concurrent")
  script_path <- script_runtime

  cl <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(
    cl,
    varlist = c("tmp_project", "run_id", "analysis_key", "script_path"),
    envir = environment()
  )

  res <- parallel::clusterApply(cl, 1:2, function(i) {
    local_env <- new.env(parent = baseenv())
    setwd(tmp_project)
    source(script_path, local = local_env)
    ctx <- local_env$.analysis_run_context(
      analysis_key = analysis_key,
      run_id = run_id,
      sim_grid_chunk_index = as.integer(i),
      sim_grid_n_chunks = 2L
    )
    local_env$.analysis_mark_chunk(
      run_ctx = ctx,
      total_sims = 1L,
      completed_sims = 1L,
      failed_sims = 0L,
      collate_ok = TRUE,
      validation_ok = TRUE
    )
    TRUE
  })
  expect_true(all(unlist(res)))

  ctx_main <- env$.analysis_run_context(
    analysis_key = analysis_key,
    run_id = run_id,
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 2L
  )
  status <- env$.analysis_read_status(ctx_main)

  expect_setequal(names(status$chunks), c("001-of-002", "002-of-002"))
  expect_identical(status$n_completed, 2L)
  expect_identical(status$n_failed, 0L)
  expect_identical(status$n_outstanding, 0L)
})

test_that("resuming from existing output reconciles missing completion markers", {
  env <- .load_runtime_env()

  tmp_dir <- withr::local_tempdir()
  file_completed <- file.path(tmp_dir, "completed-1")
  file_error <- file.path(tmp_dir, "error-1")
  file_running <- file.path(tmp_dir, "running-1")
  file.create(file_running)

  ok_output <- tibble::tibble(sim_id = 1L, error_message = NA_character_)
  env$.analysis_reconcile_resume_markers(
    existing_output = ok_output,
    file_completed = file_completed,
    file_error = file_error,
    file_running = file_running,
    error_col = "error_message"
  )

  expect_true(file.exists(file_completed))
  expect_false(file.exists(file_error))
  expect_false(file.exists(file_running))
})

test_that("concurrent promotion attempts are serialised safely", {
  env <- .load_runtime_env()

  tmp_project <- tempfile("analysis-runtime-promote-concurrent-")
  dir.create(tmp_project, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  setwd(tmp_project)
  on.exit(setwd(old_wd), add = TRUE)
  on.exit(unlink(tmp_project, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(c("directories:", "  docs:", "    path: docs"), file.path(tmp_project, "_projr.yml"))

  run_id <- "promote-concurrent-run"
  analysis_key <- c("sim", "analysis-runtime-promote-concurrent")
  script_path <- script_runtime

  ctx_1 <- env$.analysis_run_context(
    analysis_key = analysis_key,
    run_id = run_id,
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 2L
  )
  ctx_2 <- env$.analysis_run_context(
    analysis_key = analysis_key,
    run_id = run_id,
    sim_grid_chunk_index = 2L,
    sim_grid_n_chunks = 2L
  )

  saveRDS(tibble::tibble(x = 1L), file.path(ctx_1$chunk_output_dir, "bw_list_raw-chunk_001-of_002-sim_id_000001.rds"))
  saveRDS(tibble::tibble(x = 2L), file.path(ctx_2$chunk_output_dir, "bw_list_raw-chunk_002-of_002-sim_id_000002.rds"))

  env$.analysis_mark_chunk(
    run_ctx = ctx_1,
    total_sims = 1L,
    completed_sims = 1L,
    failed_sims = 0L,
    collate_ok = TRUE,
    validation_ok = TRUE
  )
  env$.analysis_mark_chunk(
    run_ctx = ctx_2,
    total_sims = 1L,
    completed_sims = 1L,
    failed_sims = 0L,
    collate_ok = TRUE,
    validation_ok = TRUE
  )

  expect_true(env$.analysis_can_promote(ctx_1))

  cl <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(
    cl,
    varlist = c("tmp_project", "run_id", "analysis_key", "script_path"),
    envir = environment()
  )

  res <- parallel::clusterApply(cl, 1:2, function(i) {
    local_env <- new.env(parent = baseenv())
    setwd(tmp_project)
    source(script_path, local = local_env)
    ctx <- local_env$.analysis_run_context(
      analysis_key = analysis_key,
      run_id = run_id,
      sim_grid_chunk_index = as.integer(i),
      sim_grid_n_chunks = 2L
    )
    tryCatch(
      as.logical(local_env$.analysis_promote_run(ctx)),
      error = function(e) FALSE
    )
  })

  expect_true(any(unlist(res)))
  expect_true(dir.exists(ctx_1$current_dir))
  current_manifest <- readRDS(file.path(ctx_1$current_dir, "manifest.rds"))
  expect_identical(current_manifest$run_id, run_id)

  status <- env$.analysis_read_status(ctx_1)
  expect_true(isTRUE(status$promotion_done))
  expect_identical(status$status, "completed")
})

test_that("explicit run ID reuses the original dated run directory", {
  env <- .load_runtime_env()

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  run_id <- "resume-cross-date"
  ctx_initial <- env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-date-reuse"),
    run_id = run_id,
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 2L
  )

  target_date <- "1999-12-31"
  moved_staging_dir <- file.path(ctx_initial$staging_root, target_date, run_id)
  moved_progress_dir <- file.path(ctx_initial$log_root, target_date, run_id)
  dir.create(dirname(moved_staging_dir), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(moved_progress_dir), recursive = TRUE, showWarnings = FALSE)
  expect_true(file.rename(ctx_initial$staging_run_dir, moved_staging_dir))
  expect_true(file.rename(ctx_initial$progress_run_dir, moved_progress_dir))

  moved_manifest_path <- file.path(moved_staging_dir, "manifest.rds")
  moved_manifest <- readRDS(moved_manifest_path)
  moved_manifest$run_date <- target_date
  moved_manifest$path_staging_run <- moved_staging_dir
  moved_manifest$path_log_run <- moved_progress_dir
  saveRDS(moved_manifest, moved_manifest_path)
  saveRDS(moved_manifest, file.path(moved_progress_dir, "manifest.rds"))

  ctx_resume <- env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-date-reuse"),
    run_id = run_id,
    sim_grid_chunk_index = 2L,
    sim_grid_n_chunks = 2L
  )

  expect_identical(ctx_resume$run_date, target_date)
  expect_identical(ctx_resume$staging_run_dir, moved_staging_dir)
  expect_identical(ctx_resume$progress_run_dir, moved_progress_dir)
})

test_that("a promoted run cannot be reset to running or lose collation/validation state", {
  env <- .load_runtime_env()

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  ctx <- env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-no-reset"),
    run_id = "promoted-run",
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 1L
  )

  env$.analysis_mark_chunk(
    run_ctx = ctx,
    total_sims = 2L,
    completed_sims = 2L,
    failed_sims = 0L,
    collate_ok = TRUE,
    validation_ok = TRUE
  )
  env$.analysis_promote_run(ctx)

  status_after_promote <- env$.analysis_read_status(ctx)
  expect_identical(status_after_promote$status, "completed")
  expect_true(isTRUE(status_after_promote$promotion_done))

  chunk_list_before <- env$.analysis_read_chunk_statuses(ctx)
  cs_before <- chunk_list_before[[ctx$chunk_label]]
  expect_true(isTRUE(cs_before$collate_ok))
  expect_true(isTRUE(cs_before$validation_ok))

  env$.analysis_mark_chunk(
    run_ctx = ctx,
    total_sims = 2L,
    completed_sims = 2L,
    failed_sims = 0L
  )

  chunk_list_after <- env$.analysis_read_chunk_statuses(ctx)
  cs_after <- chunk_list_after[[ctx$chunk_label]]
  expect_true(isTRUE(cs_after$collate_ok))
  expect_true(isTRUE(cs_after$validation_ok))

  status_after_second <- env$.analysis_read_status(ctx)
  expect_identical(status_after_second$status, "completed")
  expect_true(isTRUE(status_after_second$promotion_done))
})

test_that("reconcile_resume_markers corrects error-path markers from durable output", {
  env <- .load_runtime_env()

  tmp_dir <- withr::local_tempdir()
  file_completed <- file.path(tmp_dir, "completed-2")
  file_error <- file.path(tmp_dir, "error-2")
  file_running <- file.path(tmp_dir, "running-2")
  file.create(file_running)
  file.create(file_completed)

  err_output <- tibble::tibble(sim_id = 2L, error_message = "something went wrong")
  env$.analysis_reconcile_resume_markers(
    existing_output = err_output,
    file_completed = file_completed,
    file_error = file_error,
    file_running = file_running,
    error_col = "error_message"
  )

  expect_false(file.exists(file_completed))
  expect_true(file.exists(file_error))
  expect_false(file.exists(file_running))
})

test_that("a live lock still excludes a second writer", {
  env <- .load_runtime_env()

  tmp_dir <- withr::local_tempdir()
  lock_path <- file.path(tmp_dir, "test.lock")

  acquired <- env$.analysis_acquire_lock(lock_path, timeout_sec = 0.2, poll_sec = 0.05)
  expect_true(acquired)
  expect_true(dir.exists(lock_path))

  second <- env$.analysis_acquire_lock(lock_path, timeout_sec = 0.2, poll_sec = 0.05)
  expect_false(second)

  env$.analysis_release_lock(lock_path)
})

test_that("a stale lock from a dead process can be recovered", {
  env <- .load_runtime_env()

  tmp_dir <- withr::local_tempdir()
  lock_path <- file.path(tmp_dir, "stale.lock")
  dir.create(lock_path)

  stale_meta <- list(
    pid = .Machine$integer.max,
    host = Sys.info()[["nodename"]],
    created_at = format(Sys.time() - 3600, "%Y-%m-%d %H:%M:%S")
  )
  saveRDS(stale_meta, file.path(lock_path, "lock-meta.rds"))

  acquired <- env$.analysis_acquire_lock(lock_path, timeout_sec = 5, poll_sec = 0.05)
  expect_true(acquired)
  env$.analysis_release_lock(lock_path)
})

test_that("a stale lock from another host can be recovered based on age threshold", {
  env <- .load_runtime_env()

  tmp_dir <- withr::local_tempdir()
  lock_path <- file.path(tmp_dir, "remote-stale.lock")
  dir.create(lock_path)

  stale_meta <- list(
    pid = 12345L,
    host = "slurm-remote-node-99",
    token = "remote-node-99-12345-stale",
    created_at = Sys.time() - 3600
  )
  saveRDS(stale_meta, file.path(lock_path, "lock-meta.rds"))

  acquired <- env$.analysis_acquire_lock(lock_path, timeout_sec = 5, poll_sec = 0.05, stale_sec = 10)
  expect_true(acquired)
  expect_true(dir.exists(lock_path))

  new_meta <- readRDS(file.path(lock_path, "lock-meta.rds"))
  expect_identical(new_meta$host, Sys.info()[["nodename"]])
  expect_identical(new_meta$pid, Sys.getpid())

  env$.analysis_release_lock(lock_path)
  expect_false(dir.exists(lock_path))
})

test_that("an orphan lock without metadata can be recovered after orphan timeout", {
  env <- .load_runtime_env()

  tmp_dir <- withr::local_tempdir()
  lock_path <- file.path(tmp_dir, "orphan.lock")
  dir.create(lock_path)

  # Directory exists but no lock-meta.rds
  acquired <- env$.analysis_acquire_lock(lock_path, timeout_sec = 5, poll_sec = 0.05, orphan_stale_sec = 0.05)
  expect_true(acquired)
  expect_true(dir.exists(lock_path))

  new_meta <- readRDS(file.path(lock_path, "lock-meta.rds"))
  expect_identical(new_meta$host, Sys.info()[["nodename"]])
  expect_identical(new_meta$pid, Sys.getpid())

  env$.analysis_release_lock(lock_path)
  expect_false(dir.exists(lock_path))
})

test_that("concurrent stale-lock recovery remains safely serialised without races", {
  env <- .load_runtime_env()

  tmp_dir <- tempfile("stale-concurrent-")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  lock_path <- file.path(tmp_dir, "shared.lock")
  dir.create(lock_path)

  # Plant an initial stale lock
  stale_meta <- list(
    pid = .Machine$integer.max,
    host = Sys.info()[["nodename"]],
    token = "initial-stale",
    created_at = Sys.time() - 3600
  )
  saveRDS(stale_meta, file.path(lock_path, "lock-meta.rds"))

  script_path <- script_runtime
  log_file <- file.path(tmp_dir, "execution.log")
  sentinel_file <- file.path(tmp_dir, "in_critical_section.txt")

  cl <- parallel::makePSOCKcluster(4L)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(
    cl,
    varlist = c("tmp_dir", "lock_path", "log_file", "sentinel_file", "script_path"),
    envir = environment()
  )

  results <- parallel::clusterApply(cl, 1:4, function(worker_id) {
    local_env <- new.env(parent = baseenv())
    source(script_path, local = local_env)

    acquired <- local_env$.analysis_acquire_lock(
      lock_path,
      timeout_sec = 15,
      poll_sec = 0.05
    )

    if (!isTRUE(acquired)) {
      return(list(worker_id = worker_id, ok = FALSE, error = "failed to acquire lock"))
    }

    # Verify mutual exclusion: sentinel file must NOT exist when entering CS
    race_detected <- file.exists(sentinel_file)
    if (!race_detected) {
      writeLines(as.character(worker_id), sentinel_file)
    }

    # Append to log
    cat(paste0("worker-", worker_id, "\n"), file = log_file, append = TRUE)
    Sys.sleep(0.05)

    if (file.exists(sentinel_file)) {
      unlink(sentinel_file, force = TRUE)
    }

    local_env$.analysis_release_lock(lock_path)

    list(worker_id = worker_id, ok = TRUE, race_detected = race_detected)
  })

  for (res in results) {
    expect_true(isTRUE(res$ok))
    expect_false(isTRUE(res$race_detected))
  }

  log_lines <- readLines(log_file)
  expect_length(log_lines, 4L)
})
