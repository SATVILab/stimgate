root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
script_runtime <- file.path(root_dir, "scripts", "r", "analysis-runtime.R")

.load_runtime_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_runtime, local = env)
  env
}

# Forward-slash absolute paths: safe to compare across separator styles and to
# embed in R code run by `Rscript -e` (Windows backslashes are escapes there).
.norm_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

test_that("QMD param lookup follows param > default precedence and env override precedence", {
  env <- .load_runtime_env()
  env$params <- list(sim_grid_chunk_index = 7L)

  expect_identical(env$.get_qmd_param("sim_grid_chunk_index", 1L), 7L)
  expect_identical(env$.get_qmd_param("missing_param", "fallback"), "fallback")

  old_chunk_index <- Sys.getenv("SIM_GRID_CHUNK_INDEX", unset = NA_character_)
  old_n_chunks <- Sys.getenv("SIM_GRID_N_CHUNKS", unset = NA_character_)
  on.exit(
    {
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
    },
    add = TRUE
  )
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

  path <- tempfile("analysis-runtime-", fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)

  obj <- list(
    x = c(1, 2, 3),
    y = data.frame(a = 1:2, b = c("p", "q"))
  )

  expect_identical(env$.write_rds_atomic(obj, path), path)
  expect_equal(readRDS(path), obj)
})

test_that("atomic RDS fallback fails explicitly and retains the pending output", {
  env <- .load_runtime_env()
  folder <- withr::local_tempdir()
  path <- file.path(folder, "result.rds")
  original <- list(value = "last good result")
  replacement <- list(value = "new result")
  saveRDS(original, path)
  env$file.rename <- function(...) FALSE

  expect_identical(env$.write_rds_atomic(replacement, path), path)
  expect_equal(readRDS(path), replacement)
  expect_length(list.files(folder), 1L)

  env$file.copy <- function(...) FALSE
  expect_error(env$.write_rds_atomic(original, path), "Temporary output retained")
  expect_equal(readRDS(path), replacement)
  pending <- list.files(folder, pattern = "\\.tmp-", full.names = TRUE)
  expect_length(pending, 1L)
  expect_equal(readRDS(pending), original)
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
  expect_identical(
    .norm_path(ctx_a$progress_run_dir),
    .norm_path(file.path(
      ctx_a$sim_root, "runs",
      ctx_a$run_date, "run-a"
    ))
  )
  expect_identical(
    readRDS(ctx_a$manifest_path)$path_log_run,
    ctx_a$progress_run_dir
  )
  expect_false(dir.exists(file.path(tmp_project, "cache", "log")))
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
  expect_true(dir.exists(ctx_2$progress_run_dir))
  expect_true(file.exists(ctx_2$status_path))
  expect_false(dir.exists(file.path(ctx_2$current_dir, "runs")))
  expect_false(file.exists(file.path(ctx_2$current_dir, "status.rds")))
  expect_false(dir.exists(file.path(ctx_2$current_dir, "jobs")))

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

test_that("explicit run ID reuse checks scientific but not operational parameters", {
  env <- .load_runtime_env()

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-scientific-params"),
    run_id = "shared-scientific-run",
    params = list(
      run_simulations = TRUE,
      run_plots = FALSE,
      sim_grid_chunk_index = 1L,
      sim_grid_n_chunks = 2L,
      comparison_semantics_version = "corrected-v1"
    ),
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 2L
  )

  expect_no_error(env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-scientific-params"),
    run_id = "shared-scientific-run",
    params = list(
      run_simulations = FALSE,
      run_plots = TRUE,
      sim_grid_chunk_index = 2L,
      sim_grid_n_chunks = 2L,
      comparison_semantics_version = "corrected-v1"
    ),
    sim_grid_chunk_index = 2L,
    sim_grid_n_chunks = 2L
  ))

  expect_error(
    env$.analysis_run_context(
      analysis_key = c("sim", "analysis-runtime-scientific-params"),
      run_id = "shared-scientific-run",
      params = list(
        run_simulations = TRUE,
        run_plots = FALSE,
        sim_grid_chunk_index = 2L,
        sim_grid_n_chunks = 2L,
        comparison_semantics_version = "corrected-v2"
      ),
      sim_grid_chunk_index = 2L,
      sim_grid_n_chunks = 2L
    ),
    "incompatible with existing manifest"
  )
})

test_that("canonical current reads require completion and compatible provenance", {
  env <- .load_runtime_env()

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  ctx <- env$.analysis_run_context(
    analysis_key = c("sim", "analysis-runtime-current"),
    run_id = "current-run",
    params = list(comparison_semantics_version = "corrected-v1")
  )
  path_staged <- file.path(ctx$staging_collated_dir, "result.rds")
  env$.write_rds_atomic(tibble::tibble(x = 1L), path_staged)

  expect_error(
    env$.analysis_current_file(
      ctx,
      c("collated", "result.rds"),
      list(comparison_semantics_version = "corrected-v1")
    ),
    "No complete canonical current result"
  )

  env$.analysis_mark_chunk(
    run_ctx = ctx,
    total_sims = 1L,
    completed_sims = 1L,
    failed_sims = 0L,
    collate_ok = TRUE,
    validation_ok = TRUE
  )
  expect_true(isTRUE(env$.analysis_promote_run(ctx)))

  expect_identical(
    env$.analysis_current_file(
      ctx,
      c("collated", "result.rds"),
      list(comparison_semantics_version = "corrected-v1")
    ),
    file.path(ctx$current_dir, "collated", "result.rds")
  )
  expect_error(
    env$.analysis_current_file(
      ctx,
      c("collated", "result.rds"),
      list(comparison_semantics_version = "legacy-v0")
    ),
    "comparison_semantics_version"
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
  # Workers need the project library, where projr resolves their paths.
  parallel::clusterCall(cl, function(paths) base::.libPaths(paths), .libPaths())

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
  # Workers need the project library, where projr resolves their paths.
  parallel::clusterCall(cl, function(paths) base::.libPaths(paths), .libPaths())
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
  moved_progress_dir <- file.path(ctx_initial$runs_root, target_date, run_id)
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
  expect_identical(
    .norm_path(ctx_resume$staging_run_dir),
    .norm_path(moved_staging_dir)
  )
  expect_identical(
    .norm_path(ctx_resume$progress_run_dir),
    .norm_path(moved_progress_dir)
  )
})

test_that("legacy run state resumes at the path recorded in the staged manifest", {
  env <- .load_runtime_env()
  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  analysis_key <- c("sim", "analysis-runtime-legacy-resume")
  ctx <- env$.analysis_run_context(
    analysis_key = analysis_key,
    run_id = "legacy-run",
    sim_grid_n_chunks = 2L
  )
  env$.analysis_mark_chunk(
    run_ctx = ctx,
    total_sims = 1L,
    completed_sims = 1L,
    collate_ok = TRUE,
    validation_ok = TRUE
  )
  writeLines("existing progress", ctx$progress_file)
  file.create(file.path(ctx$chunk_jobs_dir, "completed-1"))
  status_before <- readRDS(ctx$status_path)

  # Recreate the old runtime layout, including its leading "sim" key component.
  legacy_dir <- file.path(
    tmp_project, "cache", "log", "analysis", "sim",
    "analysis-runtime-legacy-resume", ctx$run_date, ctx$run_id
  )
  dir.create(dirname(legacy_dir), recursive = TRUE, showWarnings = FALSE)
  expect_true(file.rename(ctx$progress_run_dir, legacy_dir))
  manifest <- readRDS(ctx$manifest_path)
  manifest$path_log_run <- legacy_dir
  saveRDS(manifest, ctx$manifest_path)
  saveRDS(manifest, file.path(legacy_dir, "manifest.rds"))

  resumed <- env$.analysis_run_context(
    analysis_key = analysis_key,
    run_id = ctx$run_id,
    sim_grid_chunk_index = 2L,
    sim_grid_n_chunks = 2L
  )
  expect_identical(.norm_path(resumed$staging_run_dir), .norm_path(ctx$staging_run_dir))
  expect_identical(.norm_path(resumed$progress_run_dir), .norm_path(legacy_dir))
  expect_identical(readRDS(resumed$status_path), status_before)
  expect_identical(readLines(resumed$progress_file), "existing progress")
  expect_true(file.exists(file.path(
    legacy_dir, "jobs", ctx$chunk_label, "completed-1"
  )))
  expect_true(dir.exists(file.path(legacy_dir, "jobs", resumed$chunk_label)))
  expect_true(ctx$chunk_label %in% names(env$.analysis_read_chunk_statuses(resumed)))
  expect_identical(
    .norm_path(env$.analysis_lock_path(resumed, "promotion")),
    .norm_path(file.path(dirname(resumed$current_dir), "promotion.lock"))
  )
  expect_false(dir.exists(ctx$progress_run_dir))
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

test_that("one process holding the lock excludes another process, and unlock allows acquisition", {
  env <- .load_runtime_env()

  tmp_dir <- withr::local_tempdir()
  lock_path <- .norm_path(file.path(tmp_dir, "test.lock"))
  script_path <- .norm_path(script_runtime)

  lock1 <- env$.analysis_acquire_lock(lock_path, timeout_sec = 0.5)
  expect_false(is.null(lock1))
  withr::defer(env$.analysis_release_lock(lock1))

  # Another process attempting to acquire the same lock times out and returns NULL
  cmd_excluded <- sprintf(
    "local_env <- new.env(parent = baseenv()); source(\"%s\", local = local_env); l <- local_env$.analysis_acquire_lock(\"%s\", timeout_sec = 0.1); cat(if (is.null(l)) \"EXCLUDED\" else \"ACQUIRED\")",
    script_path, lock_path
  )
  out_excluded <- system2("Rscript", args = c("-e", shQuote(cmd_excluded)), stdout = TRUE)
  expect_true(any(grepl("EXCLUDED", out_excluded)))

  # After normal unlock, another process can acquire it
  env$.analysis_release_lock(lock1)

  cmd_acquired <- sprintf(
    "local_env <- new.env(parent = baseenv()); source(\"%s\", local = local_env); l <- local_env$.analysis_acquire_lock(\"%s\", timeout_sec = 1); if (!is.null(l)) { local_env$.analysis_release_lock(l); cat(\"ACQUIRED_OK\") }",
    script_path, lock_path
  )
  out_acquired <- system2("Rscript", args = c("-e", shQuote(cmd_acquired)), stdout = TRUE)
  expect_true(any(grepl("ACQUIRED_OK", out_acquired)))
})

test_that("a worker process terminating without unlocking leaves lock immediately acquirable", {
  # On Windows, tools::pskill() cannot terminate the background Rscript worker
  # (it reports failure and the worker keeps the lock), so this Unix
  # kill-and-release check cannot run there.
  skip_on_os("windows")
  env <- .load_runtime_env()

  tmp_dir <- withr::local_tempdir()
  lock_path <- .norm_path(file.path(tmp_dir, "termination.lock"))
  ready_file <- .norm_path(file.path(tmp_dir, "ready.txt"))
  pid_file <- .norm_path(file.path(tmp_dir, "pid.txt"))
  script_path <- .norm_path(script_runtime)

  # Launch worker in background that acquires lock and signals readiness
  worker_cmd <- sprintf(
    "local_env <- new.env(parent = baseenv()); source(\"%s\", local = local_env); writeLines(as.character(Sys.getpid()), \"%s\"); l <- local_env$.analysis_acquire_lock(\"%s\", timeout_sec = 10); writeLines(\"READY\", \"%s\"); Sys.sleep(100)",
    script_path, pid_file, lock_path, ready_file
  )
  system2("Rscript", args = c("-e", shQuote(worker_cmd)), wait = FALSE)

  # Wait until worker has locked and signaled readiness
  for (i in seq_len(300L)) {
    if (file.exists(ready_file)) break
    Sys.sleep(0.1)
  }
  expect_true(file.exists(ready_file))

  # Parent process cannot acquire while worker is alive
  l_while_alive <- env$.analysis_acquire_lock(lock_path, timeout_sec = 0.1)
  expect_null(l_while_alive)

  # Terminate worker process abruptly (SIGKILL = 9)
  worker_pid <- as.integer(readLines(pid_file)[[1]])
  expect_true(tools::pskill(worker_pid, signal = tools::SIGKILL))
  Sys.sleep(0.1)

  # The OS drops the lock when the holder dies (fcntl on Unix; on Windows the
  # release can lag slightly behind TerminateProcess), so it becomes acquirable.
  l_after_death <- env$.analysis_acquire_lock(lock_path, timeout_sec = 10)
  expect_false(is.null(l_after_death))
  env$.analysis_release_lock(l_after_death)
})

test_that("concurrent lock acquisition across workers remains safely serialised", {
  env <- .load_runtime_env()

  tmp_dir <- tempfile("flock-concurrent-")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  lock_path <- file.path(tmp_dir, "shared.lock")
  script_path <- script_runtime
  log_file <- file.path(tmp_dir, "execution.log")

  cl <- parallel::makePSOCKcluster(4L)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(
    cl,
    varlist = c("tmp_dir", "lock_path", "log_file", "script_path"),
    envir = environment()
  )

  results <- parallel::clusterApply(cl, 1:4, function(worker_id) {
    local_env <- new.env(parent = baseenv())
    source(script_path, local = local_env)

    lock <- local_env$.analysis_acquire_lock(
      lock_path,
      timeout_sec = 15
    )

    if (is.null(lock)) {
      return(list(worker_id = worker_id, ok = FALSE, error = "failed to acquire lock"))
    }

    # Append to log
    cat(paste0("worker-", worker_id, "\n"), file = log_file, append = TRUE)
    Sys.sleep(0.05)

    local_env$.analysis_release_lock(lock)

    list(worker_id = worker_id, ok = TRUE)
  })

  for (res in results) {
    expect_true(isTRUE(res$ok))
  }

  log_lines <- readLines(log_file)
  expect_length(log_lines, 4L)
})

test_that("results context reads promoted outputs without run state", {
  env <- .load_runtime_env()

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")
  key <- c("sim", "analysis-runtime-results-context")

  expect_error(
    env$.analysis_results_context(key, path_root = tmp_project),
    "No complete canonical current result"
  )
  expect_false(dir.exists(file.path(tmp_project, "cache", "log")))

  ctx <- env$.analysis_run_context(
    analysis_key = key,
    run_id = "promoted-run",
    path_root = tmp_project
  )
  env$.write_rds_atomic(
    data.frame(x = 1),
    file.path(ctx$staging_collated_dir, "result.rds")
  )
  env$.analysis_mark_chunk(
    run_ctx = ctx,
    total_sims = 1L,
    completed_sims = 1L,
    failed_sims = 0L,
    collate_ok = TRUE,
    validation_ok = TRUE
  )
  expect_true(isTRUE(env$.analysis_promote_run(ctx)))

  staging_before <- list.files(ctx$staging_root, recursive = TRUE)
  ctx_read <- env$.analysis_results_context(key, path_root = tmp_project)

  expect_true(ctx_read$read_only)
  expect_identical(
    normalizePath(ctx_read$current_dir),
    normalizePath(ctx$current_dir)
  )
  expect_identical(ctx_read$staging_run_dir, ctx_read$current_dir)
  expect_identical(
    readRDS(file.path(ctx_read$staging_collated_dir, "result.rds")),
    data.frame(x = 1)
  )
  expect_identical(
    env$.analysis_current_file(ctx_read, c("collated", "result.rds")),
    file.path(ctx_read$current_dir, "collated", "result.rds")
  )
  expect_identical(
    list.files(ctx$staging_root, recursive = TRUE),
    staging_before
  )
})

test_that("project paths resolve through projr from path_root and can be read-only", {
  env <- .load_runtime_env()
  project <- .local_projr_root()
  expected <- .projr_output_path(project, "fig")
  expect_identical(
    .norm_path(env$.analysis_project_dir(
      "output", "fig", project, create = FALSE
    )),
    .norm_path(expected)
  )
  expect_false(dir.exists(expected))
  expect_identical(
    .norm_path(env$.analysis_project_dir("output", "fig", project)),
    .norm_path(expected)
  )
  expect_true(dir.exists(expected))
  # A file path creates only its parent folder.
  file <- env$.analysis_project_dir("output", c("table", "a.csv"), project, dir = FALSE)
  expect_identical(.norm_path(file), .norm_path(.projr_output_path(project, "table", "a.csv")))
  expect_true(dir.exists(dirname(file)))
  expect_false(file.exists(file))
  cache <- env$.analysis_cache_dir(c("sim", "test"), project, create = FALSE)
  expect_identical(.norm_path(cache), .norm_path(file.path(
    project, "_tmp", "sim", "test"
  )))
  expect_false(dir.exists(cache))
  # Output goes to projr's cache, never a committable output/ folder.
  expect_false(dir.exists(file.path(project, "output")))
  expect_false(dir.exists(file.path(project, "_output")))
})

test_that("figure directories are nested under the output fig folder and created", {
  env <- .load_runtime_env()
  project <- .local_projr_root()
  expected <- .projr_output_path(project, "fig", "2a-x", "quick", "signed_error")
  fig_dir <- env$.analysis_fig_dir(
    c("2a-x", "quick", "signed_error"),
    path_root = project
  )
  expect_identical(.norm_path(fig_dir), .norm_path(expected))
  expect_true(dir.exists(expected))
})

test_that("seeded evaluation is independent of and restores caller RNG", {
  env <- .load_runtime_env()
  withr::local_preserve_seed()
  old_kind <- RNGkind()
  withr::defer(RNGkind(old_kind[[1]], old_kind[[2]], old_kind[[3]]))
  RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rejection")
  set.seed(18L)
  before_kind <- RNGkind()
  before_seed <- .Random.seed
  draws <- env$.analysis_with_seed(5L, stats::rnorm(3))
  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before_seed)
  RNGkind("default", "default", "default")
  expect_identical(env$.analysis_with_seed(5L, stats::rnorm(3)), draws)
  expect_error(env$.analysis_with_seed(5L, stop("boom")), "boom")
  rm(".Random.seed", envir = .GlobalEnv)
  env$.analysis_with_seed(5L, stats::runif(1))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("analysis profiles are read from PROJR_PROFILE", {
  env <- new.env(parent = baseenv())
  source(script_runtime, local = env)
  withr::local_envvar(PROJR_PROFILE = NA)
  expect_false(env$.analysis_is_dev())
  expect_false(env$.analysis_is_quick())
  Sys.setenv(PROJR_PROFILE = "dev, quick")
  expect_true(env$.analysis_is_dev())
  expect_true(env$.analysis_is_quick())
  Sys.setenv(PROJR_PROFILE = "default")
  expect_false(env$.analysis_is_dev())
  expect_false(env$.analysis_is_quick())
})

test_that("analysis mode keys isolate dev and quick results with dev precedence", {
  withr::local_dir(root_dir)
  env <- new.env(parent = baseenv())
  source(script_runtime, local = env)
  key <- c("sim", "test")
  # _projr.yml chooses the size; SIM_SIZE makes the test independent of it.
  withr::local_envvar(PROJR_PROFILE = NA, SIM_SIZE = "draft")
  # Draft results get their own key; final runs keep the plain key.
  expect_identical(env$.analysis_mode_key(key), c(key, "draft"))
  Sys.setenv(SIM_SIZE = "final")
  expect_identical(env$.analysis_mode_key(key), key)
  Sys.setenv(PROJR_PROFILE = "quick")
  expect_identical(env$.analysis_mode_key(key), c(key, "quick"))
  Sys.setenv(PROJR_PROFILE = "dev")
  expect_identical(env$.analysis_mode_key(key), c(key, "dev"))
  Sys.setenv(PROJR_PROFILE = "dev, quick")
  expect_identical(env$.analysis_mode_key(key), c(key, "dev"))
})

test_that("canonical cache failures give render guidance and successful reads retain data", {
  env <- .load_runtime_env()
  cache <- withr::local_tempdir()
  ctx <- list(
    analysis_key = c("sim", "test"), current_dir = cache,
    qmd_path = "analysis/test.qmd"
  )
  path <- file.path(cache, "result.rds")
  command <- "RUN_SIMULATIONS=true RUN_PLOTS=false SIM_SIZE=final quarto render analysis/test.qmd"
  expect_error(env$.analysis_read_current(ctx, "result.rds"), command, fixed = TRUE)
  file.create(file.path(cache, "COMPLETE"))
  expect_error(env$.analysis_read_current(ctx, "result.rds"), command, fixed = TRUE)
  saveRDS(list(analysis_key = ctx$analysis_key, params = list(version = 1L)),
          file.path(cache, "manifest.rds"))
  expect_error(env$.analysis_read_current(ctx, "result.rds"), command, fixed = TRUE)
  object <- data.frame(value = 1:2, row.names = c("a", "b"))
  saveRDS(object, path)
  expect_identical(env$.analysis_read_current(ctx, "result.rds", list(version = 1L)), object)
  expect_error(env$.analysis_read_current(ctx, "result.rds", list(version = 2L)), command, fixed = TRUE)
  writeLines("corrupt RDS", path)
  expect_error(suppressWarnings(env$.analysis_read_current(ctx, "result.rds")), command, fixed = TRUE)
})

test_that("the shared plot saver writes and prints the supplied plot once", {
  env <- .load_runtime_env()
  printed <- list()
  env$print <- function(x) printed[[length(printed) + 1L]] <<- x
  directory <- withr::local_tempdir()
  path <- file.path(directory, "figures", "plot.pdf")
  plot <- ggplot2::ggplot(data.frame(x = 1:2, y = 1:2), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  expect_identical(env$.analysis_save_plot(plot, path, width = 5, height = 5), path)
  expect_true(file.exists(path))
  expect_identical(printed, list(plot))
})
