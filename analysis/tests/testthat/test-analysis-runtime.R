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
