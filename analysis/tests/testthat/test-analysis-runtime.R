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
