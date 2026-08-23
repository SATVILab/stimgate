root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_bw_io <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")
qmd_6 <- file.path(root_dir, "analysis", "6-sim-bw-freq_bs-adaptive.qmd")

.load_bw_analysis_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))

  for (f in c(script_bw_io, script_bw_plot)) {
    if (!file.exists(f)) {
      stop("Expected analysis helper not found: ", f)
    }
    source(f, local = env)
  }

  env
}

test_that("legacy output path and formatting helpers match the original contracts", {
  env <- .load_bw_analysis_env()

  tmp_dir <- tempfile("bw-output")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  output_dir <- file.path(tmp_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  expect_equal(
    env$.path_sim_output(
      42L,
      dir_output = output_dir,
      sim_grid_chunk_index = 1L,
      sim_grid_n_chunks = 2L
    ),
    file.path(output_dir, "bw_list_raw-chunk_001-of_002-sim_id_000042.rds")
  )

  bw_vec <- c(0.1, 0.25, NA_real_)
  expect_equal(env$format_bw_lab(bw_vec), c("0.1", "0.25", NA_character_))
  expect_equal(env$format_bw_file(bw_vec), c("0p1", "0p25", NA_character_))
  expect_equal(env$safe_file_lab(c("bw + 1", "bw_2", "low---pass")), c("bw_1", "bw_2", "low_pass"))

  out <- env$add_bw_labs(data.frame(bw_core = c(0.1, 0.25), bw_extra = c(0.2, 0.3)))
  expect_true("bw_core_lab" %in% names(out))
  expect_equal(out$bw_core_lab, c("0.1", "0.25"))
  expect_equal(out$bw_extra_lab, c("0.2", "0.3"))
})

test_that("make_bw_colour_values preserves old ordering and palette semantics", {
  env <- .load_bw_analysis_env()
  base_col_vec <- c("#111111", "#222222", "#333333", "#444444", "#555555", "#666666")
  bw_vec <- c(0.5, 1, 0.5, NA_real_, 2)

  bw_vals <- env$make_bw_colour_values(bw_vec, base_col_vec = base_col_vec)
  expect_equal(names(bw_vals), c("0.5", "1", "2"))
  expect_equal(unname(bw_vals), c("#444444", "#555555", "#666666"))
})

test_that("make_bw_linetype_scale keeps the legacy list contract", {
  env <- .load_bw_analysis_env()

  lty <- env$make_bw_linetype_scale(c(1, 0.5, NA_real_))
  expect_equal(lty$levels, c("0.5", "1"))
  expect_equal(lty$values, c("0.5" = "solid", "1" = "84"))
})

test_that("output discovery prefers active output and does not mix with cache", {
  env <- .load_bw_analysis_env()
  tmp_dir <- tempfile("bw-discovery")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  active_dir <- file.path(tmp_dir, "active-output")
  cache_dir <- file.path(tmp_dir, "cache", "sim", "bw", "freq_bs", "adaptive", "current")
  dir.create(active_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cache_dir, "chunks", "001-of-002", "output"), recursive = TRUE, showWarnings = FALSE)

  file_1 <- file.path(active_dir, "bw_list_raw-chunk_001-of_002-sim_id_000001.rds")
  file_2 <- file.path(
    cache_dir,
    "chunks",
    "001-of-002",
    "output",
    "bw_list_raw-chunk_001-of_002-sim_id_000002.rds"
  )
  saveRDS(1, file_1)
  saveRDS(2, file_2)

  found <- env$.find_bw_list_output_files(output_dir = active_dir, cache_dir = cache_dir)
  expect_equal(basename(found), "bw_list_raw-chunk_001-of_002-sim_id_000001.rds")
})

test_that("output discovery falls back to cache only when active output is empty", {
  env <- .load_bw_analysis_env()
  tmp_dir <- tempfile("bw-discovery-fallback")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  active_dir <- file.path(tmp_dir, "active-output")
  cache_dir <- file.path(tmp_dir, "cache", "sim", "bw", "freq_bs", "adaptive", "current")
  dir.create(active_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cache_dir, "chunks", "001-of-002", "output"), recursive = TRUE, showWarnings = FALSE)

  file_2 <- file.path(
    cache_dir,
    "chunks",
    "001-of-002",
    "output",
    "bw_list_raw-chunk_001-of_002-sim_id_000002.rds"
  )
  saveRDS(2, file_2)

  found <- env$.find_bw_list_output_files(output_dir = active_dir, cache_dir = cache_dir)
  expect_equal(basename(found), "bw_list_raw-chunk_001-of_002-sim_id_000002.rds")
})

test_that("output discovery can disable cache fallback", {
  env <- .load_bw_analysis_env()
  tmp_dir <- tempfile("bw-discovery-no-cache")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  active_dir <- file.path(tmp_dir, "active-output")
  cache_dir <- file.path(tmp_dir, "cache", "sim", "bw", "freq_bs", "adaptive", "current")
  dir.create(active_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cache_dir, "chunks", "001-of-002", "output"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(2L, file.path(
    cache_dir,
    "chunks",
    "001-of-002",
    "output",
    "bw_list_raw-chunk_001-of_002-sim_id_000002.rds"
  ))

  found <- env$.find_bw_list_output_files(
    output_dir = active_dir,
    cache_dir = cache_dir,
    allow_cache_fallback = FALSE
  )
  expect_identical(found, character(0))
})

test_that("canonical current results are discoverable for run_simulations = FALSE", {
  env <- .load_bw_analysis_env()
  tmp_dir <- tempfile("bw-current-discovery")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  staging_run_dir <- file.path(tmp_dir, "staging", "2026-08-22", "run-123")
  current_dir <- file.path(tmp_dir, "current")
  dir.create(staging_run_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(current_dir, "chunks", "001-of-002", "output"), recursive = TRUE, showWarnings = FALSE)

  file_curr <- file.path(
    current_dir,
    "chunks",
    "001-of-002",
    "output",
    "bw_list_raw-chunk_001-of_002-sim_id_000001.rds"
  )
  saveRDS(1L, file_curr)

  found <- env$.find_bw_list_output_files(
    output_dir = staging_run_dir,
    cache_dir = current_dir
  )
  expect_equal(basename(found), "bw_list_raw-chunk_001-of_002-sim_id_000001.rds")
})

test_that("output discovery does not use implicit legacy fallbacks", {
  env <- .load_bw_analysis_env()
  tmp_dir <- tempfile("bw-no-fallback")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(file.path(tmp_dir, "output"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(1L, file.path(tmp_dir, "output", "bw_list_raw-chunk_001-of_001-sim_id_000001.rds"))

  expect_identical(
    env$.find_bw_list_output_files(),
    character(0)
  )
})

test_that("QMD 6 does not retain the stale adaptive shim reference", {
  qmd_lines <- readLines(qmd_6, warn = FALSE)
  expect_false(any(grepl("\\.run_sim_bandwidth_bs_freq_adaptive", qmd_lines)))
})

test_that(".update_progress_summary() works without chunk/output metadata for non-chunked analyses", {
  env <- .load_bw_analysis_env()

  tmp_dir <- tempfile("bw-progress")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  dir_jobs <- file.path(tmp_dir, "jobs")
  dir.create(dir_jobs, recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(dir_jobs, "running-1"))
  file.create(file.path(dir_jobs, "completed-2"))
  path_progress_file <- file.path(tmp_dir, "progress.txt")

  summary_text <- env$.update_progress_summary(
    path_progress_file = path_progress_file,
    dir_jobs_chunk = dir_jobs,
    total_sims = 2L,
    heading = "SIMULATION PROGRESS DASHBOARD"
  )

  expect_true(file.exists(path_progress_file))
  expect_false(grepl("Chunk", summary_text))
  expect_false(grepl("Output directory", summary_text))
  expect_true(grepl("SIMULATION PROGRESS DASHBOARD", summary_text))
  expect_true(grepl("Total Simulations  : 2", summary_text, fixed = TRUE))
})

test_that(".update_progress_summary() still reports chunk/output metadata when supplied", {
  env <- .load_bw_analysis_env()

  tmp_dir <- tempfile("bw-progress")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  dir_jobs <- file.path(tmp_dir, "jobs")
  dir.create(dir_jobs, recursive = TRUE, showWarnings = FALSE)
  path_progress_file <- file.path(tmp_dir, "progress.txt")

  summary_text <- env$.update_progress_summary(
    path_progress_file = path_progress_file,
    dir_jobs_chunk = dir_jobs,
    total_sims = 5L,
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 3L,
    dir_output = tmp_dir
  )

  expect_true(grepl("Chunk              : 1 of 3", summary_text, fixed = TRUE))
  expect_true(grepl(paste0("Output directory   : ", tmp_dir), summary_text, fixed = TRUE))
  expect_true(grepl("ADAPTIVE SIMULATION PROGRESS DASHBOARD", summary_text))
})

test_that(".update_progress_summary() tolerates a concurrent write failure", {
  env <- .load_bw_analysis_env()

  tmp_dir <- tempfile("bw-progress")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  dir_jobs <- file.path(tmp_dir, "jobs")
  dir.create(dir_jobs, recursive = TRUE, showWarnings = FALSE)
  # A path under a non-existent directory forces the write to fail.
  path_progress_file <- file.path(tmp_dir, "missing-parent", "progress.txt")

  expect_no_error(
    suppressWarnings(
      env$.update_progress_summary(
        path_progress_file = path_progress_file,
        dir_jobs_chunk = dir_jobs,
        total_sims = 1L
      )
    )
  )
  expect_false(file.exists(path_progress_file))
})

test_that("QMDs 3 and 4 define chunk controls and use the run-scoped output dir", {
  qmd_3 <- readLines(file.path(root_dir, "analysis", "3-sim-bw-est-base.qmd"), warn = FALSE)
  qmd_4 <- readLines(file.path(root_dir, "analysis", "4-sim-bw-est-norm.qmd"), warn = FALSE)

  for (qmd_lines in list(qmd_3, qmd_4)) {
    expect_true(any(grepl(
      "sim_grid_chunk_index <- as.integer",
      qmd_lines,
      fixed = TRUE
    )))
    expect_true(any(grepl(
      "sim_grid_n_chunks <- as.integer",
      qmd_lines,
      fixed = TRUE
    )))
    expect_true(any(grepl(
      "dir_output <- run_ctx$chunk_output_dir",
      qmd_lines,
      fixed = TRUE
    )))
  }
})
