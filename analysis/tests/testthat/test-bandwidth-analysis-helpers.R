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

test_that("output discovery checks both active and cache directories", {
  env <- .load_bw_analysis_env()
  tmp_dir <- tempfile("bw-discovery")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  active_dir <- file.path(tmp_dir, "active-output")
  cache_dir <- file.path(tmp_dir, "cache")
  dir.create(active_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cache_dir, "output"), recursive = TRUE, showWarnings = FALSE)

  file_1 <- file.path(active_dir, "bw_list_raw-chunk_001-of_002-sim_id_000001.rds")
  file_2 <- file.path(cache_dir, "output", "bw_list_raw-chunk_001-of_002-sim_id_000002.rds")
  saveRDS(1, file_1)
  saveRDS(2, file_2)

  found <- env$.find_bw_list_output_files(output_dir = active_dir, cache_dir = cache_dir)
  expect_setequal(basename(found), c(
    "bw_list_raw-chunk_001-of_002-sim_id_000001.rds",
    "bw_list_raw-chunk_001-of_002-sim_id_000002.rds"
  ))
})

test_that("actual QMD cache fallbacks are discoverable through the shared helper", {
  env <- .load_bw_analysis_env()
  tmp_dir <- tempfile("bw-qmd-fallbacks")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  qmd_cache_dirs <- list(
    global = file.path(tmp_dir, "cache", "log", "analysis", "sim", "bw", "freq_bs", "global"),
    est_adaptive = file.path(tmp_dir, "cache", "log", "analysis", "sim", "bw", "est", "adaptive"),
    freq_adaptive = file.path(tmp_dir, "cache", "log", "analysis", "sim", "bw", "freq_bs", "adaptive")
  )

  for (cache_dir in qmd_cache_dirs) {
    out_dir <- file.path(cache_dir, "output")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      1L,
      file.path(out_dir, "bw_list_raw-chunk_001-of_002-sim_id_000001.rds")
    )
  }

  expect_equal(
    basename(env$.find_bw_list_output_files(
      output_dir = NULL,
      cache_dir = qmd_cache_dirs[["global"]],
      cache_path = c("cache", "log", "analysis", "sim", "bw", "freq_bs", "global")
    )),
    "bw_list_raw-chunk_001-of_002-sim_id_000001.rds"
  )

  expect_equal(
    basename(env$.find_bw_list_output_files(
      output_dir = NULL,
      cache_dir = qmd_cache_dirs[["est_adaptive"]],
      cache_path = c("cache", "log", "analysis", "sim", "bw", "est", "adaptive")
    )),
    "bw_list_raw-chunk_001-of_002-sim_id_000001.rds"
  )

  expect_equal(
    basename(env$.find_bw_list_output_files(
      output_dir = NULL,
      cache_dir = qmd_cache_dirs[["freq_adaptive"]],
      cache_path = c("cache", "log", "analysis", "sim", "bw", "freq_bs", "adaptive")
    )),
    "bw_list_raw-chunk_001-of_002-sim_id_000001.rds"
  )
})

test_that("QMD 6 does not retain the stale adaptive shim reference", {
  qmd_lines <- readLines(qmd_6, warn = FALSE)
  expect_false(any(grepl("\\.run_sim_bandwidth_bs_freq_adaptive", qmd_lines)))
})
