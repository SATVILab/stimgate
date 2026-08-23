root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

read_project_file <- function(...) {
  paste(
    readLines(file.path(root_dir, ...), warn = FALSE),
    collapse = "\n"
  )
}

test_that("analyses 3 and 4 use resumable transactional run contexts", {
  qmd_specs <- list(
    list(file = "3-sim-bw-est-base.qmd", key = "base"),
    list(file = "4-sim-bw-est-norm.qmd", key = "norm")
  )

  for (spec in qmd_specs) {
    content <- read_project_file("analysis", spec$file)

    expect_true(grepl("analysis_run_id:\\s*null", content))
    expect_true(grepl("sim_grid_chunk_index:\\s*1", content))
    expect_true(grepl("sim_grid_n_chunks:\\s*1", content))
    expect_true(grepl("sim_grid_shuffle_seed", content, fixed = TRUE))
    expect_true(grepl(
      paste0('c\\("sim",\\s*"bw",\\s*"est",\\s*"', spec$key, '"\\)'),
      content
    ))
    expect_true(grepl("run_ctx$progress_file", content, fixed = TRUE))
    expect_true(grepl("run_ctx$chunk_output_dir", content, fixed = TRUE))
    expect_true(grepl(".path_sim_output", content, fixed = TRUE))
    expect_true(grepl(".analysis_reconcile_resume_markers", content, fixed = TRUE))
    expect_true(grepl(".analysis_mark_chunk", content, fixed = TRUE))
    expect_true(grepl(".analysis_can_promote", content, fixed = TRUE))
    expect_true(grepl(".analysis_promote_run", content, fixed = TRUE))
    expect_true(grepl("expected_sim_ids", content, fixed = TRUE))
    expect_true(grepl("run_ctx$current_dir", content, fixed = TRUE))
  }
})

test_that("analysis Slurm launchers render top-level QMDs with one run ID", {
  launcher_specs <- list(
    list(script = "dev-2-stim-bw-freq_bs-global.sh", qmd = "2-sim-bw-freq_bs-global.qmd"),
    list(script = "dev-3-sim-bw-est-base.sh", qmd = "3-sim-bw-est-base.qmd"),
    list(script = "dev-4-sim-bw-est-norm.sh", qmd = "4-sim-bw-est-norm.qmd"),
    list(script = "dev-5-sim-bw-est-adaptive.sh", qmd = "5-sim-bw-est-adaptive.qmd"),
    list(script = "dev-6-sim-bw-freq_bs-adaptive.sh", qmd = "6-sim-bw-freq_bs-adaptive.qmd")
  )

  for (spec in launcher_specs) {
    content <- read_project_file("scripts", "slurm", spec$script)

    expect_true(grepl(
      paste0("analysis/", spec$qmd),
      content,
      fixed = TRUE
    ))
    expect_false(grepl("analysis/split/", content, fixed = TRUE))
    expect_true(grepl("export ANALYSIS_RUN_ID", content, fixed = TRUE))
    expect_true(grepl("SIM_GRID_CHUNK_INDEX", content, fixed = TRUE))
    expect_true(grepl("SIM_GRID_N_CHUNKS", content, fixed = TRUE))
    expect_true(grepl("RUN_SIMULATIONS", content, fixed = TRUE))
    expect_true(grepl("RUN_PLOTS", content, fixed = TRUE))
  }

  dev_launcher <- read_project_file("scripts", "slurm", "dev.sh")
  expect_false(grepl("prepare_split_qmds", dev_launcher, fixed = TRUE))
  expect_true(grepl("ANALYSIS_RUN_ID", dev_launcher, fixed = TRUE))
})

test_that("obsolete analysis 2 and 6 split QMDs are removed", {
  split_paths <- c(
    file.path(
      root_dir,
      "analysis",
      "split",
      paste0("2-sim-bw-freq_bs-global-", 1:4, ".qmd")
    ),
    file.path(
      root_dir,
      "analysis",
      "split",
      paste0("6-sim-bw-freq_bs-adaptive-", 1:4, ".qmd")
    )
  )

  expect_false(any(file.exists(split_paths)))
})

test_that("promoted bandwidth outputs are discoverable without a simulation run", {
  runtime_env <- new.env(parent = getNamespace("stimgate"))
  source(
    file.path(root_dir, "scripts", "r", "analysis-runtime.R"),
    local = runtime_env
  )
  source(
    file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R"),
    local = runtime_env
  )

  tmp_project <- withr::local_tempdir()
  withr::local_dir(tmp_project)
  writeLines(c("directories:", "  docs:", "    path: docs"), "_projr.yml")

  run_ctx <- runtime_env$.analysis_run_context(
    analysis_key = c("sim", "bw", "est", "base"),
    run_id = "analysis-339-current",
    params = list(
      sim_grid_n_chunks = 1L,
      sim_grid_shuffle_seed = 20260707L
    )
  )
  path_output <- runtime_env$.path_sim_output(
    sim_id = 1L,
    dir_output = run_ctx$chunk_output_dir,
    sim_grid_chunk_index = 1L,
    sim_grid_n_chunks = 1L
  )
  runtime_env$.write_rds_atomic(
    tibble::tibble(sim_id = 1L, bw = 0.2, error = NA_character_),
    path_output
  )
  runtime_env$.analysis_mark_chunk(
    run_ctx = run_ctx,
    total_sims = 1L,
    completed_sims = 1L,
    failed_sims = 0L,
    collate_ok = TRUE,
    validation_ok = TRUE
  )
  expect_true(isTRUE(runtime_env$.analysis_promote_run(run_ctx)))

  read_ctx <- runtime_env$.analysis_run_context(
    analysis_key = c("sim", "bw", "est", "base"),
    run_id = "analysis-339-read-only",
    params = list(
      sim_grid_n_chunks = 1L,
      sim_grid_shuffle_seed = 20260707L
    )
  )
  discovered <- runtime_env$.find_bw_list_output_files(
    output_dir = read_ctx$staging_run_dir,
    cache_dir = read_ctx$current_dir,
    allow_cache_fallback = TRUE
  )

  expect_length(discovered, 1L)
  expect_equal(readRDS(discovered[[1]])$sim_id, 1L)
})
