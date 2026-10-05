.hpc_safety_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "sim-trans.R", "sim-compare-freq_bs.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

test_that("submitted canonical reads reject a previous run while manual reads work", {
  env <- .hpc_safety_env()
  root <- withr::local_tempdir()
  current <- file.path(root, "current")
  dir.create(current)
  ctx <- list(current_dir = current, analysis_key = "test", qmd_path = "analysis/test.qmd")
  saveRDS(list(run_id = "old-run", analysis_key = "test"), file.path(current, "manifest.rds"))
  writeLines("complete", file.path(current, "COMPLETE"))
  saveRDS(42L, file.path(current, "result.rds"))
  withr::local_envvar(ANALYSIS_EXPECTED_RUN_ID = "new-run")
  expect_error(env$.analysis_read_current(ctx, "result.rds"), "Expected run_id 'new-run'.*old-run")
  Sys.setenv(ANALYSIS_EXPECTED_RUN_ID = "old-run")
  expect_identical(env$.analysis_read_current(ctx, "result.rds"), 42L)
  Sys.unsetenv("ANALYSIS_EXPECTED_RUN_ID")
  expect_identical(env$.analysis_read_current(ctx, "result.rds"), 42L)
})

test_that("transformation cache records submission provenance without changing scientific settings", {
  env <- .hpc_safety_env()
  path <- file.path(withr::local_tempdir(), "table.rds")
  settings <- list(seed = 42L)
  withr::local_envvar(c(ANALYSIS_RUN_ID = "new-run", ANALYSIS_EXPECTED_RUN_ID = "new-run"))
  env$sim_trans_write_cache(42L, settings, path)
  expect_identical(env$sim_trans_read_cache(path, settings), 42L)
  Sys.setenv(ANALYSIS_EXPECTED_RUN_ID = "another-run")
  expect_error(env$sim_trans_read_cache(path, settings), "Refusing to plot another run")
  Sys.unsetenv("ANALYSIS_EXPECTED_RUN_ID")
  expect_identical(env$sim_trans_read_cache(path, settings), 42L)
})

test_that("distinct runs of an analysis share the promotion lock but keep status locks separate", {
  env <- .hpc_safety_env()
  root <- normalizePath(withr::local_tempdir(), winslash = "/")
  withr::local_dir(root)
  contexts <- lapply(c("run-a", "run-b"), function(id) {
    env$.analysis_run_context("lock-test", run_id = id, path_root = root)
  })
  paths <- lapply(contexts, env$.analysis_lock_path, lock_name = "promotion")
  expect_identical(paths[[1]], paths[[2]])
  expect_identical(dirname(paths[[1]]), dirname(contexts[[1]]$current_dir))
  expect_false(identical(
    env$.analysis_lock_path(contexts[[1]], "status-update"),
    env$.analysis_lock_path(contexts[[2]], "status-update")
  ))
  # A second process using the other run context cannot acquire the shared lock.
  lock <- env$.analysis_acquire_lock(paths[[1]], timeout_sec = 1)
  expect_false(is.null(lock))
  withr::defer(env$.analysis_release_lock(lock))
  cl <- parallel::makePSOCKcluster(1L)
  withr::defer(parallel::stopCluster(cl))
  blocked <- parallel::clusterCall(cl, function(path) {
    lock <- filelock::lock(path, timeout = 100)
    if (!is.null(lock)) filelock::unlock(lock)
    is.null(lock)
  }, paths[[2]])
  expect_true(blocked[[1]])
  env$.analysis_release_lock(lock)
})

test_that("full comparison validation happens before promotion and preserves last good results", {
  env <- .hpc_safety_env()
  root <- normalizePath(withr::local_tempdir(), winslash = "/")
  withr::local_dir(root)
  ctx <- env$.analysis_run_context("validation-test", run_id = "new-run", path_root = root)
  dir.create(ctx$current_dir)
  saveRDS("last good", file.path(ctx$current_dir, "result.rds"))
  for (id in 1:2) {
    saveRDS(id, file.path(ctx$chunk_output_dir, sprintf("sim_scenario-sim_id_%06d.rds", id)))
  }
  full <- data.frame(sim_id = 1:2)
  env$.simCompareCollateScenarioOutputs <- function(...) full
  env$.simCompareGridOutputStatus <- function(...) list(collate_ok = TRUE, validation_ok = TRUE)
  env$.analysis_mark_chunk(ctx, 2L, 2L, 0L, collate_ok = TRUE, validation_ok = TRUE)
  expect_error(env$.simComparePromoteIfReady(
    ctx, full, 2L, 2L, 0L, 1L, 1L,
    validate_full = function(tbl) {
      expect_identical(tbl, full)
      stop("Cross-setting pairing failed")
    }
  ), "Cross-setting pairing failed")
  expect_identical(readRDS(file.path(ctx$current_dir, "result.rds")), "last good")
  expect_false(env$.analysis_can_promote(ctx))
  expect_false(file.exists(file.path(ctx$staging_collated_dir, "compare_raw.rds")))
})

test_that("Analysis 8 uses one complete settings list and pre-promotion mismatch validation", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  lines <- readLines(file.path(root, "analysis", "8-sim-compare-freq_bs-batch.qmd"))
  content <- paste(lines, collapse = "\n")
  expect_equal(sum(grepl("^analysis_result_params <- list", lines)), 1L)
  expect_true(grepl("required_params = analysis_result_params", content, fixed = TRUE))
  expect_true(grepl("validate_full = .simCompareValidateMismatch", content, fixed = TRUE))
  expect_true(grepl("batch-mismatch-comparison-v12", content, fixed = TRUE))
  start <- which(grepl("^analysis_result_params <- list", lines))
  end <- start + which(lines[(start + 1L):length(lines)] == ")")[[1L]]
  expr <- parse(text = lines[start:end])[[1L]][[3L]]
  expected <- c(
    "analysis_grid_spec", "gate_diagnostic_spec", "cluster_gates", "stimgate_bw_mtd",
    "stimgate_bw_scope", "stimgate_bw_ncell_max", "stimgate_bw_fallback", "stimgate_bw_min", "stimgate_bw_max",
    "loc_enforce_shape_threshold", "calc_cyt_pos_gates", "stimgate_bias_uns_factor",
    "fbeta_beta", "fbeta_theta", "fbeta_width", "tailgate_adjust", "tailgate_method",
    "tailgate_tol", "tailgate_auto_tol", "tailgate_x"
  )
  expect_true(all(expected %in% names(as.list(expr))))
  # Every required setting is enforced by canonical reads, including settings
  # whose changes do not alter the scenario IDs.
  env <- .hpc_safety_env()
  current <- file.path(withr::local_tempdir(), "current")
  dir.create(current)
  writeLines("complete", file.path(current, "COMPLETE"))
  saveRDS(42L, file.path(current, "result.rds"))
  required <- stats::setNames(as.list(seq_along(expected)), expected)
  manifest <- list(analysis_key = "test", params = required)
  ctx <- list(current_dir = current, analysis_key = "test")
  withr::local_envvar(ANALYSIS_EXPECTED_RUN_ID = NA_character_)
  for (nm in expected) {
    changed <- manifest
    changed$params[[nm]] <- -1L
    saveRDS(changed, file.path(current, "manifest.rds"))
    expect_error(env$.analysis_current_file(ctx, "result.rds", required), nm)
  }
})

test_that("mismatch validation rejects unpaired data and different zero-shift gates", {
  env <- .hpc_safety_env()
  full <- data.frame(
    sim_id = 1:2, base_scenario_id = 1L, iter = 1L, sample = "stim",
    method = "stimgate", mismatch_type = c("mean_shift_all", "mean_shift_negative"),
    mismatch_val = 0, unsExprSum = 42, nTruePos = 1L, nFalsePos = 1L,
    nFalseNeg = 1L, nTrueNeg = 7L, nPosStim = 2L, nCellStim = 10L,
    threshold = 3
  )
  expect_no_error(env$.simCompareValidateMismatch(full))
  unpaired <- full
  unpaired$unsExprSum[[2L]] <- 43
  expect_error(env$.simCompareValidateMismatch(unpaired), "settings are not paired")
  different_gate <- full
  different_gate$threshold[[2L]] <- 4
  expect_error(env$.simCompareValidateMismatch(different_gate), "did not.*reproduce")
  different_counts <- full
  different_counts$nTruePos[[2L]] <- 2L
  different_counts$nFalsePos[[2L]] <- 0L
  different_counts$nFalseNeg[[2L]] <- 0L
  different_counts$nTrueNeg[[2L]] <- 8L
  expect_error(env$.simCompareValidateMismatch(different_counts), "did not.*reproduce")
})

test_that("new unchunked cache schemas require their semantics identifier", {
  env <- .hpc_safety_env()
  withr::local_envvar(ANALYSIS_EXPECTED_RUN_ID = NA_character_)
  expect_error(
    env$.analysis_check_expected_run(list(), "acs_cytof",
      semantics_version = "acs-cytof-v2"),
    "analysis_semantics_version"
  )
  expect_error(env$.analysis_check_expected_run(
    list(analysis_semantics_version = "acs-cytof-v1"), "acs_cytof",
    semantics_version = "acs-cytof-v2"
  ), "analysis_semantics_version")
  expect_no_error(env$.analysis_check_expected_run(
    list(analysis_semantics_version = "acs-cytof-v2"), "acs_cytof",
    semantics_version = "acs-cytof-v2"
  ))
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  for (file in c("9-real-compare-acs-cytof.qmd", "10-real-compare-acs-cytof-validation.qmd")) {
    text <- paste(readLines(file.path(root, "analysis", file)), collapse = "\n")
    expect_match(text, 'semantics_version = "acs-cytof-v2"', fixed = TRUE)
  }
})
