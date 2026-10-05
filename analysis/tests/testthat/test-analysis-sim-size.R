.sim_size_root <- function() {
  normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
}

.sim_size_runtime_env <- function() {
  env <- new.env(parent = baseenv())
  source(
    file.path(.sim_size_root(), "scripts", "r", "analysis-runtime.R"),
    local = env
  )
  env
}

test_that("sim_size defaults to draft and follows the QMD param and SIM_SIZE", {
  env <- .sim_size_runtime_env()
  withr::local_envvar(PROJR_PROFILE = NA, SIM_SIZE = NA)
  expect_identical(env$.analysis_sim_size(), "draft")
  expect_identical(env$.analysis_mode_key(c("sim", "default")), c("sim", "default", "draft"))
  env$params <- list(sim_size = "draft")
  expect_identical(env$.analysis_sim_size(), "draft")
  Sys.setenv(SIM_SIZE = "final")
  expect_identical(env$.analysis_sim_size(), "final")
  env$params <- list(sim_size = "final")
  Sys.setenv(SIM_SIZE = " Draft ")
  expect_identical(env$.analysis_sim_size(), "draft")
})

test_that("invalid sim_size values give a clear error", {
  env <- .sim_size_runtime_env()
  withr::local_envvar(PROJR_PROFILE = NA, SIM_SIZE = "medium")
  expect_error(env$.analysis_sim_size(), "SIM_SIZE", fixed = TRUE)
  expect_error(env$.analysis_sim_size(), "\"final\" or \"draft\"", fixed = TRUE)
  Sys.unsetenv("SIM_SIZE")
  env$params <- list(sim_size = "full")
  expect_error(env$.analysis_sim_size(), "not \"full\"", fixed = TRUE)
})

test_that("ACS figure paths ignore simulation sizes and retain explicit profiles", {
  env <- .sim_size_runtime_env()
  withr::local_envvar(PROJR_PROFILE = NA, SIM_SIZE = NA)
  stems <- c("9-real-compare-acs-cytof", "10-real-compare-acs-cytof-validation")
  for (stem in stems) {
    line <- grep("^fig_key <-", readLines(file.path(.sim_size_root(), "analysis", paste0(stem, ".qmd"))), value = TRUE)
    for (size in c(NA_character_, "draft", "final")) {
      if (is.na(size)) Sys.unsetenv("SIM_SIZE") else Sys.setenv(SIM_SIZE = size)
      eval(parse(text = line), envir = env)
      expect_identical(env$fig_key, stem)
    }
    for (profile in c("quick", "dev")) {
      Sys.setenv(PROJR_PROFILE = profile)
      eval(parse(text = line), envir = env)
      expect_identical(env$fig_key, c(stem, profile))
    }
    Sys.unsetenv("PROJR_PROFILE")
  }
})

test_that("dev and quick runs take precedence over draft sizes", {
  env <- .sim_size_runtime_env()
  key <- c("sim", "test")
  withr::local_envvar(PROJR_PROFILE = NA, SIM_SIZE = "draft")
  expect_identical(env$.analysis_sim_size(), "draft")
  expect_identical(env$.analysis_mode_key(key), c(key, "draft"))
  expect_identical(env$.analysis_mode_key(key, sized = FALSE), key)
  expect_identical(env$.analysis_mode_key("2a-test"), c("2a-test", "draft"))
  Sys.setenv(PROJR_PROFILE = "quick")
  expect_identical(env$.analysis_sim_size(), "final")
  expect_identical(env$.analysis_mode_key(key), c(key, "quick"))
  Sys.setenv(PROJR_PROFILE = "dev")
  expect_identical(env$.analysis_sim_size(), "final")
  expect_identical(env$.analysis_mode_key(key), c(key, "dev"))
  Sys.unsetenv("PROJR_PROFILE")
  Sys.setenv(SIM_SIZE = "final")
  expect_identical(env$.analysis_mode_key(key), key)
  # Invalid values fail even when the key would not use them.
  Sys.setenv(SIM_SIZE = "bad")
  expect_error(env$.analysis_mode_key(key), "SIM_SIZE", fixed = TRUE)
})

test_that("cache errors for draft results name the SIM_SIZE setting", {
  env <- .sim_size_runtime_env()
  expect_error(
    env$.analysis_cache_error(c("sim", "test", "draft"), "Missing.", "analysis/x.qmd"),
    "RUN_SIMULATIONS=true RUN_PLOTS=false SIM_SIZE=draft quarto render analysis/x.qmd",
    fixed = TRUE
  )
  expect_error(
    env$.analysis_cache_error(c("sim", "test"), "Missing.", "analysis/x.qmd"),
    "RUN_SIMULATIONS=true RUN_PLOTS=false SIM_SIZE=final quarto render analysis/x.qmd",
    fixed = TRUE
  )
})

test_that("a final render rejects draft results and accepts older final results", {
  env <- .sim_size_runtime_env()
  current <- withr::local_tempdir()
  file.create(file.path(current, "COMPLETE"))
  saveRDS(1L, file.path(current, "result.rds"))
  ctx <- list(analysis_key = c("sim", "test"), current_dir = current,
              qmd_path = "analysis/test.qmd")
  write_manifest <- function(params) {
    saveRDS(list(analysis_key = ctx$analysis_key, params = params),
            file.path(current, "manifest.rds"))
  }
  final <- list(version = 1L, sim_size = "final")
  draft <- list(version = 1L, sim_size = "draft")

  write_manifest(draft)
  expect_error(
    env$.analysis_current_file(ctx, "result.rds", final),
    "required manifest parameters: sim_size", fixed = TRUE
  )
  expect_true(file.exists(env$.analysis_current_file(ctx, "result.rds", draft)))

  write_manifest(final)
  expect_error(
    env$.analysis_current_file(ctx, "result.rds", draft),
    "required manifest parameters: sim_size", fixed = TRUE
  )
  # Results saved before sim_size was recorded were full-size runs.
  write_manifest(list(version = 1L))
  expect_true(file.exists(env$.analysis_current_file(ctx, "result.rds", final)))
  expect_error(
    env$.analysis_current_file(ctx, "result.rds", draft),
    "sim_size", fixed = TRUE
  )
})

test_that("each simulation QMD sets final and draft sample counts in one place", {
  root <- .sim_size_root()
  # Expected counts: samples per dataset and datasets per scenario.
  expected <- list(
    `2a` = list(final = c(200L, NA), draft = c(50L, NA), quick = c(1L, NA)),
    `2b` = list(final = c(25L, NA), draft = c(10L, NA), quick = c(1L, NA)),
    `3` = list(final = c(50L, 1L), draft = c(15L, 1L), quick = c(1L, 1L)),
    `4` = list(final = c(10L, 1L), draft = c(5L, 1L), quick = c(1L, 1L)),
    `5` = list(final = c(5L, 5L), draft = c(5L, 2L), quick = c(5L, 5L)),
    `6` = list(final = c(5L, 5L), draft = c(5L, 2L), quick = c(5L, 5L)),
    `7` = list(final = c(20L, 20L), draft = c(20L, 5L), quick = c(1L, 1L)),
    `8` = list(final = c(20L, 20L), draft = c(20L, 5L), quick = c(1L, 1L))
  )
  for (id in names(expected)) {
    file <- list.files(file.path(root, "analysis"), paste0("^", id, "-.*qmd$"),
                       full.names = TRUE)
    expect_length(file, 1L)
    lines <- readLines(file, warn = FALSE)
    content <- paste(lines, collapse = "\n")
    expect_true(grepl("\n  sim_size: draft\n", content, fixed = TRUE), info = id)
    expect_true(grepl("sim_size <- .analysis_sim_size()", content, fixed = TRUE),
                info = id)
    expect_true(grepl("  sim_size = sim_size,", content, fixed = TRUE), info = id)
    expect_true(grepl("SIM_SIZE=draft", content, fixed = TRUE), info = id)
    size_lines <- grep("^n_(sample|iter)_sim <- ", lines, value = TRUE)
    expect_true(any(grepl("[[sim_size]]", size_lines, fixed = TRUE)), info = id)
    counts <- function(size, quick = FALSE) {
      env <- new.env(parent = baseenv())
      env$sim_size <- size
      env$analysis_quick <- quick
      env$analysis_dev <- FALSE
      eval(parse(text = size_lines), env)
      c(
        if (exists("n_sample_sim", env, inherits = FALSE)) env$n_sample_sim else NA,
        if (exists("n_iter_sim", env, inherits = FALSE)) env$n_iter_sim else NA
      )
    }
    for (size in c("final", "draft")) {
      expect_equal(counts(size), expected[[id]][[size]], info = paste(id, size))
    }
    # Under quick, `.analysis_sim_size()` returns "final" (tested above).
    expect_equal(counts("final", quick = TRUE), expected[[id]]$quick, info = id)
  }
  # Analyses 7 and 8 gate each dataset together, so draft reduces datasets:
  # a quarter of the final number, at least 5 and never more than final.
  for (id in c("7", "8")) {
    final <- expected[[id]]$final[[2L]]
    expect_identical(expected[[id]]$draft[[2L]],
                     as.integer(min(final, max(5L, round(final / 4)))))
  }
})

test_that("analysis 1 has no sample count, so draft runs share its results", {
  content <- paste(readLines(file.path(
    .sim_size_root(), "analysis", "1-sim-trans.qmd"
  ), warn = FALSE), collapse = "\n")
  expect_true(grepl(".analysis_mode_key(analysis_key, sized = FALSE)", content,
                    fixed = TRUE))
  expect_true(grepl('.analysis_mode_key("1-sim-trans", sized = FALSE)', content,
                    fixed = TRUE))
})

test_that("a draft render of a simulation QMD reads its own folders", {
  root <- .sim_size_root()
  lines <- readLines(file.path(root, "analysis", "2a-sim-bw-freq_bs-global.qmd"))
  start <- grep("^analysis_key <- c\\(", lines)
  end <- grep("^fig_key <- ", lines)
  code <- lines[start:end]
  env <- .sim_size_runtime_env()
  withr::local_envvar(PROJR_PROFILE = NA, SIM_SIZE = "draft")
  eval(parse(text = code), env)
  expect_identical(env$sim_size, "draft")
  expect_identical(env$analysis_key, c("sim", "bw", "freq_bs", "global", "draft"))
  expect_identical(env$fig_key, c("2a-sim-bw-freq_bs-global", "draft"))
  Sys.setenv(PROJR_PROFILE = "quick")
  eval(parse(text = code), env)
  expect_identical(env$sim_size, "final")
  expect_identical(env$analysis_key, c("sim", "bw", "freq_bs", "global", "quick"))
})
