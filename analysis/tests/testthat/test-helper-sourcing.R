root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_cyt <- file.path(root_dir, "scripts", "r", "functionsForBenchmarking-Cyt.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")
script_trans <- file.path(root_dir, "scripts", "r", "sim-trans.R")

test_that("scripts/r helpers source without error in dependency order", {
  for (f in c(script_misc, script_cyt, script_bw, script_comp, script_trans)) {
    if (!file.exists(f)) {
      stop("Expected analysis helper not found: ", f)
    }
  }

  env <- new.env(parent = getNamespace("stimgate"))
  expect_no_error(source(script_misc, local = env))
  expect_no_error(source(script_cyt, local = env))
  expect_no_error(source(script_bw, local = env))
  expect_no_error(source(script_comp, local = env))
  expect_no_error(source(script_trans, local = env))
})

test_that("QMD analysis scripts do not call scripts/r helpers via stimgate:::", {
  qmd_dir <- file.path(root_dir, "analysis")
  if (!dir.exists(qmd_dir)) stop("analysis/ directory not found at: ", qmd_dir)

  qmd_files <- list.files(qmd_dir, pattern = "\\.qmd$", full.names = TRUE)
  if (length(qmd_files) == 0) stop("No .qmd files found under: ", qmd_dir)

  # Derive the full set of helper function names from the sourced scripts
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_cyt, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)
  source(script_trans, local = env)
  helper_names <- sub("^\\.", "", ls(env, all.names = TRUE, pattern = "^\\.sim"))

  if (length(helper_names) == 0) {
    stop("No .sim* helpers found in sourced scripts — check sourcing order")
  }

  pattern <- paste0(
    "stimgate:::\\.(", paste(helper_names, collapse = "|"), ")"
  )

  violations <- character(0)
  for (f in qmd_files) {
    lines <- readLines(f, warn = FALSE)
    hits <- grep(pattern, lines, value = TRUE)
    if (length(hits) > 0) {
      violations <- c(violations, paste0(basename(f), ": ", trimws(hits)))
    }
  }
  expect_identical(
    violations, character(0),
    info = paste(
      "stimgate::: calls for scripts/r helpers in QMDs:",
      paste(violations, collapse = "; ")
    )
  )
})
