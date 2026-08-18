root_dir <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = FALSE)
if (!file.exists(file.path(root_dir, "scripts", "r", "sim-bandwidth.R"))) {
  root_dir <- normalizePath(
    file.path(testthat::test_path(), "../../.."),
    mustWork = FALSE
  )
}

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

test_that("scripts/r helpers source without error in dependency order", {
  skip_if_not(file.exists(script_misc), "sim-misc.R not found")
  skip_if_not(file.exists(script_bw), "sim-bandwidth.R not found")
  skip_if_not(file.exists(script_comp), "sim-compare-freq_bs.R not found")

  env <- new.env(parent = baseenv())
  expect_no_error(source(script_misc, local = env))
  expect_no_error(source(script_bw, local = env))
  expect_no_error(source(script_comp, local = env))
})

test_that("QMD analysis scripts do not call stimgate::: for scripts/r helpers", {
  qmd_dir <- file.path(root_dir, "analysis")
  skip_if_not(dir.exists(qmd_dir), "analysis/ directory not found")

  qmd_files <- list.files(qmd_dir, pattern = "\\.qmd$", full.names = TRUE)
  skip_if(length(qmd_files) == 0, "No .qmd files found")

  helper_fns <- c(
    "simBandwidthBwOne", "simBandwidthBsFreq", "simCompareStimgateRows",
    "simCompareFreqBs", "simBandwidthBwAll"
  )
  pattern <- paste0(
    "stimgate:::\\.(", paste(helper_fns, collapse = "|"), ")"
  )

  violations <- character(0)
  for (f in qmd_files) {
    lines <- readLines(f, warn = FALSE)
    hits <- grep(pattern, lines, value = TRUE)
    if (length(hits) > 0) {
      violations <- c(violations, paste0(basename(f), ": ", hits))
    }
  }
  expect_length(
    violations, 0L,
    label = paste("stimgate::: calls for scripts/r helpers in QMDs:", paste(violations, collapse = "; "))
  )
})
