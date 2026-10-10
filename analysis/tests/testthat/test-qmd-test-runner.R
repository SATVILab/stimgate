.load_qmd_test_targets <- function() {
  root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  env <- new.env(parent = baseenv())
  source(file.path(root_dir, "analysis", "tests", "qmd-test-targets.R"), local = env)
  list(env = env, root_dir = root_dir)
}

test_that("QMD registry covers all top-level documents with real test files", {
  loaded <- .load_qmd_test_targets()
  targets <- loaded$env$.qmd_test_targets(loaded$root_dir)
  expect_length(targets, 21L)
  expect_setequal(names(targets), list.files(
    file.path(loaded$root_dir, "analysis"), pattern = "[.]qmd$"
  ))
  expect_identical(targets[[9]], c(
    "test-signed-percentile-plots.R",
    "test-acs-cytof-gate.R", "test-acs-cytof-methods.R", "test-ratio-companion-plots.R",
    "test-acs-cytof-paths.R"
  ))
  for (target in c(2L, 3L, 7L, 8L, 9L)) {
    expect_true("test-ratio-companion-plots.R" %in% targets[[target]])
  }
  expect_true("test-sim-bandwidth-analysis-run.R" %in% targets[[2]])
  expect_true("test-sim-bandwidth-analysis-run.R" %in% targets[[3]])
  expect_true("test-bandwidth-bias-plots-by-cells.R" %in% targets[[2]])
  expect_true("test-bandwidth-bias-plots-by-cells.R" %in% targets[[3]])
  expect_true("test-sim-bw-est-base-run.R" %in% targets[[4]])
  expect_true("test-analysis-7-transactional-multichunk.R" %in% targets[[7]])
})

test_that("QMD selection accepts defaults, aliases, sets and deduplicates", {
  loaded <- .load_qmd_test_targets()
  targets <- loaded$env$.qmd_test_targets(loaded$root_dir)
  select <- loaded$env$.select_qmd_test_targets
  expect_identical(select(character(), targets), targets)
  expect_identical(select("all", targets), targets)
  expect_identical(select(c("1,3", "9"), targets), targets[c(1, 4, 9)])
  expect_identical(select("1 3 9", targets), targets[c(1, 4, 9)])
  expect_identical(select("2a,2b", targets), targets[2:3])
  expect_identical(select("analysis/2b-sim-bias_uns-freq_bs.qmd", targets), targets[3])
  expect_identical(select("2a", targets), targets[2])
  expect_identical(select("2b", targets), targets[3])
  expect_identical(select("analysis\\2a-sim-bw-freq_bs-global.qmd", targets), targets[2])
  expect_identical(select("analysis\\2b-sim-bias_uns-freq_bs.qmd", targets), targets[3])
  expect_identical(select(c("1", "1-sim-trans", "analysis/1-sim-trans.qmd"), targets), targets[1])
  expect_identical(select("analysis\\10-real-compare-acs-cytof-validation.qmd", targets), targets[10])
  expect_identical(select("2c", targets), targets[11])
  expect_identical(select("11", targets), targets[12])
  expect_identical(select("12", targets), targets[13])
  expect_identical(select("13", targets), targets[14])
})

test_that("QMD selection rejects explicit empty and unknown requests", {
  loaded <- .load_qmd_test_targets()
  targets <- loaded$env$.qmd_test_targets(loaded$root_dir)
  select <- loaded$env$.select_qmd_test_targets
  for (selection in list("", "  ", ",,,", c("1", ""))) {
    expect_error(select(selection, targets), "must not be empty")
  }
  for (selection in list("99", "missing.qmd", c("all", "1"), c("--list", "1"))) {
    expect_error(select(selection, targets), "Unknown QMD selection")
  }
})

test_that("QMD registry rejects missing files and unregistered documents", {
  loaded <- .load_qmd_test_targets()
  targets <- loaded$env$.qmd_test_targets(loaded$root_dir)
  fixture <- tempfile("qmd-registry-")
  withr::defer(unlink(fixture, recursive = TRUE))
  analysis_dir <- file.path(fixture, "analysis")
  test_dir <- file.path(analysis_dir, "tests", "testthat")
  dir.create(test_dir, recursive = TRUE)
  file.create(file.path(analysis_dir, names(targets)))
  expect_error(loaded$env$.qmd_test_targets(fixture), "Missing QMD test files")
  file.create(file.path(test_dir, unique(unlist(targets, use.names = FALSE))))
  expect_identical(loaded$env$.qmd_test_targets(fixture), targets)
  extra <- file.path(analysis_dir, "11-new-analysis.qmd")
  file.create(extra)
  expect_error(loaded$env$.qmd_test_targets(fixture), "unregistered: 11-new-analysis.qmd")
  unlink(extra)
  unlink(file.path(analysis_dir, names(targets)[1]))
  expect_error(loaded$env$.qmd_test_targets(fixture), "Missing: 1-sim-trans.qmd")
})
