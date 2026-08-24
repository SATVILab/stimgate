root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
script_gate <- file.path(root_dir, "scripts", "r", "acs_cytof-gate.R")
qmd_path <- file.path(root_dir, "analysis", "9-real-compare-acs-cytof.qmd")

.load_acs_gate_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_gate, local = env)
  env
}

test_that("ACS CyTOF batches contain one unstimulated and four stimulated samples", {
  env <- .load_acs_gate_env()

  expect_equal(
    env$.acsCytofBatchList(20L),
    list(1:5, 6:10, 11:15, 16:20)
  )
  expect_error(env$.acsCytofBatchList(19L), "multiple of five")
  expect_error(env$.acsCytofBatchList(0L), "at least 5")
})

test_that("the TCRgd tester and full population use separate output paths", {
  env <- .load_acs_gate_env()

  full <- env$.acsCytofPopulationPaths(
    pop = "tcrgd",
    pathFcsBase = "raw",
    pathGsBase = "gs",
    pathScratchBase = "scratch"
  )
  tester <- env$.acsCytofPopulationPaths(
    pop = "tcrgd",
    pathFcsBase = "raw",
    pathGsBase = "gs",
    pathScratchBase = "scratch",
    outputGroup = "tester"
  )

  expect_equal(full$fcs, tester$fcs)
  expect_equal(full$gs, file.path("gs", "tcrgd"))
  expect_equal(tester$gs, file.path("gs", "tester", "tcrgd"))
  expect_equal(
    tester$stimgate,
    file.path("scratch", "tester", "tcrgd", "stimgate")
  )
  expect_equal(
    tester$tailgate,
    file.path("scratch", "tester", "tcrgd", "tailgate", "result.rds")
  )
  expect_equal(
    tester$fbeta,
    file.path("scratch", "tester", "tcrgd", "fbeta", "result.rds")
  )
  expect_false(identical(full$stimgate, tester$stimgate))
})

test_that("the reusable population runner preserves the ACS StimGate contract", {
  env <- .load_acs_gate_env()
  runner_formals <- names(formals(env$.acsCytofRunPopulation))
  runner_body <- paste(deparse(body(env$.acsCytofRunPopulation)), collapse = "\n")

  expect_true(all(c(
    "pop", "runPreprocessing", "runMethods", "runPlots",
    "nSample", "biasUns", "outputGroup"
  ) %in% runner_formals))
  expect_true(grepl('gateCombn = "min"', runner_body, fixed = TRUE))
  expect_true(grepl('bwMtd = "hpi1"', runner_body, fixed = TRUE))
  expect_true(grepl("calcCytPosGates = TRUE", runner_body, fixed = TRUE))
})

test_that("analysis 9 uses one runner for the tester and configured populations", {
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl(
    'pop_vec <- c("tcrgd", "cd4", "cd8", "nk", "nk_pre", "b")',
    content,
    fixed = TRUE
  ))
  expect_gte(
    lengths(regmatches(
      content,
      gregexpr(".acsCytofRunPopulation(", content, fixed = TRUE)
    )),
    1L
  )
  expect_true(grepl(".acsCytofRunPopulationSafe(", content, fixed = TRUE))
  expect_true(grepl('outputGroup = "tester"', content, fixed = TRUE))
  expect_true(grepl("furrr::future_map", content, fixed = TRUE))
  expect_true(grepl(
    ".acsCytofRunComparisonMethods(",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "run-comparison-methods-sequentially",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    'methods = c("stimgate", "tailgate", "fbeta")',
    content,
    fixed = TRUE
  ))
})
