root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

test_that("Analysis 2c runs chosen settings through the 2a code path", {
  testthat::skip_if_not_installed("simcyto")
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "analysis-plot-style.R", "sim-misc.R",
    "sim-bandwidth.R", "sim-bandwidth-analysis-io.R",
    "sim-bandwidth-analysis-run.R", "sim-compare-freq_bs.R",
    "sim-debug-loc.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  lines <- readLines(file.path(root_dir, "analysis", "2c-sim-test.qmd"))
  starts <- which(startsWith(lines, "```{r"))
  for (start in starts) {
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1L]
    expect_no_error(parse(text = lines[(start + 1L):(end - 1L)]))
  }
  chunk <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1L]
    parse(text = lines[(start + 1L):(end - 1L)])
  }
  quiet <- function(code) {
    utils::capture.output(suppressMessages(code))
    invisible()
  }

  eval(chunk("test-settings"), envir = env)
  expect_identical(nrow(env$test_grid), 1L)
  expect_identical(env$test_settings$nSample, 200L)

  env$analysis_quick <- TRUE
  env$analysis_dev <- FALSE
  quiet(eval(chunk("test-grid"), envir = env))
  expect_identical(env$test_settings$nSample, 2L)
  expect_identical(env$test_grid$mean_pos, 8.5)

  env$run_simulations <- TRUE
  quiet(eval(chunk("test-run"), envir = env))
  expect_length(env$test_results, 1L)
  dbg <- env$test_results[[1L]]$dbg
  expect_identical(dbg$ind, "2")

  ref <- NULL
  quiet(ref <- env$.simBandwidthRunRow(
    env$test_grid[1, ],
    env$.simBandwidthFreqBsGlobalScenario,
    env$test_settings
  ))
  expect_equal(
    dbg$cp$cp,
    ref$threshold[ref$method == "loc_sample" & ref$ind == "2"]
  )
  expect_s3_class(eval(chunk("test-summary"), envir = env), "knitr_kable")

  # A stimulated negative-cell shift and a per-row override change only the
  # stimulated tube; the unstimulated cells are the same draws.
  eval(chunk("test-settings"), envir = env)
  quiet(eval(chunk("test-grid"), envir = env))
  env$test_grid <- dplyr::bind_rows(
    env$test_grid,
    dplyr::mutate(env$test_grid, stim_mean_shift = 0.3, locMinPeakProb = 0.1)
  ) |>
    dplyr::mutate(sim_id = dplyr::row_number())
  quiet(eval(chunk("test-run"), envir = env))
  expect_length(env$test_results, 2L)
  base <- env$test_results[[1L]]$dbg
  shifted <- env$test_results[[2L]]$dbg
  expect_equal(base$cp$cp, dbg$cp$cp)
  expect_equal(shifted$inputs$chnlSettings$locMinPeakProb, 0.1)
  expect_equal(base$inputs$chnlSettings$locMinPeakProb, 0.25)
  expect_equal(
    shifted$inputs$exTblUnsOrig$F1,
    base$inputs$exTblUnsOrig$F1
  )
  expect_false(isTRUE(all.equal(
    shifted$inputs$exTblStimOrig$F1,
    base$inputs$exTblStimOrig$F1
  )))

  env$test_settings$notAnArgument <- 1
  expect_error(
    quiet(eval(chunk("test-grid"), envir = env)),
    "notAnArgument"
  )
})
