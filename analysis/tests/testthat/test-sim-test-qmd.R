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
  expect_gt(nrow(env$test_grid), 0L)
  expect_true(is.numeric(env$test_settings$nSample))
  expect_gte(env$test_settings$nSample, 1L)
  configured_n_sample <- env$test_settings$nSample

  env$analysis_quick <- TRUE
  env$analysis_dev <- FALSE
  quiet(eval(chunk("test-grid"), envir = env))
  expect_identical(env$test_settings$nSample, min(configured_n_sample, 2L))
  expect_identical(nrow(env$test_grid), 1L)
  mean_settings <- env$.simMiscGetMeanPosTbl()
  expected_mean <- mean_settings$mean_pos[
    mean_settings$transformation == env$test_grid$transformation &
      mean_settings$mean_pos_setting == env$test_grid$mean_pos_setting
  ]
  expect_length(expected_mean, 1L)
  expect_equal(env$test_grid$mean_pos, expected_mean)

  env$run_simulations <- TRUE
  quiet(eval(chunk("test-run"), envir = env))
  expected_samples <- if (env$test_subsequent) {
    seq.int(env$test_sample, env$test_settings$nSample)
  } else {
    env$test_sample
  }
  expect_length(env$test_results, length(expected_samples))
  expect_identical(
    unname(vapply(env$test_results, function(x) x$dbg$sample, integer(1L))),
    expected_samples
  )
  dbg <- env$test_results[[1L]]$dbg
  expect_identical(dbg$ind, as.character(2L * env$test_sample))

  ref <- NULL
  quiet(ref <- env$.simBandwidthRunRow(
    env$test_grid[1, ],
    env$.simBandwidthTestScenario,
    env$test_settings
  ))
  for (result in env$test_results) {
    expect_equal(
      result$dbg$cp$cp,
      ref$threshold[ref$method == "loc_sample" & ref$ind == result$dbg$ind]
    )
  }
  table_root <- .local_projr_root()
  env$root_dir <- table_root
  env$fig_key <- c("2c-sim-test", "quick")
  env$run_plots <- TRUE
  output <- paste(utils::capture.output(summary <- eval(chunk("test-summary"), envir = env)), collapse = "\n")
  expect_s3_class(summary, "data.frame")
  expect_equal(nrow(summary), length(env$test_results))
  path_summary <- .projr_output_path(table_root, "table", "2c-sim-test", "quick", "test-summary.csv")
  expect_true(file.exists(path_summary))
  expect_equal(nrow(readr::read_csv(path_summary, show_col_types = FALSE)), nrow(summary))
  expect_match(output, "output/table/2c-sim-test/quick/test-summary.csv", fixed = TRUE)
  env$run_plots <- FALSE
  env$.analysis_report_table <- function(...) stop("Disabled summary wrote a CSV")
  expect_no_error(eval(chunk("test-summary"), envir = env))

  # A stimulated negative-cell shift and a per-row override change only the
  # stimulated tube; the unstimulated cells are the same draws.
  eval(chunk("test-settings"), envir = env)
  quiet(eval(chunk("test-grid"), envir = env))
  # Compare the same chosen sample across the baseline and shifted rows.
  env$test_subsequent <- FALSE
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
