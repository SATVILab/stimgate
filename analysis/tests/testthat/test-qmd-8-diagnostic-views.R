root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."), mustWork = TRUE
)

.qmd8_view_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  env
}

.qmd8_view_text <- function() {
  paste(readLines(file.path(root_dir, "analysis",
    "8-sim-compare-freq_bs-batch.qmd"), warn = FALSE), collapse = "\n")
}

.qmd8_view_extract <- function(text, pattern) {
  regmatches(text, regexpr(pattern, text, perl = TRUE))
}

test_that("QMD 8 removes only views duplicated by its baseline grid", {
  text <- .qmd8_view_text()
  env <- .qmd8_view_env()
  grid_code <- .qmd8_view_extract(text,
    "(?s)targeted_scenarios_tbl <- tibble::tribble\\(.*?\n\\)")
  expect_true(nzchar(grid_code))
  eval(parse(text = grid_code), env)
  grid <- env$targeted_scenarios_tbl
  counts <- grid |>
    dplyr::count(transformation, mean_pos_setting, n_cell)
  expect_true(all(counts$n == 1L))
  # Averaging within high skew/gamma combines distinct baselines.
  high <- grid[grid$mean_pos_setting == "high", ] |>
    dplyr::count(transformation)
  expect_equal(high$n[match(c("skew", "gamma"), high$transformation)], c(2L, 2L))
  for (label in c("relative-error", "signed-error")) {
    expect_false(grepl(paste0("label: ", label, "-by-n-cell-prob"), text, fixed = TRUE))
    expect_true(grepl(paste0("label: ", label, "-by-n-cell\n"), text, fixed = TRUE))
  }
  # Every shared loop opts into identifying its settings.
  loops <- strsplit(text, ".simCompareFigureLoop(", fixed = TRUE)[[1]][-1]
  expect_gt(length(loops), 0L)
  for (loop in loops) {
    chunk <- strsplit(loop, "```", fixed = TRUE)[[1]][1]
    expect_true(grepl("subtitle = .compare_loop_subtitle", chunk, fixed = TRUE))
  }
})

test_that("figure loop subtitles identify each subset and preserve the default", {
  env <- .qmd8_view_env()
  plots <- list()
  env$.analysis_save_fig <- function(...) invisible(NULL)
  env$.analysis_print_fig <- function(p, ...) {
    plots[[length(plots) + 1L]] <<- p
  }
  data <- tidyr::expand_grid(
    method = c("stimgate", "fbeta", "tailgate"),
    mean_pos_setting = c("low", "high"),
    mismatch_type = c("mean_shift_all", "sd_inflation")
  )
  make_plot <- function(d) ggplot2::ggplot(d) + ggplot2::labs(subtitle = "Original")
  run <- function(subtitle = NULL) {
    utils::capture.output(env$.simCompareFigureLoop(
      data, make_plot, dir = tempdir(), file_fn = function(...) "unused.pdf",
      height = 6, level = 4L, extra_col = "mismatch_type", subtitle = subtitle
    ))
  }
  run()
  expect_length(plots, 8L)
  expect_true(all(vapply(plots, function(p) identical(p$labels$subtitle, "Original"), logical(1))))
  plots <- list()
  run(function(set, pos, extra, d) {
    expect_equal(unique(d$mean_pos_setting), pos)
    expect_equal(unique(d$mismatch_type), extra)
    expect_true(all(d$method %in% set$methods))
    paste(set$heading, pos, extra, sep = "; ")
  })
  subtitles <- vapply(plots, function(p) p$labels$subtitle, character(1))
  expect_length(unique(subtitles), 8L)
  expect_true("Without Tailgate; high; sd_inflation" %in% subtitles)
})

test_that("QMD 8 subtitles name the actual baseline cohort", {
  env <- .qmd8_view_env()
  text <- .qmd8_view_text()
  eval(parse(text = .qmd8_view_extract(text,
    "(?s)targeted_scenarios_tbl <- tibble::tribble\\(.*?\n\\)")), env)
  env$targeted_scenarios_tbl$base_scenario_id <- seq_len(nrow(env$targeted_scenarios_tbl))
  env$compare_raw <- data.frame(base_scenario_id = 1:7)
  env$run_plots <- TRUE
  env$results_available <- TRUE
  chunk <- .qmd8_view_extract(text, "(?s)#\\| label: figure-helpers.*?```")
  eval(parse(text = sub("```$", "", chunk)), env)
  set <- list(heading = "All methods")
  data <- data.frame(transformation = "gamma", mismatch_type = "mean_shift_negative")
  subtitle <- env$.compare_loop_subtitle(set, "high", NA, data)
  expect_match(subtitle, "All methods; mean position: high", fixed = TRUE)
  expect_match(subtitle, "Shift stimulated negatives only", fixed = TRUE)
  expect_match(subtitle, "3 (gamma", fixed = TRUE)
  expect_match(subtitle, "7 (gamma", fixed = TRUE)
  data$n_cell <- 5000
  subtitle <- env$.compare_loop_subtitle(set, "high", NA, data)
  expect_match(subtitle, "3 (gamma", fixed = TRUE)
  expect_false(grepl("7 (gamma", subtitle, fixed = TRUE))
  # A dev run can omit a baseline even when its averaged columns match.
  env$compare_raw <- data.frame(base_scenario_id = 3L)
  data$n_cell <- NULL
  subtitle <- env$.compare_loop_subtitle(set, "high", NA, data)
  expect_false(grepl("7 (gamma", subtitle, fixed = TRUE))
})

test_that("gate diagnostic IDs stay character and sort in numeric order", {
  text <- .qmd8_view_text()
  env <- .qmd8_view_env()
  env$gate_diagnostic <- list(cells = data.frame(sample = c("10", "2", "1", "11")))
  code <- .qmd8_view_extract(text,
    "(?s)diagnostic_samples <- unique.*?diagnostic_samples\\[order\\(as.integer\\(diagnostic_samples\\)\\)\\]")
  eval(parse(text = code), env)
  expect_identical(env$diagnostic_samples, c("1", "2", "10", "11"))
  env$diag_primary <- data.frame(
    sim_id = 1L, sample = c("10", "2", "1", "11"), method = "stimgate",
    mismatch_type = "mean_shift_all", mismatch_val = 0.05, threshold = 1.5,
    gateReturnPoint = "calculated", fdp = 0, sensitivity = 1,
    selected_fraction = 0.2, gate_status = "selected"
  )
  env$diag_sample_gate <- data.frame(
    sim_id = 1L, sample = env$diag_primary$sample,
    sample_gate = 1.4, sample_gate_source = "calculated"
  )
  code <- .qmd8_view_extract(text, "(?s)diag_stimgate <-.*?(?=\n  print\\(knitr::kable)")
  eval(parse(text = code), env)
  expect_identical(env$diag_stimgate$sample, c("1", "2", "10", "11"))
})

test_that("gate diagnostic figures name their sample and mismatch settings", {
  env <- .qmd8_view_env()
  cells <- tidyr::expand_grid(
    sample = "10", mismatch_type = c("mean_shift_all", "mean_shift_negative"),
    mismatch_val = c(0, 0.05, 0.1), condition = c("stim", "unstim"), i = 1:10
  ) |>
    dplyr::mutate(expr = i / 10, label = "gn")
  gates <- tidyr::expand_grid(
    mismatch_type = c("mean_shift_all", "mean_shift_negative"),
    mismatch_val = c(0, 0.05, 0.1), method = c("stimgate", "fbeta", "tailgate")
  ) |>
    dplyr::mutate(threshold = 0.8)
  p <- env$.simComparePlotGateDiagnostic(cells, gates)
  expect_match(p$labels$title, "sample 10", fixed = TRUE)
  expect_match(p$labels$subtitle, "Shift all stimulated cells", fixed = TRUE)
  expect_match(p$labels$subtitle, "Shift stimulated negatives only", fixed = TRUE)
  expect_match(p$labels$subtitle, "shifts: 0, 0.05, 0.1", fixed = TRUE)
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("SD inflation routing uses the grid name and respects explicit clusters", {
  env <- .qmd8_view_env()
  env$.simCompareEnsureCurrentCheckout <- function(...) invisible(TRUE)
  captured <- NULL
  env$.simCompareFreqBs <- function(..., stimSdMultiplierClusters = NULL) {
    captured <<- stimSdMultiplierClusters
    stop("Captured routing")
  }
  row <- data.frame(sim_id = 1L, mismatch_type = "sd_inflation")
  suppressWarnings(env$.simCompareRunScenario(row, nSample = 1, nIter = 1))
  expect_identical(captured, "gn")
  row$stim_sd_multiplier_clusters <- "gn"
  suppressWarnings(env$.simCompareRunScenario(row, nSample = 1, nIter = 1))
  expect_identical(captured, "gn")
  row$stim_sd_multiplier_clusters <- NA_character_
  suppressWarnings(env$.simCompareRunScenario(row, nSample = 1, nIter = 1))
  expect_null(captured)
})
