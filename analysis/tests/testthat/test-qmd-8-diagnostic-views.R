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

# Baseline scenarios from QMD 8's settings chunk.
.qmd8_view_baselines <- function(text, env) {
  source(file.path(root_dir, "scripts", "r", "sim-misc.R"), local = env)
  env$root_dir <- root_dir
  env$analysis_quick <- FALSE
  env$run_simulations <- FALSE
  env$sim_size <- "final"
  chunk <- .qmd8_view_extract(text, "(?s)#\\| label: actual-settings.*?```")
  expect_true(nzchar(chunk))
  eval(parse(text = sub("```$", "", chunk)), env)
  env$baseline_scenarios_tbl
}

test_that("QMD 8 crosses every baseline setting for every transformation", {
  text <- .qmd8_view_text()
  env <- .qmd8_view_env()
  grid <- .qmd8_view_baselines(text, env)
  expect_identical(grid$base_scenario_id, seq_len(nrow(grid)))
  expect_equal(nrow(grid), 42L)
  expect_setequal(unique(grid$transformation), c("gaussian", "skew", "gamma"))
  expect_setequal(unique(grid$mean_pos_setting), c("high", "low"))
  # Every transformation and separation has the same cells and responses.
  settings <- split(
    paste(grid$n_cell, grid$prob_response),
    paste(grid$transformation, grid$mean_pos_setting)
  )
  expect_length(settings, 6L)
  for (x in settings) expect_identical(x, settings[[1]])
  expect_setequal(settings[[1]], c(
    paste(5e3, c(0.002, 0.01, 0.2)), paste(1e5, c(0.0005, 0.002, 0.01, 0.2))
  ))
  expect_identical(
    grid$scenario_desc[grid$transformation == "gaussian" &
      grid$mean_pos_setting == "high" & grid$n_cell == 5e3 &
      grid$prob_response == 0.01],
    "high separation, 5,000 cells, 1% response"
  )
  # Per-scenario views keep response probabilities apart; heat maps exist.
  for (label in c("relative-error", "signed-error")) {
    expect_true(grepl(paste0("label: ", label, "-by-n-cell\n"), text, fixed = TRUE))
  }
  for (label in c("f1-heatmap", "frequency-heatmap")) {
    expect_true(grepl(paste0("label: ", label, "\n"), text, fixed = TRUE))
  }
  expect_false(grepl("group_cell", text, fixed = TRUE))
  expect_false(grepl(".compare_loop_subtitle", text, fixed = TRUE))
})

test_that("figure loops put settings in headings and distinct save paths", {
  env <- .qmd8_view_env()
  plots <- list()
  paths <- character()
  env$.analysis_save_fig <- function(plot, path, ...) {
    paths <<- c(paths, path)
  }
  env$.analysis_print_fig <- function(p, ...) {
    plots[[length(plots) + 1L]] <<- p
  }
  data <- tidyr::expand_grid(
    method = c("stimgate", "fbeta", "tailgate"),
    mean_pos_setting = c("low", "high"),
    mismatch_type = c("mean_shift_all", "sd_inflation")
  )
  headings <- utils::capture.output(env$.simCompareFigureLoop(
    data, function(d) ggplot2::ggplot(d), dir = tempdir(),
    file_fn = function(pos, extra) paste0(pos, "_", extra, ".pdf"),
    height = 6, level = 4L, extra_col = "mismatch_type"
  ))
  # One method set (all methods): the data have no default-settings Tailgate.
  expect_length(plots, 4L)
  expect_length(unique(paths), 4L)
  expect_true(all(vapply(plots, function(p) {
    is.null(p$labels$title) && is.null(p$labels$subtitle)
  }, logical(1))))
  expect_false(any(grepl("Without Tailgate", headings, fixed = TRUE)))
  expect_true(any(grepl("All methods", headings, fixed = TRUE)))
  expect_true(any(grepl("Mean position: high", headings, fixed = TRUE)))
  expect_true(any(grepl("sd_inflation", headings, fixed = TRUE)))
  expect_true(any(grepl("high_sd_inflation.pdf", paths, fixed = TRUE)))
})

test_that("QMD 8 prose names the actual baseline cohort", {
  env <- .qmd8_view_env()
  text <- .qmd8_view_text()
  grid <- .qmd8_view_baselines(text, env)
  env$compare_raw <- data.frame(base_scenario_id = grid$base_scenario_id)
  env$run_plots <- TRUE
  env$results_available <- TRUE
  chunk <- .qmd8_view_extract(text, "(?s)#\\| label: figure-helpers.*?```")
  eval(parse(text = sub("```$", "", chunk)), env)
  gamma_high <- grid[grid$transformation == "gamma" & grid$mean_pos_setting == "high", ]
  small <- gamma_high$base_scenario_id[gamma_high$n_cell == 5e3][[1]]
  large <- gamma_high$base_scenario_id[gamma_high$n_cell == 1e5][[1]]
  data <- data.frame(transformation = "gamma", mean_pos_setting = "high",
    mismatch_type = "mean_shift_negative")
  note <- paste(utils::capture.output(env$.compare_baseline_note(data)), collapse = "\n")
  expect_match(note, "Baselines:", fixed = TRUE)
  expect_identical(env$.compare_mismatch_file("high",
    list(mismatch_type = "mean_shift_negative", n_cell = 5000)),
    "plot_mean_shift_negative_mean_pos_setting_high_n_cell_5000.pdf")
  expect_match(env$.compare_mismatch_heading(
    list(mismatch_type = "mean_shift_negative", n_cell = 5000)),
    "Shift stimulated negatives only", fixed = TRUE)
  expect_match(note, paste0(small, " (gamma"), fixed = TRUE)
  expect_match(note, paste0(large, " (gamma"), fixed = TRUE)
  data$n_cell <- 5000
  note <- paste(utils::capture.output(env$.compare_baseline_note(data)), collapse = "\n")
  expect_match(note, paste0(small, " (gamma"), fixed = TRUE)
  expect_false(grepl(paste0(large, " (gamma"), note, fixed = TRUE))
  # A dev run can omit a baseline even when its averaged columns match.
  env$compare_raw <- data.frame(base_scenario_id = small)
  data$n_cell <- NULL
  note <- paste(utils::capture.output(env$.compare_baseline_note(data)), collapse = "\n")
  expect_false(grepl(paste0(large, " (gamma"), note, fixed = TRUE))
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
  code <- .qmd8_view_extract(text, "(?s)diag_stimgate <-.*?(?=\n  \\.analysis_report_table)")
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
  expect_null(p$labels$title)
  expect_null(p$labels$subtitle)
  text <- .qmd8_view_text()
  expect_match(text, 'paste0("Sample ", s)', fixed = TRUE)
  expect_match(text, 'paste0("plot_sample_", s, ".pdf")', fixed = TRUE)
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

test_that("QMD 8 heat maps summarise each scenario and keep missing tiles grey", {
  env <- .qmd8_view_env()
  raw <- tidyr::expand_grid(
    method = c("stimgate", "fbeta"), transformation = c("gaussian", "gamma"),
    prob_response = c(0.01, 0.2), mismatch_val = c(0, 0.05), sample = 1:3
  ) |>
    dplyr::mutate(
      mean_pos_setting = "high", n_cell = 5000, mismatch_type = "mean_shift_all",
      propRespTruth = .data$prob_response,
      propRespEst = .data$prob_response * c(0.5, 1, 1.5)[.data$sample],
      error = NA_character_
    )
  # A run error has no estimate; a negative estimate stays in the median.
  raw$error[raw$method == "fbeta" & raw$sample == 3L] <- "failed"
  raw$propRespEst[raw$method == "fbeta" & raw$sample == 1L] <- -0.001
  cols <- c("transformation", "mean_pos_setting", "prob_response", "n_cell",
    "mismatch_type", "mismatch_val", "method")
  freq <- env$.simCompareHeatmapFrequencySummary(raw, cols)
  expect_equal(nrow(freq), 16L)
  sg <- freq[freq$method == "stimgate", ]
  expect_equal(sg$est_median, sg$prob_response)
  expect_equal(sg$ratio_median, rep(1, nrow(sg)))
  fb <- freq[freq$method == "fbeta", ]
  expect_equal(fb$n_est, rep(2L, nrow(fb)))
  expect_equal(fb$est_median, (-0.001 + fb$prob_response) / 2)
  p <- env$.simCompareHeatmapPlot(freq, "frequency")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 4L)
  expect_null(p$labels$title)
  expect_setequal(levels(p$data$y), c("1%", "20%"))
  # The F1 map prints values and leaves undefined tiles grey.
  f1 <- freq |> dplyr::mutate(f1_median = ifelse(.data$method == "fbeta", NA, 0.9))
  p1 <- env$.simCompareHeatmapPlot(f1, "f1")
  expect_true(all(p1$data$label[p1$data$method == "F-beta"] == "NA"))
  expect_true(all(p1$data$label[p1$data$method == "StimGate"] == "0.90"))
  expect_no_error(ggplot2::ggplot_build(p1))
})
