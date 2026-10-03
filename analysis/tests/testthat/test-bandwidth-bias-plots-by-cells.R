.bandwidth_cell_plot_chunk <- function(document, label) {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  lines <- readLines(file.path(root, "analysis", document), warn = FALSE)
  start <- which(lines == paste0("#| label: ", label))
  testthat::expect_length(start, 1L)
  end <- which(lines == "```" & seq_along(lines) > start)[1L]
  parse(text = lines[seq.int(start + 1L, end - 1L)])
}

.bandwidth_cell_plot_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  source(
    file.path(testthat::test_path(), "../../../scripts/r/sim-bandwidth-analysis-plot.R"),
    local = env
  )
  env$root_dir <- tempfile("cell-plots-")
  env$analysis_key <- "bias_uns"
  env$run_plots <- TRUE
  env$saved <- list()
  env$.analysis_project_dir <- function(type, subdir, root_dir) {
    path <- file.path(env$root_dir, "fig")
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
  }
  env$.analysis_cache_dir <- function(parts, path_root) {
    file.path(path_root, paste(parts, collapse = "/"))
  }
  env$ggsave <- function(filename, plot, ...) {
    env$saved[[length(env$saved) + 1L]] <- list(path = filename, plot = plot)
    invisible(NULL)
  }
  env$print <- function(x, ...) invisible(x)
  env
}

test_that("2a per-cell plots average scenario errors before combining probabilities", {
  env <- .bandwidth_cell_plot_env()
  on.exit(unlink(env$root_dir, recursive = TRUE), add = TRUE)
  env$bw_tbl_results_raw <- tidyr::expand_grid(
    transformation = c("gaussian", "gamma"),
    mean_pos_setting = c("low", "high"),
    n_cell = c(100, 1000), bw = c(0.1, 0.2),
    bias_uns_setting = c("low", "high"), condition_perturbation_sd = c(0, 0.5)
  ) |>
    dplyr::cross_join(tibble::tibble(
      prob_response = c(0.01, 0.01, 0.1), err = c(1, 3, 9)
    )) |>
    dplyr::mutate(
      propRespTruth = prob_response,
      propRespEst = propRespTruth * (1 + err * ifelse(n_cell == 100, 1, 2)),
      propRespEst = ifelse(
        bias_uns_setting != "low" | condition_perturbation_sd != 0,
        propRespTruth * 1000, propRespEst
      )
    )
  code <- .bandwidth_cell_plot_chunk(
    "2a-sim-bw-freq_bs-global.qmd", "plot-relative-error-by-n-cell"
  )
  eval(code, env)
  expect_length(env$saved, 4L)
  paths <- vapply(env$saved, `[[`, character(1), "path")
  expect_equal(length(unique(paths)), 4L)
  for (saved in env$saved) {
    data <- saved$plot$data
    expect_length(unique(data$n_cell), 1L)
    expect_length(unique(data$mean_pos_setting), 1L)
    expect_setequal(unique(data$transformation), c("gaussian", "gamma"))
    expect_true(grepl(paste0("_n_cell_", data$n_cell[1L], ".pdf"), saved$path, fixed = TRUE))
    expected <- c(err_rel_median_avg = 550, err_rel_95_avg = 595, err_rel_max_avg = 600)
    multiplier <- if (data$n_cell[1L] == 100) 1 else 2
    expect_equal(data$err_value, unname(expected[data$err_type]) * multiplier)
    expect_named(saved$plot$facet$params$facets, "transformation")
    expect_identical(rlang::as_label(saved$plot$mapping$group), "err_type")
  }
  unlink(env$root_dir, recursive = TRUE)
  env$saved <- list()
  env$run_plots <- FALSE
  eval(code, env)
  expect_length(env$saved, 0L)
  expect_false(dir.exists(env$root_dir))
})

test_that("2b averaged and per-cell plots preserve all bias scenario dimensions", {
  env <- .bandwidth_cell_plot_env()
  on.exit(unlink(env$root_dir, recursive = TRUE), add = TRUE)
  env$bias_uns_results_summary <- tidyr::expand_grid(
    transformation = c("gaussian", "gamma"),
    mean_pos_setting = c("low", "high"), prob_response = c(0.01, 0.1),
    n_cell = c(100, 1000), bw = c(0.1, 0.2),
    bias_uns_basis = c("bandwidth", "negative_width"), bias_uns_multiplier = c(0, 1),
    mismatch_val = c(0, 0.1)
  ) |>
    dplyr::mutate(
      mismatch_type = "mean_shift",
      mismatch_label = paste0("mean shift ", mismatch_val),
      median_abs_rel_error = n_cell / 100 + bw + bias_uns_multiplier +
        mismatch_val + prob_response + ifelse(bias_uns_basis == "bandwidth", 0, 2)
    )
  average_code <- .bandwidth_cell_plot_chunk(
    "2b-sim-bias_uns-freq_bs.qmd", "plot-relative-error-averaged-n-cell"
  )
  cell_code <- .bandwidth_cell_plot_chunk(
    "2b-sim-bias_uns-freq_bs.qmd", "plot-relative-error-by-n-cell"
  )
  eval(average_code, env)
  expect_length(env$saved, 8L)
  for (saved in env$saved) {
    data <- saved$plot$data
    expect_false("n_cell" %in% names(data))
    expect_equal(data$median_abs_rel_error,
      5.5 + data$bw + data$bias_uns_multiplier + data$mismatch_val +
        data$prob_response + ifelse(data$bias_uns_basis == "bandwidth", 0, 2)
    )
  }
  eval(cell_code, env)
  expect_length(env$saved, 24L)
  paths <- vapply(env$saved, `[[`, character(1), "path")
  expect_equal(length(unique(paths)), 24L)
  for (saved in env$saved) {
    data <- saved$plot$data
    for (key in c("transformation", "mean_pos_setting", "prob_response")) {
      expect_length(unique(data[[key]]), 1L)
    }
    expect_setequal(unique(data$mismatch_val), c(0, 0.1))
    expect_setequal(unique(data$bw), c(0.1, 0.2))
    expect_setequal(unique(data$bias_uns_basis), c("bandwidth", "negative_width"))
    expect_named(saved$plot$facet$params$facets, "mismatch_label")
    expect_identical(rlang::as_label(saved$plot$mapping$group), "interaction(bw, bias_uns_basis)")
    if ("n_cell" %in% names(data)) {
      expect_length(unique(data$n_cell), 1L)
      expect_true(grepl(paste0("_n_cell_", data$n_cell[1L], ".pdf"), saved$path, fixed = TRUE))
      expect_equal(data$median_abs_rel_error,
        data$n_cell / 100 + data$bw + data$bias_uns_multiplier + data$mismatch_val +
          data$prob_response + ifelse(data$bias_uns_basis == "bandwidth", 0, 2)
      )
    }
  }
  unlink(env$root_dir, recursive = TRUE)
  env$saved <- list()
  env$run_plots <- FALSE
  eval(average_code, env)
  eval(cell_code, env)
  expect_length(env$saved, 0L)
  expect_false(dir.exists(env$root_dir))
})
