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
  env$printed <- list()
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
  env$print <- function(x, ...) {
    env$printed[[length(env$printed) + 1L]] <- x
    invisible(x)
  }
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
    "2a-sim-bw-freq_bs-global.qmd", "fig-relative-error-by-n-cell"
  )
  eval(code, env)
  expect_length(env$saved, 4L)
  expect_identical(env$printed, lapply(env$saved, `[[`, "plot"))
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
  env$printed <- list()
  env$run_plots <- FALSE
  eval(code, env)
  expect_length(env$saved, 0L)
  expect_length(env$printed, 0L)
  expect_false(dir.exists(env$root_dir))
})

test_that("2b averaged and per-cell plots preserve all bias scenario dimensions", {
  stat_mult <- c(Median = 1, "90th percentile" = 2, Maximum = 3)
  env <- .bandwidth_cell_plot_env()
  on.exit(unlink(env$root_dir, recursive = TRUE), add = TRUE)
  env$bias_uns_abs_error <- tidyr::expand_grid(
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
        mismatch_val + prob_response + ifelse(bias_uns_basis == "bandwidth", 0, 2),
      q90_abs_rel_error = 2 * median_abs_rel_error,
      max_abs_rel_error = 3 * median_abs_rel_error
    )
  average_code <- .bandwidth_cell_plot_chunk(
    "2b-sim-bias_uns-freq_bs.qmd", "fig-relative-error-averaged-n-cell"
  )
  cell_code <- .bandwidth_cell_plot_chunk(
    "2b-sim-bias_uns-freq_bs.qmd", "fig-relative-error-by-n-cell"
  )
  eval(average_code, env)
  expect_length(env$saved, 8L)
  for (saved in env$saved) {
    data <- saved$plot$data
    expect_false("n_cell" %in% names(data))
    expect_equal(data$value,
      unname(stat_mult[as.character(data$statistic)]) *
        (5.5 + data$bw + data$bias_uns_multiplier + data$mismatch_val +
          data$prob_response + ifelse(data$bias_uns_basis == "bandwidth", 0, 2))
    )
  }
  eval(cell_code, env)
  expect_length(env$saved, 24L)
  expect_identical(env$printed, lapply(env$saved, `[[`, "plot"))
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
    expect_named(saved$plot$facet$params$rows, "statistic")
    expect_named(saved$plot$facet$params$cols, "mismatch_label")
    expect_setequal(as.character(data$statistic), names(stat_mult))
    expect_identical(rlang::as_label(saved$plot$mapping$group), "interaction(bw, bias_uns_basis)")
    if ("n_cell" %in% names(data)) {
      expect_length(unique(data$n_cell), 1L)
      expect_true(grepl(paste0("_n_cell_", data$n_cell[1L], ".pdf"), saved$path, fixed = TRUE))
      expect_equal(data$value,
        unname(stat_mult[as.character(data$statistic)]) *
          (data$n_cell / 100 + data$bw + data$bias_uns_multiplier +
            data$mismatch_val + data$prob_response +
            ifelse(data$bias_uns_basis == "bandwidth", 0, 2))
      )
    }
  }
  unlink(env$root_dir, recursive = TRUE)
  env$saved <- list()
  env$printed <- list()
  env$run_plots <- FALSE
  eval(average_code, env)
  eval(cell_code, env)
  expect_length(env$saved, 0L)
  expect_length(env$printed, 0L)
  expect_false(dir.exists(env$root_dir))
})

test_that("2a per-cell signed-error plots average scenario errors by direction", {
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
    "2a-sim-bw-freq_bs-global.qmd", "fig-signed-error-by-n-cell"
  )
  eval(code, env)
  expect_length(env$saved, 4L)
  expect_identical(env$printed, lapply(env$saved, `[[`, "plot"))
  paths <- vapply(env$saved, `[[`, character(1), "path")
  expect_equal(length(unique(paths)), 4L)
  for (saved in env$saved) {
    data <- saved$plot$data
    expect_length(unique(data$n_cell), 1L)
    expect_length(unique(data$mean_pos_setting), 1L)
    expect_setequal(unique(data$transformation), c("gaussian", "gamma"))
    expect_true(grepl(paste0("_n_cell_", data$n_cell[1L], ".pdf"), saved$path, fixed = TRUE))
    # All errors are over-estimates, so no under-estimate curves are drawn.
    expect_true(all(data$direction == "over"))
    expect_equal(data$prop, rep(1, nrow(data)))
    expected <- c(median = 5.5, q95 = 5.95, max = 6)
    multiplier <- if (data$n_cell[1L] == 100) 1 else 2
    expect_equal(
      data$err_value,
      unname(expected[as.character(data$err_type)]) * multiplier
    )
    expect_named(saved$plot$facet$params$facets, "transformation")
    expect_identical(
      rlang::as_label(saved$plot$mapping$group), "interaction(err_type, direction)"
    )
  }
  unlink(env$root_dir, recursive = TRUE)
  env$saved <- list()
  env$printed <- list()
  env$run_plots <- FALSE
  eval(code, env)
  expect_length(env$saved, 0L)
  expect_length(env$printed, 0L)
  expect_false(dir.exists(env$root_dir))
})

test_that("2b signed-error plots preserve all bias scenario dimensions", {
  stat_mult <- c(Median = 1, "90th percentile" = 2, Maximum = 3)
  env <- .bandwidth_cell_plot_env()
  on.exit(unlink(env$root_dir, recursive = TRUE), add = TRUE)
  env$bias_uns_signed_error <- tidyr::expand_grid(
    transformation = c("gaussian", "gamma"),
    mean_pos_setting = c("low", "high"), prob_response = c(0.01, 0.1),
    n_cell = c(100, 1000), bw = c(0.1, 0.2),
    bias_uns_basis = c("bandwidth", "negative_width"), bias_uns_multiplier = c(0, 1),
    mismatch_val = c(0, 0.1), direction = c("over", "under")
  ) |>
    dplyr::mutate(
      mismatch_type = "mean_shift",
      mismatch_label = paste0("mean shift ", mismatch_val),
      sign = ifelse(direction == "over", 1, -0.01),
      prop = ifelse(direction == "over", 0.7, 0.3),
      median = sign * (n_cell / 100 + bw + bias_uns_multiplier +
        mismatch_val + prob_response + ifelse(bias_uns_basis == "bandwidth", 0, 2)),
      q90 = 2 * median,
      max = 3 * median
    ) |>
    dplyr::select(-sign)
  average_code <- .bandwidth_cell_plot_chunk(
    "2b-sim-bias_uns-freq_bs.qmd", "fig-signed-error-averaged-n-cell"
  )
  cell_code <- .bandwidth_cell_plot_chunk(
    "2b-sim-bias_uns-freq_bs.qmd", "fig-signed-error-by-n-cell"
  )
  eval(average_code, env)
  expect_length(env$saved, 8L)
  for (saved in env$saved) {
    data <- saved$plot$data
    expect_false("n_cell" %in% names(data))
    expect_setequal(data$direction, c("over", "under"))
    expect_equal(data$prop, ifelse(data$direction == "over", 0.7, 0.3))
    expect_equal(data$value,
      unname(stat_mult[as.character(data$statistic)]) *
        ifelse(data$direction == "over", 1, -0.01) *
        (5.5 + data$bw + data$bias_uns_multiplier + data$mismatch_val +
          data$prob_response + ifelse(data$bias_uns_basis == "bandwidth", 0, 2))
    )
  }
  eval(cell_code, env)
  expect_length(env$saved, 24L)
  expect_identical(env$printed, lapply(env$saved, `[[`, "plot"))
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
    expect_named(saved$plot$facet$params$rows, "statistic")
    expect_named(saved$plot$facet$params$cols, "mismatch_label")
    expect_setequal(as.character(data$statistic), names(stat_mult))
    expect_identical(
      rlang::as_label(saved$plot$mapping$group),
      "interaction(bw, bias_uns_basis, direction)"
    )
    if ("n_cell" %in% names(data)) {
      expect_length(unique(data$n_cell), 1L)
      expect_true(grepl(paste0("_n_cell_", data$n_cell[1L], ".pdf"), saved$path, fixed = TRUE))
      expect_equal(data$value,
        unname(stat_mult[as.character(data$statistic)]) *
          ifelse(data$direction == "over", 1, -0.01) *
          (data$n_cell / 100 + data$bw + data$bias_uns_multiplier +
            data$mismatch_val + data$prob_response +
            ifelse(data$bias_uns_basis == "bandwidth", 0, 2))
      )
    }
  }
  unlink(env$root_dir, recursive = TRUE)
  env$saved <- list()
  env$printed <- list()
  env$run_plots <- FALSE
  eval(average_code, env)
  eval(cell_code, env)
  expect_length(env$saved, 0L)
  expect_length(env$printed, 0L)
  expect_false(dir.exists(env$root_dir))
})

test_that("signed relative errors are summarised separately by direction", {
  env <- .bandwidth_cell_plot_env()
  sides <- env$.simBandwidthSignedErrorSides(c(-1, -0.5, 0.2, 1, NA, 0))
  expect_identical(sides$direction, c("over", "under"))
  # Zero errors count towards the total but neither direction.
  expect_equal(sides$prop, c(0.4, 0.4))
  expect_equal(sides$median, c(0.6, -0.75))
  expect_equal(sides$q90, c(0.92, -0.95))
  expect_equal(sides$q95, c(0.96, -0.975))
  expect_equal(sides$max, c(1, -1))

  empty <- env$.simBandwidthSignedErrorSides(c(0.5, 1))
  expect_equal(empty$prop, c(1, 0))
  expect_true(all(is.na(empty[empty$direction == "under", c("median", "q95", "max")])))

  summary <- env$.simBandwidthSignedErrorSummary(
    tibble::tibble(g = c(1, 1, 2, 2), rel_error = c(-0.5, 1, 3, 1)),
    "g"
  )
  avg <- env$.simBandwidthSignedErrorAverage(summary, character(0))
  over <- avg[avg$direction == "over", ]
  under <- avg[avg$direction == "under", ]
  expect_equal(over$prop, 0.75)
  expect_equal(over$median, 1.5)
  # Groups with no under-estimates do not dilute the under-estimate size.
  expect_equal(under$prop, 0.25)
  expect_equal(under$median, -0.5)
})

test_that("signed error scale puts nothing gated and two-fold equally far from zero", {
  env <- .bandwidth_cell_plot_env()
  trans <- env$.simBandwidthSignedErrorTrans()
  x <- c(-1, -0.5, 0, 1, 3)
  expect_equal(trans$transform(x), c(-1, -0.5, 0, 1, 2))
  expect_equal(trans$inverse(trans$transform(x)), x)
  # Values below -100% (from axis expansion) stay linear without warnings.
  expect_no_warning(expect_equal(trans$transform(c(-2, 1)), c(-2, 1)))
  expect_equal(trans$breaks(c(-1, 2.5)), c(-1, -0.5, 0, 1, 3))
  # Small errors get ordinary breaks rather than only zero.
  expect_equal(trans$breaks(c(-0.02, 0.03)), pretty(c(-0.02, 0.03)))
  expect_identical(
    env$.simBandwidthSignedErrorLabel(c(-1, 0, 1)),
    c("-100%", "0%", "+100% (2x)")
  )
})

test_that("2b error summary keeps the median, 90th percentile and maximum", {
  env <- .bandwidth_cell_plot_env()
  env$bias_uns_results_raw <- tidyr::expand_grid(
    transformation = "gaussian", mean_pos_setting = "high",
    prob_response = 0.05, n_cell = 100, mismatch_label = "mean shift 0",
    mismatch_type = "mean_shift", mismatch_val = 0, bw = 0.1,
    bias_uns_basis = "bandwidth", bias_uns_multiplier = c(0, 1),
    rel_error = c(-0.1, 0.05, 0.2, 0.4, NA)
  ) |>
    dplyr::mutate(abs_rel_error = abs(rel_error))
  eval(.bandwidth_cell_plot_chunk(
    "2b-sim-bias_uns-freq_bs.qmd", "bias-uns-error-summary"
  ), env)
  abs_error <- env$bias_uns_abs_error
  expect_equal(nrow(abs_error), 2L)
  expect_equal(abs_error$median_abs_rel_error, c(0.15, 0.15))
  expect_equal(
    abs_error$q90_abs_rel_error,
    rep(stats::quantile(c(0.1, 0.05, 0.2, 0.4), 0.9, names = FALSE), 2)
  )
  expect_equal(abs_error$max_abs_rel_error, c(0.4, 0.4))
  signed <- env$bias_uns_signed_error
  expect_equal(nrow(signed), 4L)
  expect_equal(signed$max[signed$direction == "over"], c(0.4, 0.4))
  expect_equal(signed$max[signed$direction == "under"], c(-0.1, -0.1))
})

test_that("2b leaves negative-width results out of figures and tables", {
  env <- .bandwidth_cell_plot_env()
  results <- tibble::tibble(
    bias_uns_basis = c("bandwidth", "negative_width", "bandwidth"),
    bias_uns_multiplier = c(0, 0.1, 1)
  )
  env$bias_uns_results_raw <- results
  env$bias_uns_results_summary <- results
  eval(.bandwidth_cell_plot_chunk(
    "2b-sim-bias_uns-freq_bs.qmd", "bias-uns-hide-negative-width"
  ), env)
  expect_identical(env$bias_uns_results_raw$bias_uns_multiplier, c(0, 1))
  expect_identical(env$bias_uns_results_summary$bias_uns_multiplier, c(0, 1))
})

test_that("signed-error plots draw dashed lines whose weight varies", {
  env <- .bandwidth_cell_plot_env()
  sides <- tidyr::expand_grid(
    bw = c(0.1, 0.2, 0.5), direction = c("over", "under")
  ) |>
    dplyr::mutate(
      prop = ifelse(direction == "over", 1, 0) + c(0.2, 0.8, 0.4, 0.6, 0.7, 0.3) *
        ifelse(direction == "over", -1, 1),
      sign = ifelse(direction == "over", 1, -1),
      median = sign * 0.05 * seq_len(dplyr::n()),
      q90 = 2 * median, q95 = 3 * median, max = 4 * median
    ) |>
    dplyr::select(-sign)

  global <- env$.simBandwidthGlobalSignedErrorPlot(
    dplyr::mutate(sides, transformation = "gaussian")
  )
  expect_no_error(ggplot2::ggplotGrob(global))

  bias <- env$.simBandwidthBiasSignedErrorPlot(
    sides |>
      dplyr::mutate(
        bias_uns_multiplier = bw, bw = 0.1, mismatch_label = "mean shift 0",
        bias_uns_basis = rep(c("bandwidth", "bandwidth", "negative_width"), each = 2)
      ),
    title = "test"
  )
  expect_no_error(ggplot2::ggplotGrob(bias))
})
