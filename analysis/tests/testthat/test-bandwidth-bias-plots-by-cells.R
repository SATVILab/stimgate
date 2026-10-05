.bandwidth_cell_plot_chunk <- function(document, label) {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  lines <- readLines(file.path(root, "analysis", document), warn = FALSE)
  start <- which(lines == paste0("#| label: ", label))
  testthat::expect_length(start, 1L)
  end <- which(lines == "```" & seq_along(lines) > start)[1L]
  parse(text = lines[seq.int(start + 1L, end - 1L)])
}

# Evaluate a plot chunk, returning the Markdown it writes (headings).
.bandwidth_cell_plot_eval <- function(code, env) {
  paste(utils::capture.output(eval(code, env)), collapse = "\n")
}

.bandwidth_cell_plot_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c("analysis-plot-style.R", "sim-bandwidth-analysis-plot.R")) {
    source(file.path(testthat::test_path(), "../../../scripts/r", fn), local = env)
  }
  env$root_dir <- tempfile("cell-plots-")
  env$analysis_key <- "bias_uns"
  env$run_plots <- TRUE
  env$saved <- list()
  env$printed <- list()
  env$fig_key <- "2a-test"
  # Coverage-table behaviour is tested below independently of these plot fixtures.
  env$bw_tbl_results_summary <- tibble::tibble()
  env$bias_uns_results_summary <- tibble::tibble()
  env$.simBandwidthPrintCoverage <- function(plot, summary) invisible(NULL)
  env$.analysis_fig_dir <- function(path_parts, path_root = NULL, create = TRUE) {
    file.path(path_root, "output", "fig", paste(path_parts, collapse = "/"))
  }
  env$.analysis_cache_dir <- function(parts, path_root) {
    file.path(path_root, paste(parts, collapse = "/"))
  }
  env$.analysis_save_fig <- function(plot, path, height = 12, width = 16,
                                     allow_tall = FALSE) {
    env$saved[[length(env$saved) + 1L]] <- list(
      path = path, plot = plot, height = height
    )
    invisible(path)
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
    "2a-sim-bw-freq_bs-global.qmd", "relative-error-by-n-cell"
  )
  md <- .bandwidth_cell_plot_eval(code, env)
  expect_match(md, "#### Mean position: low", fixed = TRUE)
  expect_match(md, "##### Cells: 1,000", fixed = TRUE)
  expect_length(env$saved, 4L)
  expect_identical(env$printed, lapply(env$saved, `[[`, "plot"))
  paths <- vapply(env$saved, `[[`, character(1), "path")
  expect_equal(length(unique(paths)), 4L)
  for (saved in env$saved) {
    data <- saved$plot$data
    expect_length(unique(data$n_cell), 1L)
    expect_length(unique(data$mean_pos_setting), 1L)
    # Transformations are labelled and ordered Gaussian, Skew, Gamma.
    expect_identical(levels(data$transformation)[1:3], c("Gaussian", "Skew", "Gamma"))
    expect_setequal(as.character(unique(data$transformation)), c("Gaussian", "Gamma"))
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
  expect_identical(.bandwidth_cell_plot_eval(code, env), "")
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
    "2b-sim-bias_uns-freq_bs.qmd", "relative-error-averaged-n-cell"
  )
  cell_code <- .bandwidth_cell_plot_chunk(
    "2b-sim-bias_uns-freq_bs.qmd", "relative-error-by-n-cell"
  )
  md <- .bandwidth_cell_plot_eval(average_code, env)
  # Gaussian comes before Gamma; headings nest one level per loop.
  expect_lt(
    regexpr("#### Transformation: Gaussian", md, fixed = TRUE),
    regexpr("#### Transformation: Gamma", md, fixed = TRUE)
  )
  expect_match(md, "##### Mean position: low", fixed = TRUE)
  expect_match(md, "###### Response probability: 1%", fixed = TRUE)
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
  md <- .bandwidth_cell_plot_eval(cell_code, env)
  expect_match(md, "###### Response probability: 10%; cells: 1,000", fixed = TRUE)
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
    # Bandwidth is the only colour legend; there is no bias-scale line type.
    expect_identical(rlang::as_label(saved$plot$mapping$colour), "bw_lab")
    expect_null(saved$plot$mapping$linetype)
    expect_null(saved$plot$labels$title)
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
    "2a-sim-bw-freq_bs-global.qmd", "signed-error-by-n-cell"
  )
  md <- .bandwidth_cell_plot_eval(code, env)
  expect_match(md, "##### Cells: 100", fixed = TRUE)
  expect_length(env$saved, 4L)
  expect_identical(env$printed, lapply(env$saved, `[[`, "plot"))
  paths <- vapply(env$saved, `[[`, character(1), "path")
  expect_equal(length(unique(paths)), 4L)
  for (saved in env$saved) {
    data <- saved$plot$data
    expect_length(unique(data$n_cell), 1L)
    expect_length(unique(data$mean_pos_setting), 1L)
    expect_setequal(as.character(unique(data$transformation)), c("Gaussian", "Gamma"))
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
    "2b-sim-bias_uns-freq_bs.qmd", "signed-error-averaged-n-cell"
  )
  cell_code <- .bandwidth_cell_plot_chunk(
    "2b-sim-bias_uns-freq_bs.qmd", "signed-error-by-n-cell"
  )
  .bandwidth_cell_plot_eval(average_code, env)
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
  .bandwidth_cell_plot_eval(cell_code, env)
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
    expect_null(saved$plot$mapping$linetype)
    # Errors above +1500% are drawn at the cap; `value` keeps the actual error.
    expect_equal(data$value_shown, pmin(data$value, 15))
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

test_that("signed error scale puts a zero estimate and two-fold equally far from zero", {
  env <- .bandwidth_cell_plot_env()
  trans <- env$.simBandwidthSignedErrorTrans()
  x <- c(-2, -1.5, -1, -0.5, 0, 1, 3)
  expect_equal(trans$transform(x), c(-2, -1.5, -1, -0.5, 0, 1, 2))
  expect_equal(trans$inverse(trans$transform(x)), x)
  expect_equal(trans$domain, c(-Inf, Inf))
  expect_true(all(c(-2, -1) %in% trans$breaks(c(-2, 3))))
  expect_identical(env$.simBandwidthSignedErrorLabel(-1.5), "-150%")
  # Negative background-subtracted response estimates stay linear without warnings.
  expect_no_warning(expect_equal(trans$transform(c(-2, 1)), c(-2, 1)))
  expect_equal(trans$breaks(c(-1, 2.5)), c(-1, -0.5, 0, 1, 3))
  # Small errors get ordinary breaks rather than only zero.
  expect_equal(trans$breaks(c(-0.02, 0.03)), pretty(c(-0.02, 0.03)))
  expect_identical(
    env$.simBandwidthSignedErrorLabel(c(-1, -0.5, 0, 0.5, 1, 3, 7)),
    c(
      "-100% (0x)", "-50%", "0%", "+50%", "+100% (2x)", "+300% (4x)",
      "+700% (8x)"
    )
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
  # Negative-width results and bias multiplier 0 are both left out.
  expect_identical(env$bias_uns_results_raw$bias_uns_multiplier, 1)
  expect_identical(env$bias_uns_results_summary$bias_uns_multiplier, 1)
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
  # Small errors still show the full -100% to +100% range.
  y_range <- ggplot2::layer_scales(global)$y$range$range
  expect_lte(y_range[[1]], -1)
  expect_gte(y_range[[2]], 1)

  by_prob <- env$.simBandwidthGlobalSignedErrorPlot(
    dplyr::bind_rows(
      dplyr::mutate(sides, prob_response = 0.001),
      dplyr::mutate(sides, prob_response = 0.01)
    ) |>
      dplyr::mutate(transformation = "gaussian"),
    by_prob = TRUE
  )
  expect_named(by_prob$facet$params$rows, "prob_response")
  expect_named(by_prob$facet$params$cols, "transformation")
  expect_no_error(ggplot2::ggplotGrob(by_prob))

  bias <- env$.simBandwidthBiasSignedErrorPlot(
    sides |>
      dplyr::mutate(
        bias_uns_multiplier = bw, bw = 0.1, mismatch_label = "mean shift 0",
        bias_uns_basis = rep(c("bandwidth", "bandwidth", "negative_width"), each = 2)
      ),
    title = "test"
  )
  expect_no_error(ggplot2::ggplotGrob(bias))
  # Titles are not drawn; headings in the QMD carry that information.
  expect_null(bias$labels$title)
})

test_that("signed relative errors above +1500% are drawn at a labelled cap", {
  env <- .bandwidth_cell_plot_env()
  expect_equal(env$.simBandwidthSignedErrorCap, 15)
  expect_equal(
    env$.simBandwidthSignedErrorSquish(c(-1, 3, 15, 40, NA)),
    c(-1, 3, 15, 15, NA)
  )
  expect_true(env$.simBandwidthSignedErrorIsCapped(c(1, 40, NA)))
  expect_false(env$.simBandwidthSignedErrorIsCapped(c(1, 15, NA)))
  # Without a cap the +1500% tick is a plain doubling; with it, "at least".
  expect_identical(
    env$.simBandwidthSignedErrorLabel(c(-1, 1, 7, 15)),
    c("-100% (0x)", "+100% (2x)", "+700% (8x)", "+1500% (16x)")
  )
  expect_identical(
    env$.simBandwidthSignedErrorLabel(c(-1, 0.5, 1, 15), cap = 15),
    c("-100% (0x)", "+50%", "+100% (2x)", "\u2265 +1500% (16x)")
  )
  # The cap is one of the breaks whenever the data reach it.
  expect_true(15 %in% env$.simBandwidthSignedErrorTrans()$breaks(c(-1, 15)))

  sides <- tidyr::expand_grid(
    bw = c(0.1, 0.2), direction = c("over", "under")
  ) |>
    dplyr::mutate(
      transformation = "gaussian",
      prop = 0.5,
      over = direction == "over",
      median = ifelse(over, 2, -0.5),
      q95 = ifelse(over, 10, -0.8),
      max = ifelse(over, 40, -0.9)
    ) |>
    dplyr::select(-over)
  global <- env$.simBandwidthGlobalSignedErrorPlot(sides)
  expect_equal(max(global$data$err_value), 40)
  expect_equal(max(global$data$err_value_shown), 15)
  expect_identical(rlang::as_label(global$mapping$y), "err_value_shown")
  y_labels <- global$scales$get_scales("y")$labels
  expect_identical(y_labels(c(1, 15)), c("+100% (2x)", "\u2265 +1500% (16x)"))
  expect_no_error(ggplot2::ggplotGrob(global))

  # Without any capped value, the top tick keeps its plain label.
  uncapped <- env$.simBandwidthGlobalSignedErrorPlot(
    dplyr::mutate(sides, max = pmin(max, 15))
  )
  expect_identical(
    uncapped$scales$get_scales("y")$labels(15), "+1500% (16x)"
  )

  bias <- env$.simBandwidthBiasSignedErrorPlot(
    sides |>
      dplyr::mutate(
        bias_uns_multiplier = bw, bw = 0.1, mismatch_label = "mean shift 0",
        bias_uns_basis = "bandwidth", q90 = q95
      )
  )
  expect_equal(max(bias$data$value_shown), 15)
  expect_identical(
    bias$scales$get_scales("y")$labels(15), "\u2265 +1500% (16x)"
  )
  # Over/under lines and points are slightly transparent.
  alphas <- vapply(bias$layers, function(l) {
    if (is.null(l$aes_params$alpha)) NA_real_ else l$aes_params$alpha
  }, numeric(1))
  expect_true(all(alphas[!is.na(alphas)] == 0.75))
  expect_gte(sum(!is.na(alphas)), 2L)
  expect_no_error(ggplot2::ggplotGrob(bias))
})


test_that("coverage companions preserve finite fallbacks and all failed samples", {
  env <- .bandwidth_cell_plot_env()
  on.exit(unlink(env$root_dir, recursive = TRUE), add = TRUE)
  summary <- tibble::tibble(
    transformation = "gaussian", mean_pos_setting = c("low", "low", "high"),
    bw = 0.1, n_cell = c(100, 1000, 100),
    n_sample = 4L, n_valid = c(3L, 0L, 4L), n_failed = c(1L, 4L, 0L),
    n_provenance = 4L, n_fallback = c(2L, 4L, 0L)
  )
  p <- ggplot2::ggplot(tibble::tibble(
    transformation = env$.analysis_trans_factor("gaussian"),
    mean_pos_setting = "low", bw = 0.1
  ))
  companion <- env$.simBandwidthCoverageForPlot(p, summary)
  expect_equal(companion$n_scenario, 2L)
  expect_equal(companion$n_scenario_valid, 1L)
  expect_equal(companion$n_sample, 8L)
  expect_equal(companion$n_valid, 3L)
  expect_equal(companion$failure_fraction, 5 / 8)
  expect_equal(companion$fallback_fraction, 6 / 8)
  # Per-cell figures retain only the corresponding expected scenario.
  p$data$n_cell <- factor(1000)
  companion <- env$.simBandwidthCoverageForPlot(p, summary)
  expect_equal(companion$n_scenario_valid, 0L)
  expect_equal(companion$failure_fraction, 1)
})

test_that("negative background-subtracted estimates survive signed plot building", {
  env <- .bandwidth_cell_plot_env()
  on.exit(unlink(env$root_dir, recursive = TRUE), add = TRUE)
  data <- tibble::tibble(
    transformation = "gaussian", bw = c(0.1, 0.2), direction = "under",
    prop = 1, median = -1.5, q95 = -2, max = -3
  )
  plot <- env$.simBandwidthGlobalSignedErrorPlot(data)
  built <- ggplot2::ggplot_build(plot)
  expect_true(any(vapply(built$data, function(layer) {
    "y" %in% names(layer) && any(layer$y < -1, na.rm = TRUE)
  }, logical(1))))
  summary <- env$.simBandwidthSignedErrorAverage(
    dplyr::bind_rows(data, dplyr::mutate(data, median = NA_real_)),
    c("transformation", "bw")
  )
  expect_equal(summary$n_scenario_median, c(1L, 1L))
  expect_equal(summary$n_scenario_max, c(2L, 2L))
})
