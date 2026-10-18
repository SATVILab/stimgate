# Thesis figures for the fixed-bandwidth simulation (Analysis 2a): StimGate's
# F1 and background-subtracted frequency against the bandwidth, as median
# lines over 10th-90th percentile bands (`.simCompareThesisBandPlot()`), for
# a focused set of scenarios. Needs analysis-plot-style.R and
# sim-compare-thesis.R.

# The focused scenarios, as in Analysis 6: response probability 0.2%, low and
# high separation, 5,000 and 100,000 cells, every transformation, at the
# low unstimulated-sample bias of 2a's other summary figures.
.simBandwidthThesisScenarios <- list(
  prob_response = 0.002, n_cell = c(5e3, 1e5), mean_pos_setting = c("low", "high"),
  bias_uns_setting = "low", condition_perturbation_sd = 0
)

# Rows of a 2a grid or result table in the focused scenarios.
.simBandwidthThesisGrid <- function(tbl, scenarios = .simBandwidthThesisScenarios) {
  keep <- rep(TRUE, nrow(tbl))
  for (col in names(scenarios)) {
    keep <- keep & as.character(tbl[[col]]) %in% as.character(scenarios[[col]])
  }
  tbl[keep, , drop = FALSE]
}

.simBandwidthThesisCountCols <- c("nTruePos", "nFalsePos", "nFalseNeg", "nTrueNeg")

# Per-scenario percentiles over the samples of the focused scenarios: the
# estimated background-subtracted frequency in percent over valid estimates
# (freq_q10, freq_median, freq_q90, with n_valid) and, when the results carry
# confusion-matrix counts, F1 = 2TP / (2TP + FP + FN) over the samples where
# it is defined (f1_q10, f1_median, f1_q90, with n_f1_defined; NA otherwise).
.simBandwidthThesisSummary <- function(raw, scenarios = .simBandwidthThesisScenarios) {
  raw <- .simBandwidthThesisGrid(raw, scenarios)
  has_counts <- all(.simBandwidthThesisCountCols %in% names(raw))
  raw$f1 <- if (has_counts) {
    den <- 2 * raw$nTruePos + raw$nFalsePos + raw$nFalseNeg
    ifelse(!is.na(den) & den > 0, 2 * raw$nTruePos / den, NA_real_)
  } else {
    NA_real_
  }
  pct <- function(x, p) {
    x <- x[is.finite(x)]
    if (length(x)) stats::quantile(x, p, names = FALSE, type = 7) else NA_real_
  }
  raw |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(
      "transformation", "mean_pos_setting", "n_cell", "prob_response", "bw"
    )))) |>
    dplyr::summarise(
      n_sample = dplyr::n(),
      n_valid = sum(.data$valid_estimate),
      freq_q10 = 100 * pct(.data$propRespEst[.data$valid_estimate], 0.1),
      freq_median = 100 * pct(.data$propRespEst[.data$valid_estimate], 0.5),
      freq_q90 = 100 * pct(.data$propRespEst[.data$valid_estimate], 0.9),
      n_f1_defined = sum(is.finite(.data$f1)),
      f1_q10 = pct(.data$f1, 0.1), f1_median = pct(.data$f1, 0.5), f1_q90 = pct(.data$f1, 0.9),
      .groups = "drop"
    )
}

# F1 or estimated frequency (percent, truth as a grey line) against the
# bandwidth (log scale, each transformation's own range), columns for
# transformation and rows for number of cells and then separation, so the low
# and high separation of each cell count sit on adjacent rows.
.simBandwidthThesisPlot <- function(summary, measure = c("frequency", "f1"),
                                    band_lines = c("dotted", "thin", "none")) {
  measure <- match.arg(measure)
  stat <- if (identical(measure, "f1")) "f1" else "freq"
  cells <- sort(unique(summary$n_cell))
  sep <- intersect(c("low", "high"), as.character(summary$mean_pos_setting))
  d <- data.frame(
    method = "stimgate", transformation = summary$transformation,
    cells = factor(summary$n_cell, cells, paste(scales::comma(cells, accuracy = 1), "cells")),
    separation = factor(as.character(summary$mean_pos_setting), sep,
      paste(tools::toTitleCase(sep), "separation")),
    x = summary$bw,
    lo = summary[[paste0(stat, "_q10")]], mid = summary[[paste0(stat, "_median")]],
    hi = summary[[paste0(stat, "_q90")]]
  )
  truth <- if (identical(measure, "frequency")) {
    dplyr::distinct(data.frame(cells = d$cells, separation = d$separation,
      transformation = d$transformation, y = 100 * summary$prob_response))
  }
  .simCompareThesisBandPlot(d, rows = c("cells", "separation"),
    x_scale = ggplot2::scale_x_log10(labels = .analysis_label_number,
      breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 1.5),
      expand = ggplot2::expansion(mult = 0.04)),
    x_lab = "Bandwidth", measure = measure, band_lines = band_lines, truth = truth,
    single_dodge = 0.05)
}
