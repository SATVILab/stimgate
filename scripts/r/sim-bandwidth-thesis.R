# Thesis figures for the fixed-bandwidth simulation (Analysis 2a): StimGate's
# background-subtracted frequency against the bandwidth, as median lines over
# 10th-90th percentile bands (`.simCompareThesisBandPlot()`), for a focused
# set of scenarios. Needs analysis-plot-style.R and sim-compare-thesis.R.

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

# Per-scenario percentiles of the estimated background-subtracted frequency
# in percent over the valid estimates of the focused scenarios' samples
# (freq_q10, freq_median, freq_q90, with n_sample and n_valid).
.simBandwidthThesisSummary <- function(raw, scenarios = .simBandwidthThesisScenarios) {
  raw <- .simBandwidthThesisGrid(raw, scenarios)
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
      .groups = "drop"
    )
}

# Estimated frequency (percent, truth as a grey line; y from 0 to at least
# `truth_upper` times the truth, further only if the data reach higher)
# against the bandwidth (log scale, each transformation's own range).
# Transformations are columns (top strips); rows combine the number of cells
# and the separation in one strip on the right, as Analysis 6's figures label
# scenarios (e.g. "5,000 cells," over "low separation"), so the low and high
# separation of each cell count sit on adjacent rows.
.simBandwidthThesisPlot <- function(summary, band_lines = c("dotted", "thin", "none"),
                                    truth_upper = 1.5) {
  cells <- sort(unique(summary$n_cell))
  sep <- intersect(c("low", "high"), as.character(summary$mean_pos_setting))
  d <- data.frame(
    method = "stimgate", transformation = summary$transformation,
    cells = factor(summary$n_cell, cells, paste(scales::comma(cells, accuracy = 1), "cells")),
    separation = factor(as.character(summary$mean_pos_setting), sep,
      paste(sep, "separation")),
    x = summary$bw,
    lo = summary$freq_q10, mid = summary$freq_median, hi = summary$freq_q90
  )
  truth <- dplyr::distinct(data.frame(cells = d$cells, separation = d$separation,
    transformation = d$transformation, y = 100 * summary$prob_response))
  .simCompareThesisBandPlot(d, rows = c("cells", "separation"),
    x_scale = ggplot2::scale_x_log10(labels = .analysis_label_number,
      breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 1.5),
      expand = ggplot2::expansion(mult = 0.04)),
    x_lab = "Bandwidth", measure = "frequency", band_lines = band_lines, truth = truth,
    truth_upper = truth_upper, free_y = "row", row_sep = ",\n", single_dodge = 0.05)
}
