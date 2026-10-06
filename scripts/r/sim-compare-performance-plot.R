# Separate occurrence and conditional-severity views of fixed-size dataset maxima.
# Source after sim-compare-freq_bs.R and sim-bandwidth-analysis-plot.R.
.simCompareDatasetMaxCoverage <- function(tbl) {
  counts <- c("n_dataset_total", "n_eligible", "n_incomplete", "n_affected",
    "n_bootstrap_finite_occurrence", "n_bootstrap_finite_severity",
    "bootstrap_coverage_occurrence", "bootstrap_coverage_severity",
    "interval_available_occurrence", "interval_available_severity")
  keys <- setdiff(names(tbl), c(names(tbl)[startsWith(names(tbl), ".boot_")],
    counts, "occurrence", "severity",
    paste0(rep(c("occurrence", "severity"), each = 3), c("_lower", "_upper", "_mcse"))))
  dplyr::select(tbl, dplyr::any_of(c(keys, counts)))
}

# Preserve group keys plus bootstrap validity/coverage diagnostics beside figures.
.simComparePerformanceCoverage <- function(tbl, keys) {
  columns <- grep("^(n_|bootstrap_coverage_|interval_available_)", names(tbl), value = TRUE)
  dplyr::distinct(dplyr::select(tbl, dplyr::any_of(c(keys, "direction", columns))))
}

.simCompareDatasetMaxPlot <- function(
    tbl, quantity = c("occurrence", "severity"), x = "n_cell",
    x_label = "Number of stimulated cells", x_log = TRUE,
    by_prob = TRUE, mcse = FALSE) {
  quantity <- match.arg(quantity)
  tbl$transformation <- .analysis_trans_factor(tbl$transformation)
  sign <- if (quantity == "severity") ifelse(tbl$direction == "over", 1, -1) else rep(1, nrow(tbl))
  tbl$value <- sign * tbl[[quantity]]
  lower <- tbl[[paste0(quantity, "_lower")]]
  upper <- tbl[[paste0(quantity, "_upper")]]
  tbl$lower <- ifelse(sign > 0, lower, -upper)
  tbl$upper <- ifelse(sign > 0, upper, -lower)
  capped <- quantity == "severity" && .simBandwidthSignedErrorIsCapped(tbl$value)
  tbl$value_shown <- if (quantity == "severity") .simBandwidthSignedErrorSquish(tbl$value) else tbl$value
  tbl$direction <- factor(tbl$direction, c("over", "under"), c("Over-estimates", "Under-estimates"))
  p <- ggplot2::ggplot(tbl, ggplot2::aes(
    x = .data[[x]], y = value_shown, colour = method, shape = method, linetype = method, group = method
  )) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.75) +
    ggplot2::geom_point(size = 1.5, alpha = 0.75) +
    ggplot2::facet_wrap(
      if (by_prob) ggplot2::vars(direction, transformation, prob_response) else
        ggplot2::vars(direction, transformation),
      ncol = min(3L, length(unique(tbl$transformation)) *
        if (by_prob) length(unique(tbl$prob_response)) else 1L),
      scales = if (quantity == "occurrence") "fixed" else "free_y",
      labeller = ggplot2::labeller(prob_response = .analysis_labeller_percent("Response probability: "))
    ) +
    (if (x_log) ggplot2::scale_x_log10(
      breaks = sort(unique(tbl[[x]])), labels = .analysis_label_number,
      guide = ggplot2::guide_axis(angle = 45)
    ) else ggplot2::scale_x_continuous(labels = .analysis_label_number)) +
    .analysis_scale_method(c("colour", "shape", "linetype")) + .analysis_theme() +
    ggplot2::labs(x = x_label, colour = "Method")
  if (quantity == "occurrence") {
    p <- p + ggplot2::scale_y_continuous(labels = .analysis_label_percent) +
      .analysis_y_floor(c(0, 1)) +
      ggplot2::labs(y = "Eligible datasets with at least one directional error")
    if (isTRUE(mcse)) p <- p + .analysis_mcse_errorbar(tbl)
  } else {
    p <- p + .simBandwidthSignedErrorLayers(
      "Conditional mean dataset maximum (signed relative error)", capped = capped
    )
    if (isTRUE(mcse)) p <- p + .simBandwidthSignedErrorBars(tbl)
  }
  span <- function(col) paste(range(tbl[[col]], na.rm = TRUE), collapse = "–")
  p + ggplot2::labs(caption = paste0(
    "Dataset counts per setting/method: eligible ", span("n_eligible"),
    "; incomplete ", span("n_incomplete"), "; affected ", span("n_affected"), ".\n",
    if (quantity == "occurrence") "Occurrence intervals require ≥5 eligible datasets and ≥95% finite bootstrap draws."
    else "Conditional intervals require ≥5 eligible and ≥5 affected datasets and ≥95% finite bootstrap draws."
  ))
}
