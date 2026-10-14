# Analysis 7b (`7b-sim-compare-null.qmd`): Analysis 7's methods when nothing
# responds. Both tubes share an absolute background share of positive cells
# (`prob_background`, possibly zero) and `prob_response` is zero, so the true
# background-subtracted frequency is zero and every estimate is an error.
# Simulation, staging and promotion reuse Analysis 7's runner
# (`sim-compare-freq_bs.R`).

# Background levels: none, low (0.01%) and high-ish (1%).
.simCompareNullBackgrounds <- c(0, 1e-4, 1e-2)

.simCompareNullBackgroundLabel <- function(x) {
  ifelse(x == 0, "No positive cells in either tube",
    paste0(.analysis_label_percent(x), " positive cells in both tubes"))
}

# Seeds per biological scenario (every grid column except method settings),
# assigned on the full grid before dev/quick filters or chunking.
.simCompareNullSeedGrid <- function(grid, simulation_seed) {
  grid |>
    dplyr::group_by(
      .data$transformation, .data$mean_pos_setting, .data$mean_pos,
      .data$prob_background, .data$n_cell, .data$sample_perturbation_sd,
      .data$condition_perturbation_sd, .data$cluster_perturbation_sd,
      .data$n_cell_uns_relative_to_stim
    ) |>
    dplyr::mutate(
      base_scenario_id = dplyr::cur_group_id(),
      sim_seed = as.integer(simulation_seed + .data$base_scenario_id - 1L)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(sim_id = dplyr::row_number())
}

# Errors per sample in percentage points: estimate minus the simulated
# (labelled) net frequency, which is exactly zero with `probExact = TRUE`.
# Comparator runtime errors keep a missing error rather than a zero.
.simCompareNullErrors <- function(data, methods = .simCompareMethods) {
  data |>
    dplyr::filter(.data$method %in% methods) |>
    dplyr::mutate(
      error_pp = 100 * (.data$propRespEst - .data$propRespTruth),
      abs_error_pp = abs(.data$error_pp)
    )
}

# Pooled sample-level error statistics per scenario and method: median and
# 10th/90th percentiles of the signed error, median and 90th percentile of the
# absolute error (percentage points), and the share of samples with a
# positive estimate. With `mcse = TRUE`, 95% intervals come from resampling
# whole datasets (`iter`), sharing each scenario's draws across methods.
.simCompareNullSummary <- function(errors, group_cols, mcse = FALSE) {
  quantile_at <- function(p) {
    function(v) {
    v <- v[is.finite(v)]
    if (length(v)) stats::quantile(v, p, names = FALSE) else NA_real_
  }
  }
  share_positive <- function(v) {
    v <- v[is.finite(v)]
    if (length(v)) mean(v > 0) else NA_real_
  }
  errors |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(rows, key) {
      family <- unique(paste0("sim_seed:", rows$sim_seed))
      if (length(family) != 1L) stop("A null scenario must have one bootstrap family.")
      unit <- as.character(rows$iter)
      stats <- list(
        error_p10 = list(rows$error_pp, quantile_at(0.1)),
        error_median = list(rows$error_pp, quantile_at(0.5)),
        error_p90 = list(rows$error_pp, quantile_at(0.9)),
        abs_error_median = list(rows$abs_error_pp, quantile_at(0.5)),
        abs_error_p90 = list(rows$abs_error_pp, quantile_at(0.9)),
        share_positive = list(rows$error_pp, share_positive)
      )
      out <- tibble::tibble(
        n_sample = nrow(rows),
        n_dataset = dplyr::n_distinct(rows$iter),
        n_finite = sum(is.finite(rows$error_pp))
      )
      for (nm in names(stats)) {
        cols <- .analysis_mcse_pooled_cols(stats[[nm]][[1]], unit, stats[[nm]][[2]],
          nm, family, mcse = mcse)
        out <- dplyr::bind_cols(out, cols[!startsWith(names(cols), ".boot_")])
      }
      out
    }) |>
    dplyr::ungroup()
}

# How StimGate's gate was set per scenario: samples whose gate came from local
# FDR directly, from a cluster gate, or from a fallback (no response found).
.simCompareNullProvenance <- function(errors, group_cols) {
  errors |>
    dplyr::filter(.data$method == "stimgate") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_sample = dplyr::n(),
      n_direct = sum(.data$locGeneratedDirect %in% TRUE),
      n_shared = sum(.data$locGenerated %in% TRUE & !(.data$locGeneratedDirect %in% TRUE)),
      n_fallback = sum(!(.data$locGenerated %in% TRUE)),
      .groups = "drop"
    )
}

# One statistic against cell count, a line per method and a panel per
# transformation and (with `by_separation`) separation. Error bars (Monte Carlo intervals) only when
# the summary has them.
.simCompareNullPlot <- function(summary, stat, y_lab, zero_line = TRUE,
                               by_separation = TRUE, display_cap = Inf,
                                reference = NULL) {
  d <- summary |>
    dplyr::mutate(
      method = .analysis_method_factor(.data$method),
      value = .data[[stat]],
      lower = .data[[paste0(stat, "_lower")]],
      upper = .data[[paste0(stat, "_upper")]],
      transformation = .analysis_trans_factor(.data$transformation),
      separation = paste0(.data$mean_pos_setting, " separation")
    ) |>
    dplyr::filter(is.finite(.data$value))
  d$n_cell <- .analysis_cell_x(d$n_cell)
  d$capped <- d$value > display_cap
  d$cap_marker <- ifelse(d$capped, display_cap, NA_real_)
  d$value <- pmin(d$value, display_cap)
  d$lower <- pmin(d$lower, display_cap)
  d$upper <- pmin(d$upper, display_cap)
  dodge <- ggplot2::position_dodge(width = if (is.factor(d$n_cell)) 0.4 else 0.12)
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data$n_cell, y = .data$value, colour = .data$method,
    shape = .data$method, linetype = .data$method, group = .data$method
  ))
  if (zero_line) p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey50")
  if (length(reference)) {
    p <- p + ggplot2::geom_hline(
      yintercept = reference, colour = "grey50", linetype = "dotted"
    )
  }
  p +
    .analysis_mcse_layer(ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      width = 0, position = dodge, linetype = "solid", na.rm = TRUE
    )) +
    ggplot2::geom_line(position = dodge) +
    ggplot2::geom_point(position = dodge, size = 2) +
    .analysis_y_floor(if (zero_line) c(-0.01, 0.01) else c(0, 0.01)) +
    ggplot2::geom_point(ggplot2::aes(y = .data$cap_marker),
      shape = 24, fill = "white", size = 3, position = dodge, show.legend = FALSE, na.rm = TRUE
    ) +
    .analysis_scale_cells(d$n_cell) +
    .analysis_scale_method(c("colour", "shape", "linetype")) +
    (if (by_separation) {
      ggplot2::facet_wrap(~ transformation + separation, scales = "free_y", ncol = 2)
    } else {
      ggplot2::facet_wrap(~transformation, scales = "free_y", ncol = 1)
    }) +
    ggplot2::labs(x = "Parent cells per tube", y = y_lab,
      caption = if (is.finite(display_cap)) {
        paste0(
          "Hollow upward triangles: summaries above display cap ",
          .analysis_label_number(display_cap), ". Intervals capped for display only."
        )
      } else {
        NULL
      }
    ) +
    .analysis_theme()
}
# FPR curves: split background and separation in the QMD; facets keep every
# transformation/cell-count scenario distinct, with a common probability range.
.simCompareNullDetectionPlot <- function(summary) {
  d <- .analysis_method_data(summary)
  d$transformation <- .analysis_trans_factor(d$transformation)
  d$cell_label <- factor(.analysis_label_cells(d$n_cell),
    levels = .analysis_label_cells(sort(unique(d$n_cell)))
  )
  ggplot2::ggplot(d, ggplot2::aes(
    x = alpha, y = detection_rate,
    colour = method, shape = method, linetype = method, group = method
  )) +
    ggplot2::geom_hline(yintercept = 0.05, colour = "grey60", linetype = "dotted") +
    ggplot2::geom_vline(xintercept = 0.05, colour = "grey60", linetype = "dotted") +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    .analysis_mcse_errorbar(d, "detection_rate_lower", "detection_rate_upper") +
    ggplot2::scale_x_log10(breaks = c(0.001, 0.01, 0.05, 0.2), labels = .analysis_label_number) +
    ggplot2::scale_y_continuous(labels = .analysis_label_percent) +
    .analysis_y_floor(c(0, 1)) +
    ggplot2::facet_wrap(~ transformation + cell_label, ncol = 2) +
    .analysis_scale_method(c("colour", "shape", "linetype")) +
    .analysis_theme() +
    ggplot2::labs(x = "One-sided Fisher p-value cutoff", y = "Empirical false-positive rate")
}

# Explicit rectangular completion: absent scenarios are NA tiles, never zeros.
# Separate separations to keep the method x transformation grid readable.
.simCompareNullDetectionHeatmap <- function(summary) {
  d <- summary |>
    dplyr::filter(abs(.data$alpha - 0.05) < 1e-12) |>
    dplyr::select(method, transformation, n_cell, prob_background, detection_rate) |>
    dplyr::mutate(method = as.character(.data$method)) |>
    tidyr::complete(method, transformation, n_cell,
      prob_background = .simCompareNullBackgrounds
    )
  d <- .analysis_method_data(d)
  d$transformation <- .analysis_trans_factor(d$transformation)
  d$n_cell <- factor(d$n_cell, sort(unique(d$n_cell)))
  d$background <- factor(d$prob_background, sort(unique(d$prob_background)))
  ggplot2::ggplot(d, ggplot2::aes(x = n_cell, y = background, fill = detection_rate)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(detection_rate), "NA",
      .analysis_label_percent(detection_rate)
    ), colour = ifelse(
      !is.na(detection_rate) & detection_rate < 0.5, "white", "#333333"
    )), size = 2.5) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_viridis_c(
      limits = c(0, 1), na.value = "grey85",
      labels = .analysis_label_percent, breaks = seq(0, 1, 0.25), name = "False-positive rate",
      guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(5, "cm"))
    ) +
    ggplot2::scale_x_discrete(labels = .analysis_label_cells) +
    ggplot2::scale_y_discrete(labels = function(x) .analysis_label_percent(as.numeric(x))) +
    ggplot2::facet_grid(method ~ transformation,
      labeller = ggplot2::labeller(method = .analysis_method_labels)
    ) +
    .analysis_theme(grid = "none") +
    ggplot2::labs(
      x = "Parent cells per tube", y = "Background-positive frequency",
      caption = "Grey NA tiles: absent scenarios or no valid gates. Common colour scale: 0–100%."
    )
}

.simCompareDetectionRocPlot <- function(curves) {
  d <- .analysis_method_data(curves)
  d$transformation <- .analysis_trans_factor(d$transformation)
  d$cell_label <- factor(
    .analysis_label_cells(d$n_cell),
    .analysis_label_cells(sort(unique(d$n_cell)))
  )
  marked <- d[vapply(d$alpha, function(x) any(abs(x - c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2)) < 1e-12), logical(1)), ]
  ggplot2::ggplot(d, ggplot2::aes(
    x = detection_rate_null, y = detection_rate_response,
    colour = method, shape = method, linetype = method, group = method
  )) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey65", linetype = "dotted") +
    ggplot2::geom_vline(xintercept = 0.05, colour = "grey60", linetype = "dashed") +
    # Paths follow the decision threshold, even where coordinates tie.
    ggplot2::geom_path() +
    ggplot2::geom_point(data = marked, size = 1.5) +
    .analysis_mcse_errorbar(marked, "detection_rate_lower_response", "detection_rate_upper_response") +
    .analysis_mcse_layer(ggplot2::geom_errorbar(
      ggplot2::aes(xmin = detection_rate_lower_null, xmax = detection_rate_upper_null),
      data = marked, orientation = "y", width = 0, na.rm = TRUE
    )) +
    ggplot2::facet_wrap(~ transformation + cell_label, ncol = 2) +
    ggplot2::scale_x_continuous(labels = .analysis_label_percent) +
    ggplot2::scale_y_continuous(labels = .analysis_label_percent) +
    ggplot2::expand_limits(x = c(0, 1), y = c(0, 1)) +
    .analysis_scale_method(c("colour", "shape", "linetype")) +
    .analysis_theme() +
    ggplot2::labs(
      x = "Empirical false-positive rate (null)", y = "Detection sensitivity (response)",
      caption = "Paths: all cutoffs. Symbols/intervals: 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2."
    )
}
