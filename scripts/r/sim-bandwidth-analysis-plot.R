format_bw_lab <- function(x, digits = 4) {
  x <- suppressWarnings(as.numeric(x))
  lab <- format(signif(x, digits), scientific = FALSE, trim = TRUE)
  decimal <- grepl(".", lab, fixed = TRUE)
  lab[decimal] <- sub("\\.?0+$", "", lab[decimal])
  lab[is.na(x)] <- NA_character_
  lab
}

format_bw_file <- function(x) {
  x <- format_bw_lab(x)
  x <- gsub("\\.", "p", x)
  gsub("[^A-Za-z0-9]+", "_", x)
}

safe_file_lab <- function(x) {
  x <- as.character(x)
  gsub("[^A-Za-z0-9]+", "_", x)
}

make_bw_colour_values <- function(bw_vec, base_col_vec = NULL) {
  if (is.null(base_col_vec) || length(base_col_vec) == 0L) {
    base_col_vec <- c(
      "#e0ecf4",
      "#bfd3e6",
      "#9ebcda",
      "#8c96c6",
      "#8c6bb1",
      "#88419d",
      "#810f7c"
    )
  }

  bw_num <- sort(unique(as.numeric(bw_vec)))
  bw_lab <- format_bw_lab(bw_num)
  n_bw <- length(bw_num)

  if (n_bw <= length(base_col_vec)) {
    col_vec <- base_col_vec[seq(
      length(base_col_vec) - n_bw + 1L,
      length(base_col_vec)
    )]
  } else {
    col_vec <- grDevices::colorRampPalette(base_col_vec)(n_bw)
  }

  stats::setNames(col_vec, bw_lab)
}

make_bw_linetype_scale <- function(bw_vec) {
  bw_num <- sort(unique(as.numeric(bw_vec)))
  bw_lab <- format_bw_lab(bw_num)
  n_bw <- length(bw_num)

  if (n_bw == 1L) {
    lty <- "solid"
  } else {
    dash_len <- round(seq(8, 2, length.out = n_bw - 1L))
    gap_len <- 4L

    lty <- c(
      "solid",
      paste0(
        as.character(as.hexmode(dash_len)),
        as.character(as.hexmode(gap_len))
      )
    )
  }

  list(
    levels = bw_lab,
    values = stats::setNames(lty, bw_lab)
  )
}

add_bw_labs <- function(.data) {
  if (!is.data.frame(.data)) {
    return(.data)
  }

  .data |>
    dplyr::mutate(
      bw_core_lab = format_bw_lab(.data$bw_core),
      bw_extra_lab = format_bw_lab(.data$bw_extra)
    )
}

# Pivot the wide statistic columns `stats` to long form. Monte Carlo interval
# columns (`<stat>_lower`, `<stat>_upper`; see `analysis-mcse.R`), when
# present, become `lower` and `upper` (NA for statistics without them).
.simBandwidthStatLongBounds <- function(tbl, stats, names_to, values_to) {
  bound_names <- as.vector(outer(stats, c("_lower", "_upper", "_mcse"), paste0))
  has_bounds <- any(bound_names %in% names(tbl))
  # pivot_longer() emits each input row's statistics in turn, so the bounds
  # follow the same row-major order.
  bound_vec <- function(suffix) {
    m <- do.call(cbind, lapply(stats, function(s) {
      col <- paste0(s, suffix)
      if (col %in% names(tbl)) {
        as.numeric(tbl[[col]])
      } else {
        rep(NA_real_, nrow(tbl))
      }
    }))
    as.vector(t(m))
  }
  long <- tbl |>
    dplyr::select(-dplyr::any_of(bound_names)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(stats),
      names_to = names_to,
      values_to = values_to
    )
  if (has_bounds) {
    long$lower <- bound_vec("_lower")
    long$upper <- bound_vec("_upper")
  }
  long
}

# Error statistics as rows of panels: `stat_cols` maps wide columns to labels.
.simBandwidthErrorStatLong <- function(tbl, stat_cols) {
  tbl |>
    .simBandwidthStatLongBounds(names(stat_cols), "statistic", "value") |>
    dplyr::filter(!is.na(.data$value)) |>
    dplyr::mutate(
      statistic = factor(
        unname(stat_cols[.data$statistic]),
        levels = unname(stat_cols)
      )
    )
}

# Bandwidth colours (the sequential purple ramp) with every bandwidth shown in
# the legend, which sits underneath the plot in one row.
.simBandwidthBwColourScale <- function(bw_vec) {
  col_vec <- make_bw_colour_values(bw_vec)
  ggplot2::scale_colour_manual(
    values = col_vec,
    breaks = names(col_vec),
    limits = names(col_vec),
    drop = FALSE,
    guide = ggplot2::guide_legend(nrow = 1, order = 1)
  )
}

# Bandwidth labels as a factor in increasing order, matching the colour scale.
.simBandwidthBwLabFactor <- function(bw) {
  factor(format_bw_lab(bw), levels = names(make_bw_colour_values(bw)))
}

# Bias curves share the same dimensions in all absolute-error views. Colour is
# the bandwidth; curves for different bias scales stay separate (`group`)
# although only the bandwidth rule is shown in Analysis 2b. `title` is ignored
# (figure titles go in headings) and is kept so older calls still work.
# `mcse`: draw the `<stat>_lower`/`<stat>_upper` Monte Carlo intervals.
.simBandwidthBiasRelativeErrorPlot <- function(
  tbl,
  title = NULL,
  y_label = "Absolute relative error",
  facet = ggplot2::facet_grid(statistic ~ mismatch_label, scales = "free_y"),
  stat_cols = c(
    median_abs_rel_error = "Median",
    q90_abs_rel_error = "90th percentile",
    max_abs_rel_error = "Maximum"
  ),
  mcse = FALSE
) {
  tbl <- .simBandwidthErrorStatLong(tbl, stat_cols) |>
    dplyr::mutate(bw_lab = .simBandwidthBwLabFactor(.data$bw))
  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = bias_uns_multiplier,
      y = value,
      colour = bw_lab,
      group = interaction(bw, bias_uns_basis)
    )
  ) +
    # Slight transparency shows overlapping lines.
    ggplot2::geom_line(alpha = 0.75) +
    ggplot2::geom_point(alpha = 0.75) +
    (if (isTRUE(mcse) && "lower" %in% names(tbl)) {
      .analysis_mcse_errorbar(tbl)
    }) +
    facet +
    .simBandwidthBwColourScale(tbl$bw) +
    ggplot2::scale_x_continuous(labels = .analysis_label_number) +
    .analysis_theme() +
    ggplot2::labs(
      x = "Bias multiplier",
      y = y_label,
      colour = "Bandwidth"
    )
}

# Signed-error version of the bias curves. `tbl` comes from
# `.simBandwidthSignedErrorSummary()`: over-estimates sit above zero and
# under-estimates below; line weight is each direction's share. Errors above
# +1500% are drawn at +1500% (`value_shown`); `value` keeps the actual error.
.simBandwidthBiasSignedErrorPlot <- function(
  tbl,
  title = NULL,
  y_label = "Relative error",
  facet = ggplot2::facet_grid(statistic ~ mismatch_label, scales = "free_y"),
  stat_cols = c(median = "Median", q90 = "90th percentile", max = "Maximum"),
  mcse = FALSE
) {
  tbl <- .simBandwidthErrorStatLong(tbl, stat_cols) |>
    dplyr::mutate(
      bw_lab = .simBandwidthBwLabFactor(.data$bw),
      value_shown = .simBandwidthSignedErrorSquish(.data$value)
    )
  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = bias_uns_multiplier,
      y = value_shown,
      colour = bw_lab,
      group = interaction(bw, bias_uns_basis, direction)
    )
  ) +
    .simBandwidthSignedErrorLayers(
      y_label,
      capped = .simBandwidthSignedErrorIsCapped(tbl$value)
    ) +
    .simBandwidthSignedErrorSegmentLayer(
      tbl,
      "bias_uns_multiplier",
      "value_shown",
      c(
        "statistic",
        "mismatch_label",
        "n_cell",
        "bw",
        "bias_uns_basis",
        "direction"
      ),
      alpha = 0.75
    ) +
    # Slight transparency shows overlapping lines; legend keys match.
    ggplot2::geom_point(size = 1, alpha = 0.75) +
    (if (isTRUE(mcse)) .simBandwidthSignedErrorBars(tbl)) +
    facet +
    .simBandwidthBwColourScale(tbl$bw) +
    ggplot2::scale_x_continuous(labels = .analysis_label_number) +
    .analysis_theme() +
    ggplot2::labs(x = "Bias multiplier", colour = "Bandwidth")
}

# ColorBrewer BrBG: teal for over-estimates, brown for under-estimates.
.simBandwidthSignedErrorColours <- c(
  over_median = "#80CDC1",
  over_q95 = "#35978F",
  over_max = "#01665E",
  under_median = "#DFC27D",
  under_q95 = "#BF812D",
  under_max = "#8C510A"
)

# 2a curves: colour is the direction (two halves of a diverging palette, so
# neither looks worse) and shade the statistic (darker = further from truth).
# With `by_prob`, rows of panels separate response probabilities. Errors above
# +1500% are drawn at +1500% (`err_value_shown`). `title` is ignored (figure
# titles go in headings) and is kept so older calls still work.
# `mcse`: draw the Monte Carlo intervals carried by the summary.
.simBandwidthGlobalSignedErrorPlot <- function(
  tbl,
  title = NULL,
  by_prob = FALSE,
  mcse = FALSE
) {
  bw_levels <- .analysis_label_number(sort(unique(tbl$bw)))
  tbl <- tbl |>
    .simBandwidthStatLongBounds(
      c("median", "q95", "max"), "err_type", "err_value"
    ) |>
    dplyr::filter(!is.na(.data$err_value)) |>
    dplyr::mutate(
      err_type = factor(err_type, levels = c("median", "q95", "max")),
      direction = factor(direction, levels = c("over", "under")),
      transformation = .analysis_trans_factor(.data$transformation),
      bw_fct = factor(.analysis_label_number(.data$bw), levels = bw_levels),
      err_value_shown = .simBandwidthSignedErrorSquish(.data$err_value),
      series = factor(
        paste0(direction, "_", err_type),
        levels = names(.simBandwidthSignedErrorColours)
      )
    )
  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = bw_fct,
      y = err_value_shown,
      colour = series,
      group = interaction(err_type, direction)
    )
  ) +
    .simBandwidthSignedErrorLayers(
      capped = .simBandwidthSignedErrorIsCapped(tbl$err_value)
    ) +
    .simBandwidthSignedErrorSegmentLayer(
      tbl,
      "bw_fct",
      "err_value_shown",
      c("transformation", "prob_response", "err_type", "direction"),
      alpha = 0.75
    ) +
    # Slight transparency shows overlapping lines; legend keys match.
    ggplot2::geom_point(size = 1, alpha = 0.75) +
    (if (isTRUE(mcse)) .simBandwidthSignedErrorBars(tbl, width = 0.3)) +
    (if (by_prob) {
      ggplot2::facet_grid(
        prob_response ~ transformation,
        scales = "free",
        labeller = ggplot2::labeller(
          prob_response = .analysis_labeller_percent("Response: ")
        )
      )
    } else {
      ggplot2::facet_wrap(~transformation, scales = "free")
    }) +
    .analysis_theme() +
    ggplot2::scale_colour_manual(
      values = .simBandwidthSignedErrorColours,
      labels = c(
        over_median = "Over: median",
        over_q95 = "Over: 95th percentile",
        over_max = "Over: maximum",
        under_median = "Under: median",
        under_q95 = "Under: 95th percentile",
        under_max = "Under: maximum"
      ),
      drop = FALSE,
      guide = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::labs(x = "Bandwidth", colour = NULL) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}

# Signed relative error, (estimate - truth) / truth, summarised separately for
# over- and under-estimates. `prop` is each side's share of finite errors; the
# other statistics are signed and use that side's errors only.
# With `mcse`, each percentile also gets `<stat>_lower`, `<stat>_upper` and
# `<stat>_mcse` (`analysis-mcse.R`; NA for the maximum). Without `unit`, these
# come from order statistics of that side's values, treated as independent.
# With `unit` (one value per error, e.g. the simulated dataset), the plotted
# percentile is unchanged and its MCSE is the between-unit spread of the
# per-unit percentile of that side. Bounds stay on their side of zero.
.simBandwidthSignedErrorSides <- function(rel_error, mcse = FALSE, unit = NULL) {
  keep <- is.finite(rel_error)
  rel_error <- rel_error[keep]
  if (!is.null(unit)) {
    unit <- unit[keep]
  }
  probs <- c(median = 0.5, q90 = 0.9, q95 = 0.95)
  side <- function(direction) {
    sgn <- if (direction == "over") 1 else -1
    in_side <- sgn * rel_error > 0
    x <- sgn * rel_error[in_side]
    out <- if (length(x) == 0L) {
      tibble::tibble(
        direction = direction,
        prop = if (length(rel_error)) 0 else NA_real_,
        median = NA_real_,
        q90 = NA_real_,
        q95 = NA_real_,
        max = NA_real_
      )
    } else {
      tibble::tibble(
        direction = direction,
        prop = length(x) / length(rel_error),
        median = sgn * stats::median(x),
        q90 = sgn * stats::quantile(x, probs = 0.9, names = FALSE),
        q95 = sgn * stats::quantile(x, probs = 0.95, names = FALSE),
        max = sgn * max(x)
      )
    }
    if (!isTRUE(mcse)) {
      return(out)
    }
    for (s in names(probs)) {
      if (is.null(unit)) {
        q <- .analysis_mcse_quantile(x, probs[[s]])
        bounds <- sgn * c(q$lower, q$upper)
        out[[paste0(s, "_mcse")]] <- q$mcse
        out[[paste0(s, "_lower")]] <- min(bounds)
        out[[paste0(s, "_upper")]] <- max(bounds)
      } else {
        p <- probs[[s]]
        se <- .analysis_mcse_between_units(
          x, unit[in_side], function(v) .analysis_mcse_quantile_finite(v, p)
        )
        out[[paste0(s, "_mcse")]] <- se
        out[[paste0(s, "_lower")]] <- out[[s]] - .analysis_mcse_z * se
        out[[paste0(s, "_upper")]] <- out[[s]] + .analysis_mcse_z * se
      }
    }
    out$max_mcse <- NA_real_
    out$max_lower <- NA_real_
    out$max_upper <- NA_real_
    .simBandwidthSignedErrorClipSide(out)
  }
  dplyr::bind_rows(side("over"), side("under"))
}

# Keep Monte Carlo bounds on their direction's side of zero.
.simBandwidthSignedErrorClipSide <- function(tbl) {
  over <- tbl$direction == "over"
  for (s in c("median", "q90", "q95")) {
    lo <- paste0(s, "_lower")
    hi <- paste0(s, "_upper")
    if (lo %in% names(tbl)) {
      tbl[[lo]] <- dplyr::if_else(over, pmax(tbl[[lo]], 0), tbl[[lo]])
    }
    if (hi %in% names(tbl)) {
      tbl[[hi]] <- dplyr::if_else(!over, pmin(tbl[[hi]], 0), tbl[[hi]])
    }
  }
  tbl
}

# One row per group and direction; `tbl` must have a `rel_error` column.
# `mcse` and `unit` (a column name) are passed to
# `.simBandwidthSignedErrorSides()`.
.simBandwidthSignedErrorSummary <- function(
  tbl,
  group_cols,
  mcse = FALSE,
  unit = NULL
) {
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::reframe(.simBandwidthSignedErrorSides(
      .data$rel_error,
      mcse = mcse,
      unit = if (is.null(unit)) NULL else .data[[unit]]
    ))
}

# Average side summaries equally over scenarios, ignoring empty sides. When
# the summaries carry Monte Carlo errors (`<stat>_mcse`), these are combined
# as for independent scenarios (`.analysis_mcse_average()`), with bounds at
# the average +/- 1.96 combined MCSE, kept on their side of zero.
.simBandwidthSignedErrorAverage <- function(tbl, group_cols) {
  stats <- intersect(c("prop", "median", "q90", "q95", "max"), names(tbl))
  if (any(paste0(stats, "_mcse") %in% names(tbl))) {
    out <- .analysis_mcse_average_cols(tbl, c(group_cols, "direction"), stats)
    out <- out[, setdiff(names(out), c("prop_mcse", "prop_lower", "prop_upper"))]
    return(.simBandwidthSignedErrorClipSide(out))
  }
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "direction")))) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::any_of(c("prop", "median", "q90", "q95", "max")),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(
      dplyr::any_of(c("prop", "median", "q90", "q95", "max")),
      ~ dplyr::if_else(is.nan(.x), NA_real_, .x)
    ))
}

# Analysis 2a absolute relative error: median, 95th percentile and maximum of
# |estimate - truth| / truth within each scenario (`scenario_cols`), then
# averaged equally over the scenarios in each `group_cols` group, in percent
# (three significant figures). Monte Carlo bounds (`<stat>_lower`/`_upper`,
# `analysis-mcse.R`) combine each scenario's order-statistic half-width as for
# independent scenarios; the maximum has none.
.simBandwidthGlobalAbsErrorAverage <- function(tbl, scenario_cols, group_cols) {
  stats <- c("err_rel_median_avg", "err_rel_95_avg", "err_rel_max_avg")
  bounds <- as.vector(outer(stats, c("_lower", "_upper"), paste0))
  tbl |>
    dplyr::mutate(
      .abs_rel = abs(.data$propRespEst - .data$propRespTruth) /
        .data$propRespTruth
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenario_cols))) |>
    dplyr::summarise(
      err_rel_median_avg = stats::median(.data$.abs_rel, na.rm = TRUE),
      err_rel_median_avg_mcse = .analysis_mcse_quantile(.data$.abs_rel, 0.5)$mcse,
      err_rel_95_avg = stats::quantile(
        .data$.abs_rel,
        probs = 0.95, na.rm = TRUE, names = FALSE
      ),
      err_rel_95_avg_mcse = .analysis_mcse_quantile(.data$.abs_rel, 0.95)$mcse,
      err_rel_max_avg = max(.data$.abs_rel, na.rm = TRUE),
      err_rel_max_avg_mcse = NA_real_,
      .groups = "drop"
    ) |>
    .analysis_mcse_average_cols(group_cols, stats) |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(stats), ~ signif(.x, digits = 3) * 1e2),
      dplyr::across(dplyr::all_of(bounds), ~ pmax(.x, 0) * 1e2)
    ) |>
    dplyr::select(-dplyr::all_of(paste0(stats, "_mcse")))
}

# Monte Carlo error bars for a long signed (or absolute) error table with
# `lower`/`upper` in data units, squished to the +1500% cap like the values.
.simBandwidthSignedErrorBars <- function(tbl, width = 0) {
  if (!all(c("lower", "upper") %in% names(tbl))) {
    return(NULL)
  }
  tbl$lower_shown <- .simBandwidthSignedErrorSquish(tbl$lower)
  tbl$upper_shown <- .simBandwidthSignedErrorSquish(tbl$upper)
  .analysis_mcse_errorbar(tbl, "lower_shown", "upper_shown", width = width)
}

# Under-estimates are linear down to -100% (nothing gated); over-estimates are
# on a log2 fold scale, so -100% and +100% (two-fold) are equally far from zero.
.simBandwidthSignedErrorTrans <- function() {
  scales::trans_new(
    "signed_rel_error",
    transform = function(x) {
      pos <- !is.na(x) & x > 0
      x[pos] <- log2(1 + x[pos])
      x
    },
    inverse = function(x) {
      pos <- !is.na(x) & x > 0
      x[pos] <- 2^x[pos] - 1
      x
    },
    breaks = function(limits) {
      lo <- max(limits[1], -1, na.rm = TRUE)
      hi <- max(limits[2], 0, na.rm = TRUE)
      # Close to zero the scale is near-linear, so ordinary breaks suffice.
      if (hi <= 1) {
        return(pretty(c(lo, hi)))
      }
      c(pretty(c(lo, 0), n = 3), 2^seq_len(ceiling(log2(1 + hi))) - 1)
    },
    domain = c(-1, Inf)
  )
}

# Over-estimates are capped at +1500% (16 times the truth): larger errors are
# drawn at the cap, whose tick then reads ">= +1500% (16x)". The cap is a
# doubling, so it is always one of the scale's breaks.
.simBandwidthSignedErrorCap <- 15

.simBandwidthSignedErrorSquish <- function(x, cap = .simBandwidthSignedErrorCap) {
  pmin(x, cap)
}

.simBandwidthSignedErrorIsCapped <- function(x, cap = .simBandwidthSignedErrorCap) {
  any(x > cap, na.rm = TRUE)
}

# `cap`: errors drawn at this value may be larger, so its label gets a ">=" sign.
.simBandwidthSignedErrorLabel <- function(x, cap = Inf) {
  lab <- ifelse(x > 0, sprintf("+%g%%", 100 * x), sprintf("%g%%", 100 * x))
  # Nothing gated (-100%, 0x) and each doubling (+100% 2x, +300% 4x, ...)
  # also show the multiple of the true response.
  doublings <- log2(1 + pmax(x, 0))
  fold <- is.finite(x) &
    (abs(x + 1) < 1e-8 | (x > 0 & abs(doublings - round(doublings)) < 1e-8))
  lab[fold] <- paste0(lab[fold], " (", sprintf("%g", 1 + x[fold]), "x)")
  at_cap <- is.finite(x) & is.finite(cap) & x >= cap - 1e-8
  lab[at_cap] <- paste0("\u2265 ", lab[at_cap])
  lab
}

# Absolute relative error on the over-estimate part of the signed scale (fold
# scale above zero), capped at 1500%; the cap's tick then reads ">= 1500%".
.simBandwidthAbsErrorLabel <- function(x, cap = Inf) {
  lab <- .analysis_label_percent(x)
  at_cap <- is.finite(x) & is.finite(cap) & x >= cap - 1e-8
  lab[at_cap] <- paste0("\u2265 ", lab[at_cap])
  lab
}

.simBandwidthAbsErrorLayers <- function(
  y_label = "Absolute relative error",
  capped = FALSE
) {
  cap <- if (isTRUE(capped)) .simBandwidthSignedErrorCap else Inf
  list(
    ggplot2::scale_y_continuous(
      transform = .simBandwidthSignedErrorTrans(),
      labels = function(x) .simBandwidthAbsErrorLabel(x, cap = cap)
    ),
    ggplot2::expand_limits(y = c(0, 1)),
    ggplot2::labs(y = y_label)
  )
}

# ggplot2 cannot vary line width along a dashed line, so draw each line as
# segments between neighbouring points, weighted by the mean share at their ends.
.simBandwidthSignedErrorSegmentLayer <- function(
  tbl,
  x,
  y,
  line_cols,
  alpha = 1
) {
  segments <- tbl |>
    dplyr::group_by(dplyr::across(dplyr::any_of(line_cols))) |>
    dplyr::arrange(.data[[x]], .by_group = TRUE) |>
    dplyr::mutate(
      x_end = dplyr::lead(.data[[x]]),
      y_end = dplyr::lead(.data[[y]]),
      prop_segment = (.data$prop + dplyr::lead(.data$prop)) / 2
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$x_end))
  ggplot2::geom_segment(
    data = segments,
    ggplot2::aes(xend = x_end, yend = y_end, linewidth = prop_segment),
    lineend = "round",
    alpha = alpha
  )
}

# Shared y scale, zero line and line-weight scale for signed-error plots.
# `capped`: some errors were drawn at the +1500% cap, so label that tick with ">=".
.simBandwidthSignedErrorLayers <- function(
  y_label = "Relative error",
  capped = FALSE
) {
  cap <- if (isTRUE(capped)) .simBandwidthSignedErrorCap else Inf
  list(
    ggplot2::geom_hline(yintercept = 0, colour = "grey40"),
    ggplot2::scale_y_continuous(
      transform = .simBandwidthSignedErrorTrans(),
      labels = function(x) .simBandwidthSignedErrorLabel(x, cap = cap)
    ),
    # Always show losing the whole response (-100%) and doubling it (+100%).
    ggplot2::expand_limits(y = c(-1, 1)),
    ggplot2::scale_linewidth_continuous(
      range = c(0.4, 2),
      limits = c(0, 1),
      labels = scales::percent
    ),
    ggplot2::labs(
      y = y_label,
      linewidth = "Share of estimates\nin this direction"
    )
  )
}
