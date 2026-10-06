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
  default_palette <- is.null(base_col_vec) || length(base_col_vec) == 0L
  if (is.null(base_col_vec) || length(base_col_vec) == 0L) {
    base_col_vec <- c(
      "#8c96c6", "#810f7c"
    )
  }

  bw_num <- sort(unique(as.numeric(bw_vec)))
  bw_lab <- format_bw_lab(bw_num)
  n_bw <- length(bw_num)

  if (default_palette) {
    col_vec <- grDevices::colorRampPalette(base_col_vec)(n_bw)
  } else if (n_bw <= length(base_col_vec)) {
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
  tbl <- dplyr::select(tbl, -dplyr::starts_with(".boot_"))
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
  if ("n_scenario_median" %in% names(tbl)) {
    averaged_labels <- c(
      Median = "Mean of scenario medians",
      "90th percentile" = "Mean of scenario 90th percentiles",
      "95th percentile" = "Mean of scenario 95th percentiles",
      Maximum = "Mean of scenario maxima"
    )
    matched <- stat_cols %in% names(averaged_labels)
    stat_cols[matched] <- averaged_labels[stat_cols[matched]]
  }
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

# Rank within each transformation; retain every simulated bandwidth.
.simBandwidthRankData <- function(tbl) {
  groups <- if ("transformation" %in% names(tbl)) as.character(tbl$transformation) else rep("All", nrow(tbl))
  rank <- stats::ave(tbl$bw, groups, FUN = function(x) match(x, sort(unique(x))))
  k <- if (any(is.finite(rank))) max(rank, na.rm = TRUE) else 0L
  tbl$bw_rank <- factor(rank, levels = seq_len(k))
  tbl
}

.simBandwidthRankScales <- function(tbl) {
  if (!"bw_rank" %in% names(tbl)) tbl <- .simBandwidthRankData(tbl)
  all_levels <- levels(tbl$bw_rank)
  levels <- all_levels[all_levels %in% as.character(tbl$bw_rank)]
  labels <- vapply(levels, function(rank) {
    rows <- unique(tbl[as.character(tbl$bw_rank) == rank, intersect(c("transformation", "bw"), names(tbl)), drop = FALSE])
    if (!"transformation" %in% names(rows)) return(paste0(rank, ": ", format_bw_lab(rows$bw[1])))
    values <- split(as.character(.analysis_trans_factor(rows$transformation)), format_bw_lab(rows$bw))
    paste0(rank, ": ", paste(vapply(names(values), function(bw) {
      paste0(bw, " (", paste(values[[bw]], collapse = "/"), ")")
    }, character(1)), collapse = ", "))
  }, character(1))
  guide <- ggplot2::guide_legend(ncol = 1, override.aes = list(alpha = 1, size = 2))
  list(
    ggplot2::scale_colour_manual(values = make_bw_colour_values(seq_along(all_levels))[levels],
      breaks = levels, labels = labels, name = "Bandwidth rank", guide = guide),
    ggplot2::scale_shape_manual(values = stats::setNames(rep(c(16, 15, 18, 3, 7, 8, 0, 1, 2, 5, 6, 9, 10, 12), length.out = length(all_levels)), all_levels)[levels],
      breaks = levels, labels = labels, name = "Bandwidth rank", guide = guide)
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
# although only the bandwidth rule is shown in Analysis 2b.
# `mcse`: draw the `<stat>_lower`/`<stat>_upper` Monte Carlo intervals.
.simBandwidthBiasRelativeErrorPlot <- function(
  tbl,
  y_label = "Absolute relative error",
  facet = NULL,
  stat_cols = c(
    median_abs_rel_error = "Median",
    q90_abs_rel_error = "90th percentile",
    max_abs_rel_error = "Maximum"
  ),
  mcse = FALSE
) {
  if (is.null(facet)) facet <- ggplot2::facet_wrap(ggplot2::vars(statistic, mismatch_label),
    ncol = dplyr::n_distinct(tbl$mismatch_label), scales = "free_y", labeller = ggplot2::label_both)
  tbl <- .simBandwidthRankData(tbl) |>
    .simBandwidthErrorStatLong(stat_cols) |>
    dplyr::mutate(bw_lab = .simBandwidthBwLabFactor(.data$bw))
  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = bias_uns_multiplier,
      y = value,
      colour = bw_rank,
      shape = bw_rank,
      group = interaction(bw, bias_uns_basis)
    )
  ) +
    (if (isTRUE(mcse) && "lower" %in% names(tbl)) .analysis_mcse_errorbar(tbl)) +
    ggplot2::geom_line(alpha = 0.75) +
    ggplot2::geom_point(alpha = 0.75) +
    facet +
    .simBandwidthRankScales(tbl) +
    ggplot2::scale_x_continuous(labels = .analysis_label_number) +
    ggplot2::scale_y_continuous(labels = .analysis_label_percent, guide = ggplot2::guide_axis(check.overlap = TRUE)) +
    .analysis_theme() +
    ggplot2::labs(
      x = "Bias multiplier",
      y = y_label,
      colour = "Bandwidth",
      caption = .simBandwidthScenarioCaption(tbl)
    )
}

# Signed-error version of the bias curves. `tbl` comes from
# `.simBandwidthSignedErrorSummary()`: over-estimates sit above zero and
# under-estimates below; line weight is each direction's share. Errors above
# +1500% are drawn at +1500% (`value_shown`); `value` keeps the actual error.
.simBandwidthBiasSignedErrorPlot <- function(
  tbl,
  y_label = "Relative error",
  facet = NULL,
  stat_cols = c(median = "Median", q90 = "90th percentile", max = "Maximum"),
  mcse = FALSE
) {
  if (is.null(facet)) facet <- ggplot2::facet_wrap(ggplot2::vars(statistic, mismatch_label),
    ncol = dplyr::n_distinct(tbl$mismatch_label), scales = "free_y", labeller = ggplot2::label_both)
  tbl <- .simBandwidthRankData(tbl) |>
    .simBandwidthErrorStatLong(stat_cols) |>
    dplyr::mutate(
      bw_lab = .simBandwidthBwLabFactor(.data$bw),
      value_shown = .simBandwidthSignedErrorSquish(.data$value)
    )
  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = bias_uns_multiplier,
      y = value_shown,
      colour = bw_rank,
      shape = bw_rank,
      group = interaction(bw, bias_uns_basis, direction)
    )
  ) +
    (if (isTRUE(mcse)) .simBandwidthSignedErrorBars(tbl)) +
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
    facet +
    .simBandwidthRankScales(tbl) +
    ggplot2::scale_x_continuous(labels = .analysis_label_number) +
    .analysis_theme() +
    ggplot2::labs(
      x = "Bias multiplier", colour = "Bandwidth",
      caption = .simBandwidthDisplayCaption(tbl$value, .simBandwidthScenarioCaption(tbl))
    )
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
# +1500% are drawn at +1500% (`err_value_shown`).
# `mcse`: draw the Monte Carlo intervals carried by the summary.
.simBandwidthGlobalSignedErrorPlot <- function(
  tbl,
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
    (if (isTRUE(mcse)) .simBandwidthSignedErrorBars(tbl)) +
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
    (if (by_prob) {
      ggplot2::facet_wrap(
        ggplot2::vars(prob_response, transformation),
        ncol = dplyr::n_distinct(tbl$transformation), scales = "free",
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
      labels = if ("n_scenario" %in% names(tbl)) c(
        over_median = "Over: mean median",
        over_q95 = "Over: mean 95th",
        over_max = "Over: mean maximum",
        under_median = "Under: mean median",
        under_q95 = "Under: mean 95th",
        under_max = "Under: mean maximum"
      ) else c(
        over_median = "Over: median",
        over_q95 = "Over: 95th percentile",
        over_max = "Over: maximum",
        under_median = "Under: median",
        under_q95 = "Under: 95th percentile",
        under_max = "Under: maximum"
      ),
      drop = FALSE,
      guide = ggplot2::guide_legend(ncol = 1)
    ) +
    ggplot2::labs(
      x = "Bandwidth", colour = NULL,
      y = if ("n_scenario" %in% names(tbl)) "Mean of scenario statistics (relative error)" else "Relative error",
      caption = .simBandwidthDisplayCaption(tbl$err_value, .simBandwidthScenarioCaption(tbl))
    ) +
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
# percentile is unchanged; uncertainty resamples all datasets and recomputes
# that same pooled percentile. Conditional draws without this direction remain
# undefined. Five contributing datasets and 95% finite draws are required.
.simBandwidthSignedErrorSides <- function(
    rel_error, mcse = FALSE, unit = NULL, bootstrap_family = "default") {
  if (!is.null(unit)) {
    direction_stats <- function(direction) {
      sign <- if (direction == "over") 1 else -1
      stat <- function(p) function(v) {
        x <- sign * v[is.finite(v) & sign * v > 0]
        .analysis_mcse_quantile_finite(x, p) * sign
      }
      share <- function(v) {
        x <- v[is.finite(v)]
        if (length(x)) mean(sign * x > 0) else NA_real_
      }
      specs <- list(prop = share, median = stat(0.5), q90 = stat(0.9), q95 = stat(0.95))
      out <- dplyr::bind_cols(purrr::imap(specs, function(fn, name) {
        .analysis_mcse_pooled_cols(rel_error, unit, fn, name, bootstrap_family, mcse)
      }))
      out$direction <- direction
      out$n_dataset_total <- dplyr::n_distinct(unit, na.rm = TRUE)
      .simBandwidthSignedErrorClipSide(out)
    }
    return(dplyr::bind_rows(direction_stats("over"), direction_stats("under")))
  }
  keep <- is.finite(rel_error)
  rel_error <- rel_error[keep]
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
      q <- .analysis_mcse_quantile(x, probs[[s]])
      bounds <- sgn * c(q$lower, q$upper)
      out[[paste0(s, "_mcse")]] <- q$mcse
      out[[paste0(s, "_lower")]] <- min(bounds)
      out[[paste0(s, "_upper")]] <- max(bounds)
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
  if (".boot_median" %in% names(tbl)) {
    stats <- intersect(c("prop", "median", "q90", "q95"), names(tbl))
    mcse <- any(lengths(tbl$.boot_median) > 0L)
    return(.simBandwidthSignedErrorClipSide(.analysis_mcse_bootstrap_average(
      tbl, c(group_cols, "direction"), stats, mcse)))
  }
  stats <- intersect(c("prop", "median", "q90", "q95", "max"), names(tbl))
  if (any(paste0(stats, "_mcse") %in% names(tbl))) {
    out <- .analysis_mcse_average_cols(tbl, c(group_cols, "direction"), stats)
    out <- out[, setdiff(names(out), c("prop_mcse", "prop_lower", "prop_upper"))]
    return(.simBandwidthSignedErrorClipSide(out))
  }
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "direction")))) |>
    dplyr::summarise(
      n_scenario = dplyr::n(),
      dplyr::across(
        dplyr::any_of(c("prop", "median", "q90", "q95", "max")),
        ~ sum(is.finite(.x)), .names = "n_scenario_{.col}"
      ),
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
# `abs_rel_error` (|estimate - truth| / truth) within each scenario
# (`scenario_cols`), then
# averaged equally over the scenarios in each `group_cols` group, in percent
# (three significant figures). Monte Carlo bounds (`<stat>_lower`/`_upper`,
# `analysis-mcse.R`) combine each scenario's order-statistic half-width as for
# independent scenarios; the maximum has none.
.simBandwidthGlobalAbsErrorAverage <- function(tbl, scenario_cols, group_cols) {
  stats <- c("err_rel_median_avg", "err_rel_95_avg", "err_rel_max_avg")
  bounds <- as.vector(outer(stats, c("_lower", "_upper"), paste0))
  tbl |>
    dplyr::mutate(.abs_rel = .data$abs_rel_error) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenario_cols))) |>
    dplyr::summarise(
      err_rel_median_avg = stats::median(.data$.abs_rel, na.rm = TRUE),
      err_rel_median_avg_mcse = .analysis_mcse_quantile(.data$.abs_rel, 0.5)$mcse,
      err_rel_95_avg = stats::quantile(
        .data$.abs_rel,
        probs = 0.95, na.rm = TRUE, names = FALSE
      ),
      err_rel_95_avg_mcse = .analysis_mcse_quantile(.data$.abs_rel, 0.95)$mcse,
      err_rel_max_avg = if (any(is.finite(.data$.abs_rel))) {
        max(.data$.abs_rel, na.rm = TRUE)
      } else NA_real_,
      err_rel_max_avg_mcse = NA_real_,
      .groups = "drop"
    ) |>
    .analysis_mcse_average_cols(group_cols, stats) |>
    dplyr::rename(
      n_scenario_median = n_scenario_err_rel_median_avg,
      n_scenario_q95 = n_scenario_err_rel_95_avg,
      n_scenario_max = n_scenario_err_rel_max_avg
    ) |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(stats), ~ signif(.x, digits = 3) * 1e2),
      dplyr::across(dplyr::all_of(bounds), ~ pmax(.x, 0) * 1e2)
    ) |>
    dplyr::select(-dplyr::all_of(paste0(stats, "_mcse")))
}

# Monte Carlo error bars for a long signed (or absolute) error table with
# `lower`/`upper` in data units, squished to the +1500% cap like the values.
.simBandwidthSignedErrorBars <- function(tbl, width = 0.08) {
  if (!all(c("lower", "upper") %in% names(tbl))) {
    return(NULL)
  }
  tbl$lower_shown <- .simBandwidthSignedErrorSquish(tbl$lower)
  tbl$upper_shown <- .simBandwidthSignedErrorSquish(tbl$upper)
  .analysis_mcse_errorbar(tbl, "lower_shown", "upper_shown", width = width)
}

# Under-estimates are linear from -100% to zero and logarithmic below -100%;
# over-estimates are on a log2 fold scale, so -100% and +100% (two-fold) are equally far from zero.
.simBandwidthSignedErrorTrans <- function() {
  scales::trans_new(
    "signed_rel_error",
    transform = function(x) {
      neg <- !is.na(x) & x < -1
      pos <- !is.na(x) & x > 0
      x[neg] <- -1 - log2(-x[neg])
      x[pos] <- log2(1 + x[pos])
      x
    },
    inverse = function(x) {
      neg <- !is.na(x) & x < -1
      pos <- !is.na(x) & x > 0
      x[neg] <- -2^(-1 - x[neg])
      x[pos] <- 2^x[pos] - 1
      x
    },
    breaks = function(limits) {
      lo <- min(limits[1], 0, na.rm = TRUE)
      hi <- max(limits[2], 0, na.rm = TRUE)
      below <- if (lo < -1) -2^seq_len(ceiling(log2(-lo))) else numeric()
      lo <- max(lo, -1)
      # Close to zero the scale is near-linear, so ordinary breaks suffice.
      if (hi <= 1) {
        return(sort(unique(c(below, pretty(c(lo, hi)), if (lo <= -1) -1))))
      }
      sort(unique(c(
        below, pretty(c(lo, 0), n = 3), if (lo <= -1) -1,
        2^seq_len(ceiling(log2(1 + hi))) - 1
      )))
    },
    domain = c(-Inf, Inf)
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
  # An estimate of zero (-100%, 0x) and each doubling (+100% 2x, +300% 4x, ...)
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
      breaks = function(limits) {
        b <- .simBandwidthSignedErrorTrans()$breaks(pmax(limits, 0))
        b[b >= 0]
      },
      guide = ggplot2::guide_axis(check.overlap = TRUE),
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

.simBandwidthCapPoints <- function() {
  ggplot2::geom_point(data = function(data) {
    column <- intersect(c("value", "err_value", "error", "estimate"), names(data))
    if (!length(column)) return(data[FALSE, , drop = FALSE])
    data[is.finite(data[[column[1]]]) & data[[column[1]]] > .simBandwidthSignedErrorCap, , drop = FALSE]
  }, shape = 17, size = 2, show.legend = FALSE)
}

# Printed beside the figure and logged with its full output path.
.simBandwidthDisplayNote <- function(plot, path) {
  scale <- plot$scales$get_scales("y")
  if (is.null(scale) || !identical(scale$trans$name, "signed_rel_error")) return(invisible(NULL))
  data <- .simBandwidthPlotData(plot)
  columns <- intersect(c("value", "err_value"), names(data))
  points <- unlist(data[columns], use.names = FALSE)
  intervals <- lapply(plot$layers, function(layer) {
    if (!inherits(layer$geom, "GeomErrorbar") || !is.data.frame(layer$data)) return(NULL)
    cols <- intersect(c("lower", "upper", "lower_shown", "upper_shown"), names(layer$data))
    unlist(layer$data[cols], use.names = FALSE)
  })
  values <- c(points, unlist(intervals, use.names = FALSE))
  notes <- c(if (any(points > .simBandwidthSignedErrorCap, na.rm = TRUE)) "triangle: above display cap",
    if (any(values < -1, na.rm = TRUE)) "values below -100% (negative estimates) are drawn on a compressed log scale")
  if (length(notes)) {
    note <- paste(notes, collapse = "; ")
    cat("\n", note, ".\n\n", sep = "")
    message(path, ": ", note)
  }
  invisible(NULL)
}

.simBandwidthDisplayCaption <- function(values, caption = NULL) {
  paste(c(caption,
    if (any(values > .simBandwidthSignedErrorCap, na.rm = TRUE)) "triangle: above display cap",
    if (any(values < -1, na.rm = TRUE)) "values below -100% (negative estimates) are drawn on a compressed log scale"),
    collapse = "\n")
}

# Shared y scale, zero line and line-weight scale for signed-error plots.
# `capped`: some errors were drawn at the +1500% cap, so label that tick with ">=".
.simBandwidthSignedErrorLayers <- function(
  y_label = "Relative error",
  capped = FALSE
) {
  cap <- if (isTRUE(capped)) .simBandwidthSignedErrorCap else Inf
  list(
    .simBandwidthCapPoints(),
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
      caption = if (isTRUE(capped)) "triangle: above display cap",
      linewidth = "Share of estimates\nin this direction"
    )
  )
}

# Keep signed-error coordinates and intervals; only express ticks as estimate/truth.
.simBandwidthRatioPlot <- function(plot) {
  scale <- plot$scales$get_scales("y")$clone()
  capped <- startsWith(scale$labels(.simBandwidthSignedErrorCap), "\u2265")
  scale$labels <- function(x) {
    lab <- paste0(.analysis_label_number(1 + x), "x")
    lab[is.finite(x) & abs(x) < 1e-8] <- "1x (exact)"
    at_cap <- isTRUE(capped) & is.finite(x) & x >= .simBandwidthSignedErrorCap - 1e-8
    lab[at_cap] <- paste0("\u2265 ", lab[at_cap])
    lab
  }
  plot$scales <- plot$scales$clone()
  plot$scales$scales[[which(vapply(plot$scales$scales,
    function(x) "y" %in% x$aesthetics, logical(1)))]] <- scale
  plot + ggplot2::labs(y = "Estimate / reference (multiple)")
}

# Signed-error companions are saved only, in sibling ratio folders with the same dimensions.
.simBandwidthPrintRatioTwin <- function(plot, path, height, level = 6L, allow_tall = FALSE, mcse_mode = NULL) {
  ratio <- .simBandwidthRatioPlot(plot)
  .analysis_save_fig(ratio, sub("signed_error", "ratio", path, fixed = TRUE),
    height = height, allow_tall = allow_tall, mcse_mode = mcse_mode)
  invisible(NULL)
}

# Generating distributions for threshold companions, computed once per distinct
# biological setting in a render. Cell counts, bandwidth and bias are not part
# of this reference: these are not densities of the gated simulation samples.
# By default the unstimulated reference is its negative component (bandwidth
# analyses); use FALSE for whole-tube comparison references such as QMD 7.
.simBandwidthThresholdDensities <- function(
    panels, settings, n_cell = 1e5, seed = 271L, density_n = 2048L,
    unstimulated_negative_only = TRUE) {
  keys <- c(
    "transformation", "mean_pos", "prob_response",
    "sample_perturbation_sd", "condition_perturbation_sd",
    "cluster_perturbation_sd", "background_relative_to_response"
  )
  panels <- panels |>
    dplyr::select(dplyr::all_of(keys)) |>
    dplyr::distinct()
  purrr::map_dfr(seq_len(nrow(panels)), function(i) {
    row <- panels[i, ]
    out <- .analysis_with_seed(seed, {
      prob_uns <- row$prob_response * row$background_relative_to_response
      simcyto::simCytExperiment(
        nSample = 1L, nMarker = 1L, nCondition = 2L, nCluster = 2L,
        nCellByCondition = c(n_cell, n_cell),
        transformationFunc = .simMiscGetTrans(as.character(row$transformation)),
        mixtureType = "gaussianOnly",
        meanExprMat = matrix(c(0, row$mean_pos), ncol = 1),
        clusterLabelVec = c("gn", "gp"),
        probVecUns = c(1 - prob_uns, prob_uns),
        probResponseVecByStimCondition = list(c(-row$prob_response, row$prob_response)),
        probExact = settings$probExact,
        covEvMin = settings$covEvMin, covEvMax = settings$covEvMax,
        samplePerturbationSd = row$sample_perturbation_sd,
        conditionPerturbationSd = row$condition_perturbation_sd,
        clusterPerturbationSd = row$cluster_perturbation_sd
      )
    })
    x_uns <- flowCore::exprs(out$flowFrameList[[1]])[, 1]
    if (isTRUE(unstimulated_negative_only)) {
      x_uns <- x_uns[out$labelsList[[1]] == "gn"]
    }
    x_stim <- flowCore::exprs(out$flowFrameList[[2]])[, 1]
    bw <- if (row$transformation == "gamma") 0.025 else 0.15
    limits <- range(c(x_uns, x_stim), finite = TRUE) + c(-3, 3) * bw
    purrr::map_dfr(c("unstimulated", "stimulated"), function(condition) {
      x <- if (condition == "unstimulated") x_uns else x_stim
      dens <- stats::density(
        x, bw = bw, n = density_n, from = limits[1], to = limits[2]
      )
      dplyr::bind_cols(
        row[rep(1L, length(dens$x)), ],
        tibble::tibble(condition = condition, expression = dens$x, density = dens$y)
      )
    })
  }) |>
    dplyr::mutate(transformation = .analysis_trans_factor(.data$transformation))
}

# Context is a separate figure: no threshold marks are drawn over densities.
.simBandwidthThresholdDensityPlot <- function(plot, panels, densities) {
  keys <- setdiff(intersect(names(panels), names(densities)), c("condition", "expression", "density"))
  density_panels <- dplyr::semi_join(densities, dplyr::distinct(panels[, keys, drop = FALSE]), by = keys)
  labels <- c(unstimulated = "Unstimulated negative component", stimulated = "Stimulated (all cells)")
  ggplot2::ggplot(density_panels, ggplot2::aes(x = expression, y = density, colour = condition, fill = condition)) +
    ggplot2::geom_area(alpha = 0.18, position = "identity") +
    ggplot2::geom_line(linewidth = 0.3, alpha = 0.75) +
    ggplot2::facet_wrap(ggplot2::vars(prob_response, mean_pos), scales = "free",
      labeller = ggplot2::labeller(prob_response = .analysis_labeller_percent("Response: "), mean_pos = function(x) paste0("Mean: ", .analysis_label_number(as.numeric(x))))) +
    ggplot2::scale_y_sqrt(labels = .analysis_label_number) +
    ggplot2::scale_colour_manual(values = .simMiscGetStimColVec(), labels = labels) +
    ggplot2::scale_fill_manual(values = .simMiscGetStimColVec(), labels = labels) +
    .analysis_theme() +
    ggplot2::labs(x = "Marker expression", y = "Density, square-root scale", colour = NULL, fill = NULL)
}

# Ordered bandwidth rows, horizontal IQR and median; optional raw thresholds
# jitter only vertically, preserving the threshold value on the x axis.
.simBandwidthThresholdPlot <- function(tbl, samples = FALSE) {
  keys <- intersect(c("sim_id", "transformation", "prob_response", "n_cell", "bw"), names(tbl))
  # The median view matches promoted valid-estimate summaries; the sample
  # view retains the original all-threshold median/IQR convention.
  summary_source <- tbl
  if (!isTRUE(samples) && "valid_estimate" %in% names(tbl)) {
    summary_source$threshold[!tbl$valid_estimate | is.na(tbl$valid_estimate)] <- NA_real_
  }
  summary <- summary_source |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(threshold_median = stats::median(threshold, na.rm = TRUE),
      threshold_iqr_lower = stats::quantile(threshold, 0.25, na.rm = TRUE, names = FALSE),
      threshold_iqr_upper = stats::quantile(threshold, 0.75, na.rm = TRUE, names = FALSE), .groups = "drop")
  prepare <- function(data) dplyr::mutate(data,
    transformation = .analysis_trans_factor(transformation),
    bw_lab = .simBandwidthBwLabFactor(bw))
  summary <- prepare(summary)
  tbl <- prepare(tbl)
  ggplot2::ggplot(summary, ggplot2::aes(x = threshold_median, y = bw_lab)) +
    (if (isTRUE(samples)) ggplot2::geom_point(data = tbl, ggplot2::aes(x = threshold),
      position = ggplot2::position_jitter(width = 0, height = 0.12, seed = 271L), size = 0.5, alpha = 0.25)) +
    ggplot2::geom_segment(ggplot2::aes(x = threshold_iqr_lower, xend = threshold_iqr_upper, yend = bw_lab), linewidth = 0.5) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(ggplot2::vars(n_cell, transformation), ncol = dplyr::n_distinct(tbl$transformation),
      scales = "free", labeller = ggplot2::labeller(n_cell = function(x) paste0("Cells: ", .analysis_label_number(as.numeric(x))))) +
    .analysis_theme() + ggplot2::labs(x = "Expression threshold", y = "Bandwidth")
}

# A plot's data, or its layers' data when it was built as ggplot() + layers.
.simBandwidthPlotData <- function(plot) {
  if (is.data.frame(plot$data)) return(plot$data)
  dplyr::bind_rows(lapply(plot$layers, function(layer) {
    if (is.data.frame(layer$data)) layer$data else NULL
  }))
}

# Match a figure's displayed dimensions, then show coverage beside that figure.
# Omitted scenario dimensions are pooled only for this sample-count diagnostic;
# error curves themselves average per-scenario statistics equally.
.simBandwidthCoverageForPlot <- function(plot, summary) {
  dimensions <- c(
    "mean_pos_setting", "bias_uns_setting", "transformation", "prob_response",
    "n_cell", "bw", "bias_uns_basis", "bias_uns_multiplier",
    "mismatch_label", "mismatch_type", "mismatch_val", "bw_mtd", "bw_ncell_upper"
  )
  plot_data <- .simBandwidthPlotData(plot)
  keys <- intersect(dimensions, intersect(names(plot_data), names(summary)))
  if ("transformation" %in% keys) {
    summary$transformation <- .analysis_trans_factor(summary$transformation)
  }
  selected <- plot_data |>
    dplyr::select(dplyr::all_of(keys)) |>
    dplyr::distinct()
  # Character conversion also handles a figure's numeric/factor cell-count axis.
  summary <- summary |>
    dplyr::mutate(dplyr::across(dplyr::all_of(keys), as.character))
  selected <- selected |>
    dplyr::mutate(dplyr::across(dplyr::all_of(keys), as.character))
  # Keep all expected curve settings, even when ErrorStatLong removed a failed
  # setting's NA statistics. Match only the dimensions selecting the figure.
  figure_keys <- intersect(keys, c(
    "mean_pos_setting", "bias_uns_setting", "n_cell", "bw_ncell_upper"
  ))
  if (length(figure_keys) && nrow(selected)) {
    summary <- dplyr::semi_join(
      summary, dplyr::distinct(selected, dplyr::across(dplyr::all_of(figure_keys))),
      by = figure_keys
    )
  }
  summary |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(
      n_scenario = dplyr::n(),
      n_scenario_valid = sum(.data$n_valid > 0L),
      n_sample = sum(.data$n_sample),
      n_valid = sum(.data$n_valid),
      n_failed = sum(.data$n_failed),
      n_provenance = sum(.data$n_provenance),
      n_fallback = sum(.data$n_fallback),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      failure_fraction = .data$n_failed / .data$n_sample,
      fallback_fraction = dplyr::if_else(
        .data$n_provenance > 0L, .data$n_fallback / .data$n_provenance, NA_real_
      )
    )
}

# Compact coverage beside a figure: overall denominators and failure/fallback
# totals, then only the plotted settings with failures or fallbacks (full
# per-setting tables made the reports tens of megabytes).
.simBandwidthPrintCoverage <- function(plot, summary, max_rows = 20L) {
  tbl <- .simBandwidthCoverageForPlot(plot, summary)
  pct <- function(n, d) if (d > 0) paste0(" (", .analysis_label_percent(n / d), ")") else ""
  n_sample <- sum(tbl$n_sample)
  n_provenance <- sum(tbl$n_provenance)
  per_point <- range(tbl$n_sample)
  cat("\n\n**Coverage.** ", nrow(tbl), " plotted settings; ",
    if (per_point[1] == per_point[2]) per_point[1] else paste0(per_point[1], "\u2013", per_point[2]),
    " samples each (", n_sample, " in total). Failed estimates: ", sum(tbl$n_failed),
    pct(sum(tbl$n_failed), n_sample), ". Threshold fallbacks: ", sum(tbl$n_fallback),
    " of ", n_provenance, pct(sum(tbl$n_fallback), n_provenance), ".\n\n", sep = "")
  flagged <- tbl[tbl$n_failed > 0L | tbl$n_fallback > 0L, , drop = FALSE]
  if (nrow(flagged)) {
    # Drop columns that are constant across the flagged rows.
    keep <- vapply(flagged, function(x) length(unique(x)) > 1L, logical(1)) |
      names(flagged) %in% c("n_sample", "n_failed", "n_fallback", "n_provenance")
    shown <- flagged[order(-flagged$n_failed, -flagged$n_fallback), keep, drop = FALSE]
    cat("Settings with failures or fallbacks",
      if (nrow(shown) > max_rows) paste0(" (first ", max_rows, " of ", nrow(shown), ")"),
      ":\n\n", sep = "")
    cat(knitr::kable(utils::head(shown, max_rows), digits = 3), sep = "\n")
    cat("\n\n")
  }
  plot_data <- .simBandwidthPlotData(plot)
  count_cols <- names(plot_data)[startsWith(names(plot_data), "n_scenario_")]
  if (length(count_cols)) {
    ranges <- vapply(count_cols, function(col) {
      x <- plot_data[[col]]
      paste0(sub("^n_scenario_", "", col), " ", min(x, na.rm = TRUE), "\u2013", max(x, na.rm = TRUE))
    }, character(1))
    cat("Contributing scenarios per plotted statistic: ", paste(ranges, collapse = "; "),
      ".\n\n", sep = "")
  }
  invisible(tbl)
}

# "Contributing scenarios" text for a caption: one count when every statistic
# shares it, otherwise one range per statistic.
.simBandwidthScenarioCountText <- function(tbl, cols, prefix) {
  ranges <- vapply(cols, function(col) {
    x <- tbl[[paste0("n_scenario_", col)]]
    if (min(x) == max(x)) as.character(min(x)) else paste0(min(x), "\u2013", max(x))
  }, character(1))
  if (length(unique(ranges)) == 1L) return(paste0(prefix, " per point: ", ranges[[1]], "."))
  paste0(prefix, " per point: ", paste(paste(cols, ranges), collapse = "; "), ".")
}

.simBandwidthScenarioCaption <- function(tbl) {
  if (!nrow(tbl)) return("No finite scenario statistics")
  cols <- intersect(c("median", "q90", "q95", "max"),
                    sub("^n_scenario_", "", names(tbl)[startsWith(names(tbl), "n_scenario_")]))
  if (!length(cols)) return(NULL)
  .simBandwidthScenarioCountText(tbl, cols, "Contributing scenarios")
}

# Unconditional percentiles are computed on raw signed errors, including zeros.
# Both outer pairs are retained until the complete figure chooses its pair.
.simBandwidthSignedErrorProbs <- c(
  q025 = 0.025, q05 = 0.05, q10 = 0.10, q25 = 0.25, median = 0.50,
  q75 = 0.75, q90 = 0.90, q95 = 0.95, q975 = 0.975
)

.simBandwidthSignedErrorPercentiles <- function(
    rel_error, probs = .simBandwidthSignedErrorProbs, mcse = FALSE,
    unit = NULL, bootstrap_family = "default") {
  if (is.null(names(probs)) || any(!is.finite(probs) | probs <= 0 | probs >= 1)) {
    stop("Percentile probabilities must be named and strictly between zero and one.")
  }
  # Eligibility of the outer pair must not change when interval display is off.
  # Compute intervals once; the plot's mcse argument controls their visibility.
  out <- dplyr::bind_cols(purrr::imap(probs, function(p, name) {
    if (!is.null(unit)) {
      return(.analysis_mcse_pooled_cols(rel_error, unit,
        function(v) .analysis_mcse_quantile_finite(v, p), name,
        bootstrap_family, mcse = TRUE))
    }
    point <- tibble::tibble(value = .analysis_mcse_quantile_finite(rel_error, p))
    names(point) <- name
    dplyr::bind_cols(point, .analysis_mcse_quantile_cols(rel_error, p, name))
  }))
  out$n_finite <- sum(is.finite(rel_error))
  out
}

.simBandwidthSignedPercentileSummary <- function(tbl, group_cols, mcse = FALSE) {
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::reframe(.simBandwidthSignedErrorPercentiles(.data$rel_error, mcse = mcse))
}

.simBandwidthSignedPercentileAverage <- function(tbl, group_cols) {
  .analysis_mcse_average_cols(tbl, group_cols, names(.simBandwidthSignedErrorProbs))
}

# Choose one pair for all panels and methods, checking only finite plotted points.
.simBandwidthSignedPercentileOuter <- function(tbl) {
  available <- vapply(c("q025", "q975"), function(name) {
    if (!all(paste0(name, c("_lower", "_upper")) %in% names(tbl))) return(FALSE)
    keep <- is.finite(tbl[[name]])
    all(is.finite(tbl[[paste0(name, "_lower")]][keep]) &
      is.finite(tbl[[paste0(name, "_upper")]][keep]))
  }, logical(1))
  if (all(available)) c("q025", "q975") else c("q05", "q95")
}

# All seven percentiles share a panel. Symmetric pairs share alpha/line width.
# Facets retain scenario dimensions; no direction/share conditioning is applied.
.simBandwidthSignedPercentilePlot <- function(
    tbl, x = "bw", x_label = "Bandwidth", x_log = FALSE,
    by_prob = FALSE, facet = NULL, mcse = FALSE,
    y_label = NULL,
    alphas = c(outer = 0.35, tail = 0.55, quartile = 0.75, median = 1),
    linewidths = c(outer = 0.4, tail = 0.6, quartile = 0.9, median = 1.2)) {
  outer <- .simBandwidthSignedPercentileOuter(tbl)
  fallback <- identical(outer, c("q05", "q95"))
  stats <- c(outer[1], "q10", "q25", "median", "q75", "q90", outer[2])
  labels <- c(if (fallback) "5th" else "2.5th", "10th", "25th", "50th (median)",
    "75th", "90th", if (fallback) "95th" else "97.5th")
  long <- .simBandwidthStatLongBounds(tbl, stats, "percentile", "value")
  long$percentile <- factor(long$percentile, levels = stats, labels = labels)
  multi <- "method" %in% names(long)
  if ("transformation" %in% names(long)) long$transformation <- .analysis_trans_factor(long$transformation)
  long$value_shown <- .simBandwidthSignedErrorSquish(long$value)
  # Bandwidth is an ordered discrete grid, as in the sibling 2a view.
  if (x == "bw") long$bw <- factor(.analysis_label_number(long$bw),
    levels = .analysis_label_number(sort(unique(tbl$bw))))
  line_cols <- intersect(c("method", "percentile", "transformation", "prob_response",
    "n_cell", "mismatch_label", "bw", "bias_uns_basis", "pop"), names(long))
  line_cols <- setdiff(line_cols, x)
  long$series <- interaction(long[line_cols], drop = TRUE)
  if (is.null(facet)) {
    facet <- if (x == "bias_uns_multiplier") {
      ggplot2::facet_wrap(ggplot2::vars(bw, bias_uns_basis, mismatch_label),
        ncol = dplyr::n_distinct(long$mismatch_label), scales = "free_y", labeller = ggplot2::label_both)
    } else if (by_prob) {
      ggplot2::facet_wrap(ggplot2::vars(prob_response, transformation),
        ncol = dplyr::n_distinct(long$transformation), scales = "free",
        labeller = ggplot2::labeller(prob_response = .analysis_labeller_percent("Response: ")))
    } else ggplot2::facet_wrap(~transformation, scales = "free")
  }
  scenario_counts <- grep("^n_scenario_", names(tbl), value = TRUE)
  counts <- if (length(scenario_counts)) .simBandwidthScenarioCountText(tbl,
    sub("^n_scenario_", "", scenario_counts), "Finite contributing scenarios") else NULL
  caption <- paste(counts,
    if (fallback) "5th/95th (too few samples for 2.5th/97.5th intervals)." else "Outer pair: 2.5th/97.5th.",
    "Unavailable intervals are omitted; points remain.")
  caption <- .simBandwidthDisplayCaption(long$value, caption)
  caption <- paste(strwrap(caption, width = 110), collapse = "\n")
  if (is.null(y_label)) y_label <- if ("n_scenario" %in% names(tbl)) {
    "Mean of scenario percentiles (signed relative error)"
  } else "Signed relative error"
  # Alpha, line width and line type all map to the percentile under one legend
  # title, so ggplot2 merges them. Single-method figures also colour by
  # percentile: teal above the median and brown below, matching over/under views.
  tiers <- c("outer", "tail", "quartile", "median", "quartile", "tail", "outer")
  linetypes <- stats::setNames(c(rep("dotted", 3), "solid", rep("dashed", 3)), labels)
  bars <- if (isTRUE(mcse)) .simBandwidthSignedErrorBars(long) else NULL
  # Intervals are solid regardless of their percentile's line type.
  if (!is.null(bars)) bars$aes_params$linetype <- "solid"
  p <- ggplot2::ggplot(long, ggplot2::aes(
    x = .data[[x]], y = .data$value_shown, group = .data$series,
    colour = .data[[if (multi) "method" else "percentile"]],
    alpha = .data$percentile, linewidth = .data$percentile,
    linetype = .data$percentile)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
    bars + ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    .simBandwidthCapPoints() + facet +
    ggplot2::scale_y_continuous(transform = .simBandwidthSignedErrorTrans(),
      labels = function(v) .simBandwidthSignedErrorLabel(v,
        cap = if (.simBandwidthSignedErrorIsCapped(long$value)) .simBandwidthSignedErrorCap else Inf)) +
    ggplot2::expand_limits(y = c(-1, 1)) +
    ggplot2::scale_alpha_manual(values = stats::setNames(unname(alphas[tiers]), labels),
      name = "Percentile", drop = FALSE) +
    ggplot2::scale_linewidth_manual(values = stats::setNames(unname(linewidths[tiers]), labels),
      name = "Percentile", drop = FALSE) +
    ggplot2::scale_linetype_manual(values = linetypes, name = "Percentile", drop = FALSE) +
    (if (multi) .analysis_scale_method() else ggplot2::scale_colour_manual(
      values = stats::setNames(c(rep(.simBandwidthSignedErrorColours[["under_q95"]], 3), "#333333",
        rep(.simBandwidthSignedErrorColours[["over_q95"]], 3)), labels),
      name = "Percentile", drop = FALSE)) +
    (if (x_log) ggplot2::scale_x_log10(breaks = sort(unique(tbl[[x]])), labels = .analysis_label_number)
      else if (is.numeric(long[[x]])) ggplot2::scale_x_continuous(labels = .analysis_label_number)) +
    .analysis_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::labs(x = x_label, y = y_label, caption = caption)
  p
}
