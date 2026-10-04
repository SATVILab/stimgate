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

# Error statistics as rows of panels: `stat_cols` maps wide columns to labels.
.simBandwidthErrorStatLong <- function(tbl, stat_cols) {
  tbl |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(names(stat_cols)),
      names_to = "statistic",
      values_to = "value"
    ) |>
    dplyr::filter(!is.na(.data$value)) |>
    dplyr::mutate(
      statistic = factor(
        unname(stat_cols[.data$statistic]),
        levels = unname(stat_cols)
      )
    )
}

# Bias curves share the same dimensions in all absolute-error views.
.simBandwidthBiasRelativeErrorPlot <- function(
  tbl,
  title,
  y_label = "Absolute relative error",
  facet = ggplot2::facet_grid(statistic ~ mismatch_label, scales = "free_y"),
  stat_cols = c(
    median_abs_rel_error = "Median",
    q90_abs_rel_error = "90th percentile",
    max_abs_rel_error = "Maximum"
  )
) {
  ggplot2::ggplot(
    .simBandwidthErrorStatLong(tbl, stat_cols),
    ggplot2::aes(
      x = bias_uns_multiplier,
      y = value,
      colour = factor(bw),
      linetype = bias_uns_basis,
      group = interaction(bw, bias_uns_basis)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    facet +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "xy") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title = title,
      x = "Bias multiplier",
      y = y_label,
      colour = "Bandwidth",
      linetype = "Bias scale"
    )
}

# Signed-error version of the bias curves. `tbl` comes from
# `.simBandwidthSignedErrorSummary()`: over-estimates sit above zero and
# under-estimates below; line weight is each direction's share.
.simBandwidthBiasSignedErrorPlot <- function(
  tbl,
  title,
  y_label = "Relative error",
  facet = ggplot2::facet_grid(statistic ~ mismatch_label, scales = "free_y"),
  stat_cols = c(median = "Median", q90 = "90th percentile", max = "Maximum")
) {
  tbl <- .simBandwidthErrorStatLong(tbl, stat_cols)
  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = bias_uns_multiplier,
      y = value,
      colour = factor(bw),
      linetype = bias_uns_basis,
      group = interaction(bw, bias_uns_basis, direction)
    )
  ) +
    .simBandwidthSignedErrorLayers(y_label) +
    .simBandwidthSignedErrorSegmentLayer(
      tbl,
      "bias_uns_multiplier",
      "value",
      c(
        "statistic",
        "mismatch_label",
        "n_cell",
        "bw",
        "bias_uns_basis",
        "direction"
      )
    ) +
    ggplot2::geom_point(size = 1) +
    facet +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "xy") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title = title,
      x = "Bias multiplier",
      colour = "Bandwidth",
      linetype = "Bias scale"
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
# With `by_prob`, rows of panels separate response probabilities.
.simBandwidthGlobalSignedErrorPlot <- function(
  tbl,
  title = NULL,
  by_prob = FALSE
) {
  tbl <- tbl |>
    tidyr::pivot_longer(
      cols = c("median", "q95", "max"),
      names_to = "err_type",
      values_to = "err_value"
    ) |>
    dplyr::filter(!is.na(.data$err_value)) |>
    dplyr::mutate(
      err_type = factor(err_type, levels = c("median", "q95", "max")),
      direction = factor(direction, levels = c("over", "under")),
      transformation = factor(
        transformation,
        levels = c("gaussian", "skew", "gamma")
      ),
      bw_fct = factor(bw),
      series = factor(
        paste0(direction, "_", err_type),
        levels = names(.simBandwidthSignedErrorColours)
      )
    )
  transformation_lab <- c(gamma = "Gamma", gaussian = "Gaussian", skew = "Skew")
  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = bw_fct,
      y = err_value,
      colour = series,
      group = interaction(err_type, direction)
    )
  ) +
    .simBandwidthSignedErrorLayers() +
    .simBandwidthSignedErrorSegmentLayer(
      tbl,
      "bw_fct",
      "err_value",
      c("transformation", "prob_response", "err_type", "direction"),
      alpha = 0.75
    ) +
    # Slight transparency shows overlapping lines; legend keys match.
    ggplot2::geom_point(size = 1, alpha = 0.75) +
    (if (by_prob) {
      ggplot2::facet_grid(
        prob_response ~ transformation,
        scales = "free",
        labeller = ggplot2::labeller(
          transformation = transformation_lab,
          prob_response = function(x) paste0("p = ", x)
        )
      )
    } else {
      ggplot2::facet_wrap(
        ~transformation,
        scales = "free",
        labeller = ggplot2::labeller(transformation = transformation_lab)
      )
    }) +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "xy", minor = "y") +
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
    ggplot2::labs(title = title, x = "Bandwidth", colour = NULL) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(
        fill = "white",
        colour = "black"
      ),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 10),
      legend.position = "bottom",
      legend.box = "vertical"
    )
}

# Signed relative error, (estimate - truth) / truth, summarised separately for
# over- and under-estimates. `prop` is each side's share of finite errors; the
# other statistics are signed and use that side's errors only.
.simBandwidthSignedErrorSides <- function(rel_error) {
  rel_error <- rel_error[is.finite(rel_error)]
  side <- function(direction) {
    sgn <- if (direction == "over") 1 else -1
    x <- sgn * rel_error[sgn * rel_error > 0]
    if (length(x) == 0L) {
      return(tibble::tibble(
        direction = direction,
        prop = if (length(rel_error)) 0 else NA_real_,
        median = NA_real_,
        q90 = NA_real_,
        q95 = NA_real_,
        max = NA_real_
      ))
    }
    tibble::tibble(
      direction = direction,
      prop = length(x) / length(rel_error),
      median = sgn * stats::median(x),
      q90 = sgn * stats::quantile(x, probs = 0.9, names = FALSE),
      q95 = sgn * stats::quantile(x, probs = 0.95, names = FALSE),
      max = sgn * max(x)
    )
  }
  dplyr::bind_rows(side("over"), side("under"))
}

# One row per group and direction; `tbl` must have a `rel_error` column.
.simBandwidthSignedErrorSummary <- function(tbl, group_cols) {
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::reframe(.simBandwidthSignedErrorSides(.data$rel_error))
}

# Average side summaries equally over scenarios, ignoring empty sides.
.simBandwidthSignedErrorAverage <- function(tbl, group_cols) {
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

.simBandwidthSignedErrorLabel <- function(x) {
  lab <- ifelse(x > 0, sprintf("+%g%%", 100 * x), sprintf("%g%%", 100 * x))
  # Nothing gated (-100%, 0x) and each doubling (+100% 2x, +300% 4x, ...)
  # also show the multiple of the true response.
  doublings <- log2(1 + pmax(x, 0))
  fold <- is.finite(x) &
    (abs(x + 1) < 1e-8 | (x > 0 & abs(doublings - round(doublings)) < 1e-8))
  lab[fold] <- paste0(lab[fold], " (", sprintf("%g", 1 + x[fold]), "x)")
  lab
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
.simBandwidthSignedErrorLayers <- function(
  y_label = "Relative error"
) {
  list(
    ggplot2::geom_hline(yintercept = 0, colour = "grey40"),
    ggplot2::scale_y_continuous(
      transform = .simBandwidthSignedErrorTrans(),
      labels = .simBandwidthSignedErrorLabel
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
