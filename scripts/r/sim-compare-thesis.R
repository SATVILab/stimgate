# Thesis figures for the simulated method comparison (Analysis 7). One figure
# of each kind per mean-position (separation) setting, drawn at the thesis
# text width (.analysis_save_thesis_fig()) and saved to the projr
# `output-thesis` folder's results subfolder. Needs analysis-plot-style.R and
# sim-bandwidth-analysis-plot.R (signed display scale, threshold references).

.simCompareThesisMethods <- c("stimgate", "fbeta", "tailgate")

# Method labels with F-beta as F_beta (plot maths).
.simCompareThesisMethodLabels <- function() {
  c(stimgate = "StimGate", fbeta = "F[beta]", tailgate = "Tailgate")
}

.simCompareThesisPercent <- function(x) {
  out <- paste0(formatC(100 * x, format = "fg", digits = 2, flag = "#"), "%")
  out <- sub("\\.%$", "%", sub("(\\.[0-9]*?)0+%$", "\\1%", out))
  ifelse(is.finite(x), out, NA_character_)
}

# Response probability strip labels, e.g. "Response 0.2%".
.simCompareThesisProbLab <- function(x) paste("Response", .simCompareThesisPercent(as.numeric(x)))

# Estimated response frequency: median (solid) and 10th and 90th percentiles
# (dotted) of each method's per-sample estimates against the number of cells,
# rows for response probability and columns for transformation. The y axis is
# the signed display scale of the relative error (linear from an estimate of
# zero up to the truth, doubling above it), labelled in each row's actual
# frequencies. Rows are built separately so their tick labels can differ, then
# stacked. The truth is the thick grey line; every row shows zero and twice
# the truth, so the truth sits at most halfway up its row.
.simCompareThesisFreqPlot <- function(pct, pos) {
  trans <- .simBandwidthSignedErrorTrans()
  tbl <- pct[pct$mean_pos_setting == pos & pct$method %in% .simCompareThesisMethods, , drop = FALSE]
  long <- tidyr::pivot_longer(
    tbl[c("method", "transformation", "prob_response", "n_cell", "q10", "median", "q90")],
    c("q10", "median", "q90"), names_to = "stat", values_to = "value"
  )
  long$stat <- factor(long$stat, c("q10", "median", "q90"),
    c("10th percentile", "Median", "90th percentile"))
  long$tier <- factor(ifelse(long$stat == "Median", "Median", "10th and 90th percentiles"),
    c("Median", "10th and 90th percentiles"))
  long$y <- trans$transform(.simBandwidthSignedErrorSquish(long$value))
  long$transformation <- .analysis_trans_factor(long$transformation)
  long$method <- factor(long$method, .simCompareThesisMethods)
  long$series <- interaction(long$method, long$stat, drop = TRUE)
  probs <- sort(unique(long$prob_response), decreasing = TRUE)
  cols <- .analysis_method_colours[.simCompareThesisMethods]
  shapes <- .analysis_method_shapes[.simCompareThesisMethods]
  labs <- .simCompareThesisMethodLabels()
  method_scales <- list(
    ggplot2::scale_colour_manual(values = cols, labels = scales::label_parse()(labs), name = "Method"),
    ggplot2::scale_shape_manual(values = shapes, labels = scales::label_parse()(labs), name = "Method"),
    ggplot2::scale_linetype_manual(values = c(Median = "solid", `10th and 90th percentiles` = "11"),
      name = NULL)
  )
  row_plot <- function(prob, first, last) {
    d <- long[long$prob_response == prob, , drop = FALSE]
    rng <- range(c(d$y, -1, 1), na.rm = TRUE)
    breaks_err <- trans$breaks(trans$inverse(rng))
    breaks_err <- breaks_err[breaks_err >= trans$inverse(rng[1]) - 1e-9 &
      breaks_err <= trans$inverse(rng[2]) + 1e-9]
    ggplot2::ggplot(d, ggplot2::aes(x = .data$n_cell, y = .data$y, colour = .data$method,
      group = .data$series)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 1.3) +
      ggplot2::geom_blank(data = data.frame(y = c(-1, 1)), ggplot2::aes(y = .data$y),
        inherit.aes = FALSE) +
      ggplot2::geom_line(ggplot2::aes(linetype = .data$tier,
        linewidth = .data$tier)) +
      ggplot2::geom_point(data = d[d$stat == "Median", , drop = FALSE],
        ggplot2::aes(shape = .data$method), size = 1.3) +
      ggplot2::scale_linewidth_manual(values = c(Median = 0.6, `10th and 90th percentiles` = 0.45),
        guide = "none") +
      method_scales +
      ggplot2::scale_x_log10(breaks = sort(unique(long$n_cell)),
        limits = range(long$n_cell), labels = scales::label_comma()) +
      ggplot2::scale_y_continuous(breaks = trans$transform(breaks_err),
        labels = .simCompareThesisPercent(prob * (1 + breaks_err))) +
      ggplot2::facet_grid(prob_response ~ transformation,
        labeller = ggplot2::labeller(prob_response = .simCompareThesisProbLab)) +
      ggplot2::labs(x = if (last) "Number of cells" else NULL, y = NULL) +
      .analysis_theme(grid = "y") + .analysis_thesis_theme() +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = if (first) ggplot2::element_text(size = 7.5) else ggplot2::element_blank(),
        axis.text.x = if (last) ggplot2::element_text(size = 7, angle = 45, hjust = 1) else ggplot2::element_blank(),
        axis.ticks.x = if (last) ggplot2::element_line() else ggplot2::element_blank(),
        plot.margin = ggplot2::margin(1, 2, 1, 2),
        panel.spacing.x = ggplot2::unit(5, "mm")
      )
  }
  rows <- lapply(seq_along(probs), function(i) {
    row_plot(probs[[i]], first = i == 1L, last = i == length(probs))
  })
  legend_plot <- row_plot(probs[[1]], TRUE, TRUE) +
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal") +
    ggplot2::guides(linetype = ggplot2::guide_legend(override.aes = list(colour = "grey20")))
  legend <- cowplot::get_plot_component(legend_plot, "guide-box-bottom", return_all = TRUE)
  if (is.list(legend) && !inherits(legend, "grob")) legend <- legend[[1]]
  stack <- cowplot::plot_grid(plotlist = rows, ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1.12, rep(1, length(rows) - 2L), 1.25)[seq_along(rows)])
  body <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label("Response frequency", angle = 90, size = 8.5),
    stack, nrow = 1, rel_widths = c(0.04, 1)
  )
  cowplot::plot_grid(body, legend, ncol = 1, rel_heights = c(1, 0.06)) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

# Median F1 heat maps: rows for method, columns for transformation; number of
# cells across and response probability up each panel, on one 0-1 scale.
.simCompareThesisF1Plot <- function(cls, pos) {
  tbl <- cls[cls$mean_pos_setting == pos & cls$method %in% .simCompareThesisMethods, , drop = FALSE]
  tbl$transformation <- .analysis_trans_factor(tbl$transformation)
  tbl$method <- factor(tbl$method, .simCompareThesisMethods, .simCompareThesisMethodLabels())
  tbl$cells <- factor(scales::label_comma()(tbl$n_cell),
    scales::label_comma()(sort(unique(tbl$n_cell))))
  tbl$prob <- factor(.simCompareThesisPercent(tbl$prob_response),
    .simCompareThesisPercent(sort(unique(tbl$prob_response))))
  tbl$label <- ifelse(is.finite(tbl$f1_median), sprintf("%.2f", tbl$f1_median), "NA")
  tbl$text_colour <- ifelse(is.finite(tbl$f1_median) & tbl$f1_median < 0.55, "white", "black")
  ggplot2::ggplot(tbl, ggplot2::aes(x = .data$cells, y = .data$prob, fill = .data$f1_median)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, colour = .data$text_colour), size = 2.2) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), name = expression("Median" ~ F[1]),
      na.value = "grey85") +
    ggplot2::facet_grid(method ~ transformation, labeller = ggplot2::labeller(method = ggplot2::label_parsed)) +
    ggplot2::labs(x = "Number of cells", y = "Response probability") +
    .analysis_theme(grid = "none") + .analysis_thesis_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7, angle = 45, hjust = 1),
      axis.line = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(14, "mm"),
      legend.key.height = ggplot2::unit(3, "mm"),
      legend.title = ggplot2::element_text(size = 8, vjust = 0.8)
    )
}

# Threshold distributions: outline histograms of each method's thresholds over
# the stimulated tube's density (grey, scaled so its peak matches the tallest
# threshold bin), with a rug of cells drawn from the simulated positive
# (responding) component, showing where responders lie. One panel per
# response probability and transformation; y values are not shown.
.simCompareThesisThresholdPlot <- function(data, densities, pos, bins = 30L) {
  data <- data[data$mean_pos_setting == pos & data$approach %in% .simCompareThesisMethods &
    is.finite(data$threshold), , drop = FALSE]
  keys <- c("transformation", "prob_response")
  data$transformation <- .analysis_trans_factor(data$transformation)
  densities <- dplyr::semi_join(densities,
    dplyr::distinct(data, dplyr::across(dplyr::all_of(c(keys, "mean_pos")))),
    by = c(keys, "mean_pos"))
  stim <- densities[densities$condition == "stimulated", , drop = FALSE]
  rug <- densities[densities$condition == "positive_rug", , drop = FALSE]
  panels <- dplyr::distinct(data, dplyr::across(dplyr::all_of(keys)))
  steps <- list(); dens <- list()
  for (i in seq_len(nrow(panels))) {
    key <- panels[i, ]
    in_panel <- function(d) d$transformation == key$transformation & d$prob_response == key$prob_response
    d <- data[in_panel(data), , drop = FALSE]
    s <- stim[in_panel(stim), , drop = FALSE]
    breaks <- seq(min(d$threshold), max(d$threshold), length.out = bins + 1L)
    if (diff(range(breaks)) == 0) breaks <- breaks[[1]] + seq(-0.5, 0.5, length.out = bins + 1L)
    hist_max <- 0
    for (m in .simCompareThesisMethods) {
      x <- d$threshold[d$approach == m]
      if (!length(x)) next
      h <- graphics::hist(x, breaks = breaks, plot = FALSE)
      hist_max <- max(hist_max, h$density)
      steps[[length(steps) + 1L]] <- tibble::tibble(key[rep(1, length(breaks) + 1L), ],
        approach = m, x = c(breaks, max(breaks)), y = c(0, h$density, 0))
    }
    if (nrow(s)) {
      s$y <- s$density / max(s$density, na.rm = TRUE) * hist_max
      dens[[length(dens) + 1L]] <- s
    }
  }
  steps <- dplyr::bind_rows(steps)
  steps$approach <- factor(steps$approach, .simCompareThesisMethods)
  dens <- dplyr::bind_rows(dens)
  # Highest response probability first, as in the other thesis figures.
  prob_levels <- sort(unique(data$prob_response), decreasing = TRUE)
  for (nm in c("steps", "dens", "rug")) {
    d <- get(nm)
    d$prob_response <- factor(d$prob_response, prob_levels)
    assign(nm, d)
  }
  labs <- .simCompareThesisMethodLabels()
  ggplot2::ggplot() +
    ggplot2::geom_area(data = dens, ggplot2::aes(x = .data$expression, y = .data$y),
      fill = "grey88", colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_step(data = steps, ggplot2::aes(x = .data$x, y = .data$y,
      colour = .data$approach, linetype = .data$approach), linewidth = 0.5, direction = "hv") +
    ggplot2::geom_rug(data = rug, ggplot2::aes(x = .data$expression), sides = "b",
      colour = "grey15", alpha = 0.6, length = ggplot2::unit(0.06, "npc"), linewidth = 0.25) +
    ggplot2::facet_wrap(ggplot2::vars(.data$prob_response, .data$transformation), ncol = 3,
      scales = "free", labeller = ggplot2::labeller(
        prob_response = function(x) .simCompareThesisProbLab(as.character(x)),
        .multi_line = FALSE)) +
    ggplot2::scale_colour_manual(values = .analysis_method_colours[.simCompareThesisMethods],
      labels = scales::label_parse()(labs), name = "Method") +
    ggplot2::scale_linetype_manual(values = .analysis_method_linetypes[.simCompareThesisMethods],
      labels = scales::label_parse()(labs), name = "Method") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = "Threshold (marker expression)", y = "Threshold density") +
    .analysis_theme(grid = "none") + .analysis_thesis_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 7), legend.position = "bottom")
}
