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
        limits = range(long$n_cell), labels = .analysis_label_cells) +
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

# F1 heat maps (median, or 10th percentile with stat = "q10"): rows for method,
# columns for transformation; number of cells across and response probability
# up each panel, on one 0-1 scale.
.simCompareThesisF1Plot <- function(cls, pos, stat = c("median", "q10")) {
  stat <- match.arg(stat)
  col <- paste0("f1_", stat)
  fill_name <- if (stat == "median") expression("Median" ~ F[1]) else expression("10th percentile" ~ F[1])
  tbl <- cls[cls$mean_pos_setting == pos & cls$method %in% .simCompareThesisMethods, , drop = FALSE]
  tbl <- .analysis_method_data(tbl)
  tbl$transformation <- .analysis_trans_factor(tbl$transformation)
  tbl$method <- factor(tbl$method, .simCompareThesisMethods, .simCompareThesisMethodLabels())
  tbl$cells <- factor(.analysis_label_cells(tbl$n_cell),
    .analysis_label_cells(sort(unique(tbl$n_cell))))
  tbl$prob <- factor(.simCompareThesisPercent(tbl$prob_response),
    .simCompareThesisPercent(sort(unique(tbl$prob_response))))
  tbl$value <- tbl[[col]]
  tbl$label <- ifelse(is.finite(tbl$value), sprintf("%.2f", tbl$value), "NA")
  tbl$text_colour <- ifelse(is.finite(tbl$value) & tbl$value < 0.55, "white", "#333333")
  ggplot2::ggplot(tbl, ggplot2::aes(x = .data$cells, y = .data$prob, fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, colour = .data$text_colour), size = 2.2) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), name = fill_name,
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
# response probability (rows) and `cols` (columns: transformation, or parent
# cell count when the caller has kept one transformation); y values are not
# shown. Thresholds are pooled over whatever cell counts `data` holds; with
# `split_cells`, each method gets one outline per cell count (line type).
# Panels without data are left blank so columns line up.
.simCompareThesisThresholdPlot <- function(
  data, densities, pos, bins = 30L, cols = c("transformation", "n_cell"),
  split_cells = FALSE
) {
  cols <- match.arg(cols)
  data <- data[data$mean_pos_setting == pos & data$approach %in% .simCompareThesisMethods &
    is.finite(data$threshold), , drop = FALSE]
  data <- .analysis_method_data(data)
  data$transformation <- droplevels(.analysis_trans_factor(data$transformation))
  densities$transformation <- .analysis_trans_factor(densities$transformation)
  dens_keys <- c("transformation", "prob_response")
  densities <- dplyr::semi_join(densities,
    dplyr::distinct(data, dplyr::across(dplyr::all_of(c(dens_keys, "mean_pos")))),
    by = c(dens_keys, "mean_pos"))
  stim <- densities[densities$condition == "stimulated", , drop = FALSE]
  rug <- densities[densities$condition == "positive_rug", , drop = FALSE]
  panel_keys <- unique(c(cols, "prob_response", "transformation"))
  panels <- dplyr::distinct(data, dplyr::across(dplyr::all_of(panel_keys)))
  cell_levels <- sort(unique(data$n_cell))
  groups <- if (split_cells) {
    expand.grid(approach = .simCompareThesisMethods, n_cell = cell_levels, stringsAsFactors = FALSE)
  } else {
    data.frame(approach = .simCompareThesisMethods, n_cell = NA_real_, stringsAsFactors = FALSE)
  }
  steps <- list(); dens <- list(); rugs <- list()
  for (i in seq_len(nrow(panels))) {
    key <- panels[i, ]
    in_panel <- function(d) {
      ok <- d$prob_response == key$prob_response &
        as.character(d$transformation) == as.character(key$transformation)
      if (identical(cols, "n_cell")) ok <- ok & d$n_cell == key$n_cell
      ok
    }
    in_dens <- function(d) {
      as.character(d$transformation) == as.character(key$transformation) &
        d$prob_response == key$prob_response
    }
    d <- data[in_panel(data), , drop = FALSE]
    s <- stim[in_dens(stim), , drop = FALSE]
    breaks <- seq(min(d$threshold), max(d$threshold), length.out = bins + 1L)
    if (diff(range(breaks)) == 0) breaks <- breaks[[1]] + seq(-0.5, 0.5, length.out = bins + 1L)
    hist_max <- 0
    for (g in seq_len(nrow(groups))) {
      sel <- d$approach == groups$approach[[g]]
      if (split_cells) sel <- sel & d$n_cell == groups$n_cell[[g]]
      x <- d$threshold[sel]
      if (!length(x)) next
      h <- graphics::hist(x, breaks = breaks, plot = FALSE)
      hist_max <- max(hist_max, h$density)
      steps[[length(steps) + 1L]] <- tibble::tibble(key[rep(1, length(breaks) + 1L), ],
        approach = groups$approach[[g]], cells = groups$n_cell[[g]],
        x = c(breaks, max(breaks)), y = c(0, h$density, 0))
    }
    if (nrow(s)) {
      s$y <- s$density / max(s$density, na.rm = TRUE) * hist_max
      dens[[length(dens) + 1L]] <- tibble::tibble(key, expression = s$expression, y = s$y)
    }
    r <- rug[in_dens(rug), , drop = FALSE]
    if (nrow(r)) rugs[[length(rugs) + 1L]] <- tibble::tibble(key, expression = r$expression)
  }
  steps <- dplyr::bind_rows(steps)
  steps$approach <- factor(steps$approach, .simCompareThesisMethods)
  dens <- dplyr::bind_rows(dens)
  rug <- dplyr::bind_rows(rugs)
  # Highest response probability first, as in the other thesis figures; keep
  # every level so empty panels stay in place.
  prob_levels <- sort(unique(data$prob_response), decreasing = TRUE)
  col_levels <- if (identical(cols, "n_cell")) cell_levels else levels(data$transformation)
  for (nm in c("steps", "dens", "rug")) {
    d <- get(nm)
    d$prob_response <- factor(d$prob_response, prob_levels)
    d$column <- factor(as.character(d[[cols]]), as.character(col_levels))
    assign(nm, d)
  }
  col_lab <- if (identical(cols, "n_cell")) {
    function(x) paste(.analysis_label_cells(as.numeric(x)), "cells")
  } else {
    as.character
  }
  labs <- .simCompareThesisMethodLabels()
  line_aes <- if (split_cells) {
    steps$cells <- factor(.analysis_label_cells(steps$cells), .analysis_label_cells(cell_levels))
    ggplot2::aes(x = .data$x, y = .data$y, colour = .data$approach,
      linetype = .data$cells, group = interaction(.data$approach, .data$cells))
  } else {
    ggplot2::aes(x = .data$x, y = .data$y, colour = .data$approach, linetype = .data$approach)
  }
  line_scale <- if (split_cells) {
    ggplot2::scale_linetype_manual(
      values = c("solid", "42", "22", "12", "11")[seq_along(cell_levels)],
      name = "Number of cells")
  } else {
    ggplot2::scale_linetype_manual(values = .analysis_method_linetypes[.simCompareThesisMethods],
      labels = scales::label_parse()(labs), name = "Method")
  }
  ggplot2::ggplot() +
    ggplot2::geom_area(data = dens, ggplot2::aes(x = .data$expression, y = .data$y),
      fill = "grey88", colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_step(data = steps, line_aes, linewidth = 0.5, direction = "hv") +
    ggplot2::geom_rug(data = rug, ggplot2::aes(x = .data$expression), sides = "b",
      colour = "grey15", alpha = 0.6, length = ggplot2::unit(0.06, "npc"), linewidth = 0.25) +
    ggplot2::facet_wrap(ggplot2::vars(.data$prob_response, .data$column), ncol = length(col_levels),
      scales = "free", drop = FALSE, labeller = ggplot2::labeller(
        prob_response = function(x) .simCompareThesisProbLab(as.character(x)),
        column = col_lab, .multi_line = FALSE)) +
    ggplot2::scale_colour_manual(values = .analysis_method_colours[.simCompareThesisMethods],
      labels = scales::label_parse()(labs), name = "Method") +
    line_scale +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = "Threshold (marker expression)", y = "Threshold density") +
    .analysis_theme(grid = "none") + .analysis_thesis_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 7), legend.position = "bottom",
      legend.box = "vertical")
}

# Height (mm) for a threshold figure with `n_rows` rows of panels.
.simCompareThesisThresholdHeight <- function(n_rows) min(215, 45 + 34 * n_rows)
