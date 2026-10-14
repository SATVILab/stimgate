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

# Percentile figures: each method's median (points joined by lines) with the
# 10th to 90th percentile range as error bars, over the samples of each
# scenario, dodged by method so single points still show their range. Columns
# are transformations. With rows = "prob", rows are response probabilities and
# the x axis is the number of cells; with rows = "cells", rows are numbers of
# cells and the x axis is the response probability.
# - measure = "frequency": `tbl` holds percentiles of the relative error of the
#   background-subtracted frequency (q10, median, q90). With rows = "prob" the
#   y axis is the background-subtracted frequency in percent (truth: grey line,
#   y range free per row); with rows = "cells" it is the estimated over the
#   true frequency (1 = correct), since the response probabilities differ by
#   orders of magnitude.
# - measure = "f1": `tbl` holds f1_q10, f1_median and f1_q90; the y axis is F1
#   on 0-1.
.simCompareThesisPctPlot <- function(tbl, pos, measure = c("frequency", "f1"),
                                     rows = c("prob", "cells")) {
  measure <- match.arg(measure)
  rows <- match.arg(rows)
  tbl <- tbl[tbl$mean_pos_setting == pos & tbl$method %in% .simCompareThesisMethods, , drop = FALSE]
  if (identical(measure, "f1")) {
    tbl$q10 <- tbl$f1_q10; tbl$median <- tbl$f1_median; tbl$q90 <- tbl$f1_q90
  }
  conv <- if (identical(measure, "f1")) {
    function(v, prob) v
  } else if (identical(rows, "prob")) {
    function(v, prob) 100 * prob * (1 + v)
  } else {
    function(v, prob) 1 + v
  }
  d <- data.frame(
    method = factor(tbl$method, .simCompareThesisMethods),
    transformation = .analysis_trans_factor(tbl$transformation),
    prob_response = tbl$prob_response, n_cell = tbl$n_cell,
    lo = conv(tbl$q10, tbl$prob_response), mid = conv(tbl$median, tbl$prob_response),
    hi = conv(tbl$q90, tbl$prob_response)
  )
  probs <- sort(unique(d$prob_response), decreasing = TRUE)
  cells <- sort(unique(d$n_cell))
  if (identical(rows, "prob")) {
    d$row <- factor(d$prob_response, probs, .simCompareThesisProbLab(probs))
    d$x <- d$n_cell
    x_scale <- ggplot2::scale_x_log10(breaks = cells, labels = .analysis_label_cells,
      expand = ggplot2::expansion(mult = 0.08))
    x_lab <- "Number of cells"
  } else {
    d$row <- factor(d$n_cell, cells, paste(.analysis_label_cells(cells), "cells"))
    d$x <- d$prob_response
    x_scale <- ggplot2::scale_x_log10(breaks = rev(probs), labels = .simCompareThesisPercent,
      expand = ggplot2::expansion(mult = 0.08))
    x_lab <- "Response probability"
  }
  truth <- if (identical(measure, "frequency")) {
    t <- dplyr::distinct(d, .data$row, .data$prob_response, .data$transformation)
    t$y <- if (identical(rows, "prob")) 100 * t$prob_response else 1
    dplyr::distinct(t, .data$row, .data$transformation, .data$y)
  }
  y_lab <- switch(paste(measure, rows),
    "frequency prob" = "Background-subtracted frequency (%)",
    "frequency cells" = "Estimated / true background-subtracted frequency",
    expression(F[1] ~ "score"))
  dodge <- ggplot2::position_dodge(width = 0.18)
  labs <- .simCompareThesisMethodLabels()
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, colour = .data$method, group = .data$method))
  if (!is.null(truth)) {
    p <- p + ggplot2::geom_hline(data = truth, ggplot2::aes(yintercept = .data$y),
      colour = "grey55", linewidth = 1.1)
  }
  p +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), width = 0.12,
      linewidth = 0.4, position = dodge) +
    ggplot2::geom_line(ggplot2::aes(y = .data$mid), linewidth = 0.5, position = dodge) +
    ggplot2::geom_point(ggplot2::aes(y = .data$mid, shape = .data$method), size = 1.4,
      position = dodge) +
    ggplot2::facet_grid(ggplot2::vars(.data$row), ggplot2::vars(.data$transformation),
      scales = if (identical(measure, "f1")) "fixed" else "free_y") +
    x_scale +
    (if (identical(measure, "f1")) ggplot2::scale_y_continuous(limits = c(0, 1)) else
      ggplot2::scale_y_continuous(labels = .analysis_label_number)) +
    ggplot2::scale_colour_manual(values = .analysis_method_colours[.simCompareThesisMethods],
      labels = scales::label_parse()(labs), name = "Method") +
    ggplot2::scale_shape_manual(values = .analysis_method_shapes[.simCompareThesisMethods],
      labels = scales::label_parse()(labs), name = "Method") +
    ggplot2::labs(x = x_lab, y = y_lab) +
    .analysis_theme(grid = "y") + .analysis_thesis_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7, angle = 45, hjust = 1),
      legend.position = "bottom", panel.spacing.x = ggplot2::unit(4, "mm")
    )
}

# Estimated background-subtracted frequency, rows for response probability.
.simCompareThesisFreqPlot <- function(pct, pos) {
  .simCompareThesisPctPlot(pct, pos, measure = "frequency", rows = "prob")
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
