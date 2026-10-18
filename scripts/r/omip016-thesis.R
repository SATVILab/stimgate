# Adapt the existing, count-audited OMIP-016 scores; no gates are estimated here.
.omip016ThesisRows <- function(scored) {
  required <- c("method", "cyt", "stim", "freq_bs_man", "freq_bs", "f1",
    "tp", "fp", "fn", "tn", "nCellStim", "countStim", "nManualStim")
  if (!all(required %in% names(scored))) {
    stop("OMIP-016 thesis rows require the complete existing score table.")
  }
  rows <- scored[as.character(scored$method) != "manual", , drop = FALSE]
  rows$method <- as.character(rows$method)
  counts <- as.matrix(rows[, c("tp", "fp", "fn", "tn")])
  defined <- rowSums(is.finite(counts)) == 4L
  if (any(!is.finite(as.matrix(rows[defined, c("nCellStim", "countStim", "nManualStim")]))) ||
      any(counts[defined, , drop = FALSE] < 0) ||
      any(rowSums(counts)[defined] != rows$nCellStim[defined]) ||
      any((rows$tp + rows$fp)[defined] != rows$countStim[defined]) ||
      any((rows$tp + rows$fn)[defined] != rows$nManualStim[defined])) {
    stop("OMIP-016 thesis confusion counts do not reproduce saved cell counts.")
  }
  expected <- .omipClassificationMetrics(rows)$f1
  if (any(!defined & is.finite(rows$f1)) ||
      any(is.finite(expected) != is.finite(rows$f1)) ||
      any(abs(expected[is.finite(expected)] - rows$f1[is.finite(expected)]) > 1e-12)) {
    stop("OMIP-016 thesis F1 does not reproduce its confusion counts.")
  }
  rows$population <- "CD4 T cells"
  rows$marker <- as.character(rows$cyt)
  rows$condition <- as.character(rows$stim)
  rows$reference <- rows$freq_bs_man
  rows$estimate <- rows$freq_bs
  rows$scatterEligible <- TRUE
  rows
}

# Thesis figure: for each stimulation (rows), the correlation across cytokines
# between each method's and the manual net frequencies (one tile per method,
# stacked vertically), beside a 1x3 scatter of automated against manual net
# frequencies (one panel per method, one point per cytokine). The single
# macaque gives one value per cytokine, so cytokines are the points.
.omip016ThesisMethods <- c(stimgate = "StimGate", fbeta = "F[beta]", tailgate = "Tailgate")
.omip016ThesisMethodColours <- c(stimgate = "#0072B2", fbeta = "#D55E00", tailgate = "#009E73")
.omip016ThesisMarkers <- c(CD154 = "CD154", MIP1b = "MIP-1β", IFNg = "IFNγ", IL2 = "IL-2", TNFa = "TNF")
.omip016ThesisStims <- c(gag = "Gag peptides", seb = "SEB")

.omip016ThesisCorrelations <- function(rows) {
  out <- expand.grid(method = names(.omip016ThesisMethods), condition = names(.omip016ThesisStims),
    stringsAsFactors = FALSE)
  stats <- lapply(seq_len(nrow(out)), function(i) {
    d <- rows[rows$method == out$method[[i]] & rows$condition == out$condition[[i]], ]
    ok <- is.finite(d$estimate) & is.finite(d$reference)
    x <- d$reference[ok]
    y <- d$estimate[ok]
    defined <- sum(ok) >= 3L && stats::sd(x) > 0 && stats::sd(y) > 0
    data.frame(n = sum(ok), n_intended = nrow(d),
      pearson = if (defined) stats::cor(x, y) else NA_real_,
      ccc = if (defined) .acsCytofValidationCcc(x, y) else NA_real_)
  })
  cbind(out, dplyr::bind_rows(stats))
}

.omip016ThesisTiles <- function(correlations, condition, metric) {
  d <- correlations[correlations$condition == condition, ]
  d$method <- factor(d$method, levels = rev(names(.omip016ThesisMethods)))
  d$value <- d[[metric]]
  d$label <- ifelse(is.finite(d$value), sprintf("%.2f", d$value), "NA")
  d$stim <- .omip016ThesisStims[[condition]]
  colours <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  ggplot2::ggplot(d, ggplot2::aes(x = 1, y = .data$method, fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label,
      colour = ifelse(abs(.data$value) > 0.6 & is.finite(.data$value), "white", "black")), size = 2.8) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradientn(colours = colours, limits = c(-1, 1), na.value = "white",
      breaks = seq(-1, 1, 0.5),
      name = c(pearson = "Pearson correlation", ccc = "Concordance correlation")[[metric]],
      guide = ggplot2::guide_colourbar(title.position = "top",
        barwidth = ggplot2::unit(38, "mm"), barheight = ggplot2::unit(3, "mm"))) +
    ggplot2::scale_y_discrete(labels = parse(text = .omip016ThesisMethods[rev(names(.omip016ThesisMethods))])) +
    ggplot2::facet_wrap(ggplot2::vars(.data$stim)) +
    ggplot2::labs(x = NULL, y = NULL) +
    .analysis_theme(grid = "none") + .analysis_thesis_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(), legend.position = "bottom")
}

.omip016ThesisScatter <- function(rows, condition) {
  d <- rows[rows$condition == condition & rows$method %in% names(.omip016ThesisMethods), ]
  d$method <- factor(d$method, levels = names(.omip016ThesisMethods),
    labels = .omip016ThesisMethods)
  d$cytokine <- factor(.omip016ThesisMarkers[d$marker], levels = .omip016ThesisMarkers)
  # One range for all methods of a stimulation, including the identity line.
  range <- range(c(0, d$reference, d$estimate), na.rm = TRUE)
  range <- range + c(-1, 1) * diff(range) * 0.05
  ggplot2::ggplot(d, ggplot2::aes(.data$reference, .data$estimate,
    colour = .data$cytokine, shape = .data$cytokine)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(size = 2, alpha = 0.8, na.rm = TRUE) +
    ggplot2::facet_wrap(ggplot2::vars(.data$method), nrow = 1, labeller = ggplot2::label_parsed) +
    ggplot2::scale_colour_manual(values = c("#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"),
      name = NULL, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 18, 8), name = NULL, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = range, ylim = range) +
    ggplot2::labs(x = "Background-subtracted frequency (manual gating, %)",
      y = "Background-subtracted frequency\n(automated gating, %)") +
    .analysis_theme(grid = "xy") + .analysis_thesis_theme() +
    ggplot2::theme(legend.position = "bottom", panel.spacing = ggplot2::unit(3, "mm"))
}

.omip016ThesisStimFigure <- function(rows, metric = "ccc") {
  correlations <- .omip016ThesisCorrelations(rows)
  noLegend <- ggplot2::theme(legend.position = "none")
  panels <- lapply(names(.omip016ThesisStims), function(condition) {
    list(.omip016ThesisTiles(correlations, condition, metric) + noLegend,
      .omip016ThesisScatter(rows, condition) + noLegend)
  })
  legendOf <- function(p) {
    legend <- cowplot::get_plot_component(p, "guide-box-bottom", return_all = TRUE)
    if (is.list(legend) && !inherits(legend, "grob")) legend <- legend[[1]]
    legend
  }
  legends <- cowplot::plot_grid(
    legendOf(.omip016ThesisTiles(correlations, "gag", metric)),
    legendOf(.omip016ThesisScatter(rows, "gag")),
    nrow = 1, rel_widths = c(0.45, 0.55)
  )
  grid <- cowplot::plot_grid(
    panels[[1]][[1]], panels[[1]][[2]], panels[[2]][[1]], panels[[2]][[2]],
    ncol = 2, rel_widths = c(0.24, 0.76), labels = c("A", "B", "C", "D"), label_size = 10,
    align = "h", axis = "tb"
  )
  structure(cowplot::plot_grid(grid, legends, ncol = 1, rel_heights = c(1, 0.13)),
    correlations = correlations)
}
