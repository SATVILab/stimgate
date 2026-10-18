# OMIP-111 (Analysis 14) thesis validation figure: (A) one tube's bivariate
# CD44-by-cytokine panel with the authors' manual gate and each method's gate,
# (B) the same tube's stimulated and control densities in the cytokine's CD44+
# parent with each threshold, and (C) Pearson correlations of net frequencies.
# Source after analysis-plot-style.R, omip-thesis-plots.R and omip111.R.
.omip111ThesisColours <- c(
  StimGate = "#0072B2", `F-beta` = "#D55E00", Tailgate = "#009E73", Manual = "#E41A1C"
)

# Plotmath method labels: F-beta is written as F with subscript beta.
.omip111ThesisMethodLabels <- c(StimGate = "StimGate", `F-beta` = "F[beta]", Tailgate = "Tailgate")

.omip111ThesisAxis <- function(cofactor) {
  raw <- c(-1e3, 0, 1e3, 1e4, 1e5, 1e6)
  list(breaks = asinh(raw / cofactor),
    labels = c("-10³", "0", "10³", "10⁴", "10⁵", "10⁶"))
}

# Gates of one stimulated tube, with their frequencies over every non-naive
# cell and over the cytokine's CD44+ parent (positives: CD44+ cells above the
# gate; manual: the authors' gate). Method gates come from the CD44+ fit.
.omip111ThesisGates <- function(prepared, results, sid, population, marker, methods) {
  input <- prepared$input[[sid]][[population]]
  parent <- input$parent[, marker]
  ref <- prepared$references[prepared$references$sample == sid &
    prepared$references$population == population & prepared$references$marker == marker, ]
  if (nrow(ref) != 1L) stop("Missing author gate for ", sid, " ", marker)
  methodRows <- lapply(methods, function(method) {
    row <- results[results$sampleStim == sid & results$population == population &
      results$marker == marker & results$method == method, ]
    if (nrow(row) != 1L) stop("Missing saved ", method, " result for ", sid, " ", marker)
    pos <- parent & input$expression[, marker] > row$threshold
    if (sum(pos) != row$countStim) stop(method, " gate does not reproduce the saved count.")
    data.frame(gate = method, threshold = row$threshold,
      nonNaive = 100 * sum(pos & input$nonNaive) / sum(input$nonNaive),
      parent = 100 * sum(pos) / sum(parent))
  })
  manual <- input$manual[, marker]
  list(
    gates = rbind(dplyr::bind_rows(methodRows), data.frame(
      gate = "Manual", threshold = NA_real_,
      nonNaive = 100 * sum(manual) / sum(input$nonNaive), parent = 100 * sum(manual) / sum(parent)
    )),
    reference = ref
  )
}

# Gate legend without a title: the manual gate first (drawn thick), then each
# method, with F-beta written as F with subscript beta.
.omip111ThesisGateScale <- function(methods) {
  breaks <- c("Manual", methods)
  labels <- parse(text = ifelse(breaks == "F-beta", "F[beta]", paste0("'", breaks, "'")))
  list(
    ggplot2::scale_colour_manual(values = .omip111ThesisColours[breaks], breaks = breaks,
      labels = labels, name = NULL, drop = FALSE),
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1, override.aes = list(
      linewidth = c(1.1, rep(0.6, length(methods))))))
  )
}

.omip111ThesisBivariate <- function(prepared, results, sid, population, marker, cofactor,
                                    methods = c("StimGate", "F-beta", "Tailgate"), bins = 128L) {
  info <- .omip111ThesisGates(prepared, results, sid, population, marker, methods)
  ref <- info$reference
  input <- prepared$input[[sid]][[population]]
  keep <- input$nonNaive
  cells <- data.frame(x = asinh(input$cd44Raw[keep] / cofactor), y = input$expression[keep, marker])
  key <- (cut(cells$x, bins, labels = FALSE) - 1L) * bins + cut(cells$y, bins, labels = FALSE)
  density <- tabulate(key, bins * bins)[key]
  cells <- cells[order(density), ]
  cells$grey <- grDevices::grey(0.85 - 0.6 * scales::rescale(log10(sort(density))))
  rect <- data.frame(
    xmin = asinh(ref$cd44LowerRaw / cofactor), xmax = asinh(ref$cd44UpperRaw / cofactor),
    ymin = asinh(ref$cytokineLowerRaw / cofactor), ymax = asinh(ref$cytokineUpperRaw / cofactor),
    gate = "Manual"
  )
  lines <- info$gates[info$gates$gate != "Manual", ]
  lines$xmin <- rect$xmin
  lines$xmax <- rect$xmax
  levels <- c(methods, "Manual")
  lines$gate <- factor(lines$gate, levels = levels)
  rect$gate <- factor(rect$gate, levels = levels)
  axis <- .omip111ThesisAxis(cofactor)
  ggplot2::ggplot(cells, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point(colour = cells$grey, size = 0.1, shape = 16) +
    # The thick manual gate first, so method lines that coincide stay visible.
    ggplot2::geom_rect(
      data = rect, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin,
        ymax = .data$ymax, colour = .data$gate),
      inherit.aes = FALSE, fill = NA, linewidth = 1.1, key_glyph = "path"
    ) +
    ggplot2::geom_segment(
      data = lines, ggplot2::aes(x = .data$xmin, xend = .data$xmax, y = .data$threshold,
        yend = .data$threshold, colour = .data$gate),
      inherit.aes = FALSE, linewidth = 0.6, key_glyph = "path"
    ) +
    .omip111ThesisGateScale(methods) +
    ggplot2::scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
    ggplot2::scale_y_continuous(breaks = axis$breaks, labels = axis$labels) +
    ggplot2::labs(x = "CD44", y = .omip111ThesisMarkerLabel(marker)) +
    .analysis_theme(grid = "none") + .analysis_thesis_theme() +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
}

.omip111ThesisMarkerLabel <- function(marker) {
  c(IFNg = "IFNγ", TNF = "TNF", IL2 = "IL-2", IL4_5 = "IL-4/5", IL17A = "IL-17A")[[marker]]
}

.omip111ThesisDensities <- function(prepared, results, sid, sidUns, population, marker, cofactor,
                                    methods = c("StimGate", "F-beta", "Tailgate")) {
  info <- .omip111ThesisGates(prepared, results, sid, population, marker, methods)
  tubes <- c(Stimulated = sid, Control = sidUns)
  dens <- dplyr::bind_rows(lapply(names(tubes), function(tube) {
    x <- .omip111MarkerInput(prepared$input[[tubes[[tube]]]][[population]], marker)$expression[, 1]
    d <- stats::density(x, n = 512)
    data.frame(tube = tube, x = d$x, y = d$y)
  }))
  dens$tube <- factor(dens$tube, levels = names(tubes))
  levels <- c(methods, "Manual")
  lines <- info$gates
  lines$threshold[lines$gate == "Manual"] <- asinh(info$reference$cytokineLowerRaw / cofactor)
  lines$gate <- factor(lines$gate, levels = levels)
  # Draw the thick manual line first, so method lines that coincide stay visible.
  lines <- lines[order(lines$gate != "Manual"), ]
  axis <- .omip111ThesisAxis(cofactor)
  ggplot2::ggplot(dens, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_line(ggplot2::aes(linetype = .data$tube), colour = "grey25", linewidth = 0.5) +
    ggplot2::geom_vline(
      data = lines, ggplot2::aes(xintercept = .data$threshold, colour = .data$gate,
        linewidth = .data$gate == "Manual"),
      key_glyph = "path"
    ) +
    ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.1), guide = "none") +
    ggplot2::scale_linetype_manual(values = c(Stimulated = "solid", Control = "22"), name = NULL) +
    .omip111ThesisGateScale(methods) +
    ggplot2::scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
    ggplot2::labs(x = .omip111ThesisMarkerLabel(marker), y = "Density") +
    .analysis_theme(grid = "none") + .analysis_thesis_theme() +
    ggplot2::theme(
      legend.position = "bottom", legend.direction = "vertical",
      legend.box = "horizontal", legend.box.just = "top", panel.grid = ggplot2::element_blank()
    )
}

# Pearson heatmap with readable cytokine and strain labels, cytokines in the
# paper's order and F-beta written as F with subscript beta.
.omip111ThesisHeatmap <- function(strata, metric = "pearson") {
  order <- c(IFNg = "IFNγ", TNF = "TNF", IL2 = "IL-2", IL4_5 = "IL-4/5", IL17A = "IL-17A")
  strata$marker <- unname(order[strata$marker])
  strata$condition <- c(C57 = "C57BL/6", BALB = "BALB/c")[strata$condition]
  strata <- strata[order(match(strata$condition, c("C57BL/6", "BALB/c")), match(strata$marker, order)), ]
  .omipThesisHeatmap(strata, metric, "standard", .omip111ThesisMethodLabels, parseMethods = TRUE)
}

# The three panels at the thesis text width. `strata` are the per-stratum
# correlations behind panel C (.omipThesisSummary()$strata).
.omip111ThesisFigure <- function(prepared, results, strata, strain, mouse, population, marker,
                                 cofactor) {
  pair <- prepared$samples[prepared$samples$mouse == mouse, ]
  sid <- pair$sample[pair$condition == "stim"]
  sidUns <- pair$sample[pair$condition == "uns"]
  if (length(sid) != 1L || length(sidUns) != 1L || !all(pair$strain == strain)) {
    stop("Expected one stimulated and one control tube for ", mouse, ".")
  }
  a <- .omip111ThesisBivariate(prepared, results, sid, population, marker, cofactor)
  b <- .omip111ThesisDensities(prepared, results, sid, sidUns, population, marker, cofactor)
  c <- .omip111ThesisHeatmap(strata)
  top <- cowplot::plot_grid(a, b, nrow = 1, rel_widths = c(0.45, 0.55), labels = c("A", "B"),
    label_size = 10, align = "h", axis = "tb")
  cowplot::plot_grid(top, c, ncol = 1, rel_heights = c(0.46, 0.54), labels = c("", "C"),
    label_size = 10)
}

# Appendix figure: automated against manual net frequencies for every mouse,
# cytokine and strain, one panel per strain and cytokine (free scales), for
# (A) non-naive CD4 and (B) non-naive CD8 T cells. `scored` comes from
# .omip111FitEval() (combination B: CD44+ positives over every non-naive cell);
# references outside the import tolerance are left out.
.omip111ThesisAgreementPanel <- function(scored, population,
                                         methods = c("StimGate", "F-beta", "Tailgate")) {
  markers <- c(IFNg = "IFNγ", TNF = "TNF", IL2 = "IL-2", IL4_5 = "IL-4/5", IL17A = "IL-17A")
  d <- scored[scored$population == population & scored$referenceAccepted %in% TRUE &
    scored$method %in% methods, ]
  d$strainLabel <- factor(c(C57 = "C57BL/6", BALB = "BALB/c")[d$strain], levels = c("C57BL/6", "BALB/c"))
  d$markerLabel <- factor(markers[d$marker], levels = markers[markers %in% markers[unique(d$marker)]])
  d$method <- factor(d$method, levels = methods)
  ggplot2::ggplot(d, ggplot2::aes(100 * .data$manualBs, 100 * .data$propBs,
    colour = .data$method, shape = .data$method)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.6, alpha = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(ggplot2::vars(.data$strainLabel, .data$markerLabel),
      nrow = 2, scales = "free",
      labeller = ggplot2::labeller(.multi_line = FALSE)
    ) +
    ggplot2::scale_colour_manual(values = .omip111ThesisColours[methods], name = NULL,
      labels = parse(text = .omip111ThesisMethodLabels[methods])) +
    ggplot2::scale_shape_manual(values = c(StimGate = 16, `F-beta` = 17, Tailgate = 15)[methods], name = NULL,
      labels = parse(text = .omip111ThesisMethodLabels[methods])) +
    # Few x breaks, so narrow panels' decimal labels do not overlap.
    ggplot2::scale_x_continuous(n.breaks = 3, labels = function(x) format(x, drop0trailing = TRUE, trim = TRUE)) +
    ggplot2::scale_y_continuous(n.breaks = 4, labels = function(x) format(x, drop0trailing = TRUE, trim = TRUE)) +
    ggplot2::labs(
      x = "Background-subtracted frequency (manual gating, %)",
      y = "Background-subtracted frequency\n(automated gating, %)"
    ) +
    .analysis_theme(grid = "xy") + .analysis_thesis_theme() +
    ggplot2::theme(legend.position = "bottom", panel.spacing = ggplot2::unit(2, "mm"))
}

.omip111ThesisAgreementFigure <- function(scored) {
  a <- .omip111ThesisAgreementPanel(scored, "CD4")
  b <- .omip111ThesisAgreementPanel(scored, "CD8")
  legend <- cowplot::get_plot_component(a, "guide-box-bottom", return_all = TRUE)
  if (is.list(legend) && !inherits(legend, "grob")) legend <- legend[[1]]
  noLegend <- ggplot2::theme(legend.position = "none")
  cowplot::plot_grid(
    a + noLegend, b + noLegend, legend,
    ncol = 1, rel_heights = c(1, 1, 0.08), labels = c("A", "B", ""), label_size = 10
  )
}
