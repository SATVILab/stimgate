# Analysis 15 figures in the layout of the OMIP-016 paper's Figure 1E, as
# Analysis 14 redraws OMIP-111's Figures S11/S12 (scripts/r/omip111.R).
#
# Every manual response gate is a rectangle on its marker (x) against a
# second parameter (y): IFN-gamma for CD154 and MIP-1beta, SSC-A for the
# others. The panels show CD4 T cells (the gates' parent) on the GatingSet
# scale, asinh(x / 150) of compensated fluorescence, with SSC-A linear.
# Rectangle sides drawn at or below zero on a compensated axis are extended
# to the panel edge, as the gates are applied (.omip016ExtendVertices()).
#
# Requires scripts/r/omip016-prepare.R and scripts/r/omip016-methods.R.

.omip016FigureMarkers <- c(
  CD154 = "CD154 FITC", MIP1b = "MIP-1β PE", IFNg = "IFN-γ V450",
  TNFa = "TNF-α Ax700", IL2 = "IL-2 APC"
)

.omip016FigureTubes <- c(uns = "Unstimulated", gag = "Gag peptides", seb = "SEB")

.omip016FigureParamLabel <- function(param, markerMap) {
  chnl <- .omip016ParamChannel(param)
  k <- match(chnl, markerMap$channel)
  if (is.na(k) || !nzchar(markerMap$fluorochrome[[k]])) return(chnl)
  marker <- markerMap$marker[[k]]
  if (marker %in% names(.omip016FigureMarkers)) return(.omip016FigureMarkers[[marker]])
  paste(marker, markerMap$fluorochrome[[k]])
}

.omip016FigureLabel <- function(x) {
  vapply(x, function(v) {
    if (!is.finite(v)) return("NA")
    paste0(formatC(v, digits = if (v < 1) 2L else 3L, format = "fg", flag = "#"), "%")
  }, character(1))
}

# A rectangle bound on the GatingSet scale. Bounds at or below zero on a
# compensated axis are extended to the data edge (-Inf here).
.omip016FigureBound <- function(value, param, lower) {
  if (!grepl("^<.*>$", param)) return(value)
  if (lower && value <= .omip016ExtendVal) return(-Inf)
  .omip016Transform(value)
}

# One row per marker: the gate's axes and its rectangle on the GatingSet scale.
.omip016FigureGates <- function(thresholds) {
  thresholds <- thresholds[match(names(.omip016FigureMarkers), thresholds$marker), , drop = FALSE]
  if (anyNA(thresholds$marker)) stop("Manual thresholds lack a Figure 1E marker.")
  data.frame(
    marker = thresholds$marker, channel = thresholds$channel,
    y_param = thresholds$y_param, y_channel = .omip016ParamChannel(thresholds$y_param),
    xmin = thresholds$threshold_trans,
    xmax = .omip016Transform(thresholds$x_max),
    ymin = mapply(.omip016FigureBound, thresholds$y_min, thresholds$y_param, TRUE),
    ymax = mapply(.omip016FigureBound, thresholds$y_max, thresholds$y_param, FALSE),
    stringsAsFactors = FALSE
  )
}

# CD4 T cells of every tube on each panel's axes, with FlowJo-style
# pseudocolour (cell count in the cell's bin), and the manual gate's count.
.omip016FigurePanels <- function(prep, thresholds, tubes = names(.omip016FigureTubes),
                                 bins = 128L) {
  gates <- .omip016FigureGates(thresholds)
  sm <- prep$pre$sampleMap
  cells <- list()
  author <- list()
  for (st in tubes) {
    ind <- sm$ind[sm$stim == st]
    if (length(ind) != 1L) stop("Expected one OMIP-016 tube for ", st, ".")
    ind <- as.integer(ind)
    ex <- .omip016Expr(prep$gs, ind, unique(c(gates$channel, gates$y_channel)))
    lab <- prep$labels[[sm$file[[ind]]]]
    if (nrow(lab) != nrow(ex)) stop("Manual labels do not match the GatingSet cells.")
    for (j in seq_len(nrow(gates))) {
      x <- ex[, gates$channel[[j]]]
      y <- ex[, gates$y_channel[[j]]]
      key <- (cut(x, bins, labels = FALSE) - 1L) * bins + cut(y, bins, labels = FALSE)
      density <- tabulate(key, bins * bins)[key]
      ord <- order(density)
      cells[[length(cells) + 1L]] <- data.frame(
        tube = st, marker = gates$marker[[j]], x = x[ord], y = y[ord], density = density[ord]
      )
      author[[length(author) + 1L]] <- data.frame(
        gates[j, c("marker", "xmin", "xmax", "ymin", "ymax")], tube = st, ind = ind,
        count = sum(lab[[gates$marker[[j]]]]), n_cell = nrow(ex)
      )
    }
  }
  author <- dplyr::bind_rows(author)
  author$frequency <- 100 * author$count / author$n_cell
  tube <- function(x) {
    x$tube <- factor(.omip016FigureTubes[x$tube], levels = .omip016FigureTubes[tubes])
    x
  }
  list(cells = tube(dplyr::bind_rows(cells)), reference = tube(author), gates = gates)
}

# Every method's gate on the stimulated tubes' panels, as a % of CD4 T cells
# above it. The counts must reproduce the scored stimulated counts.
.omip016MethodFigureData <- function(prep, thresholds, scored, methods, bins = 128L) {
  stims <- intersect(names(.omip016FigureTubes), as.character(unique(scored$stim)))
  figureData <- .omip016FigurePanels(prep, thresholds, stims, bins)
  author <- figureData$reference
  scored <- scored[as.character(scored$method) %in% names(methods), , drop = FALSE]
  sm <- prep$pre$sampleMap
  rows <- lapply(seq_len(nrow(author)), function(i) {
    a <- author[i, , drop = FALSE]
    chnl <- figureData$gates$channel[figureData$gates$marker == a$marker]
    x <- .omip016Expr(prep$gs, a$ind, chnl)[, 1]
    dplyr::bind_rows(lapply(names(methods), function(method) {
      s <- scored[as.character(scored$method) == method & scored$ind == a$ind &
        as.character(scored$cyt) == a$marker, , drop = FALSE]
      if (nrow(s) != 1L) stop("Missing scored ", method, " gate for ", a$tube, " ", a$marker, ".")
      count <- if (is.finite(s$gate)) sum(x > s$gate) else NA_integer_
      # Without cytokine-positive refinement, x > gate is the scored count.
      if (is.finite(s$gate) && (!is.finite(s$gateCyt) || s$gateCyt == s$gate) &&
          count != s$countStim) {
        stop(method, " gate does not reproduce the scored count: ", a$tube, " ", a$marker, ".")
      }
      data.frame(
        a[c("tube", "marker", "ymin", "ymax")], method = methods[[method]],
        x = s$gate, frequency = 100 * count / length(x)
      )
    }))
  })
  gates <- dplyr::bind_rows(rows)
  gates$method <- factor(gates$method, levels = unname(methods))
  gates$label <- ifelse(is.finite(gates$x),
    paste0(gates$method, ": ", .omip016FigureLabel(gates$frequency)),
    paste0(gates$method, ": no gate")
  )
  author$label <- paste0("Manual: ", .omip016FigureLabel(author$frequency))
  figureData$reference <- author
  figureData$methods <- gates
  figureData
}

# Compensated fluorescence on the asinh scale (ticks in raw units); scatter
# parameters linear, in thousands.
.omip016FigureScale <- function(param, axis) {
  scale <- if (axis == "x") ggplot2::scale_x_continuous else ggplot2::scale_y_continuous
  if (!grepl("^<.*>$", param)) return(scale(labels = function(v) paste0(v / 1000, "K")))
  raw <- c(0, 1e3, 1e4, 1e5)
  scale(breaks = .omip016Transform(raw), labels = c("0", "10³", "10⁴", "10⁵"))
}

.omip016FigureTheme <- function() {
  list(
    .analysis_theme(grid = "none"),
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey40"),
      axis.text = ggplot2::element_text(size = 7),
      legend.position = "none"
    )
  )
}

# One column of panels (a marker, tubes in rows) with that gate's own axes.
.omip016FigureColumn <- function(figureData, marker, markerMap, layers, strips) {
  g <- figureData$gates[figureData$gates$marker == marker, ]
  cells <- figureData$cells[figureData$cells$marker == marker, ]
  p <- ggplot2::ggplot(cells, ggplot2::aes(.data$x, .data$y)) +
    layers(marker) +
    ggplot2::facet_grid(ggplot2::vars(.data$tube)) +
    .omip016FigureScale(paste0("<", g$channel, ">"), "x") +
    .omip016FigureScale(g$y_param, "y") +
    ggplot2::labs(
      x = .omip016FigureParamLabel(g$channel, markerMap),
      y = .omip016FigureParamLabel(g$y_param, markerMap)
    ) +
    .omip016FigureTheme()
  if (!strips) p <- p + ggplot2::theme(strip.text.y = ggplot2::element_blank())
  p
}

.omip016FigureGrid <- function(figureData, markerMap, layers, legend = NULL) {
  markers <- figureData$gates$marker
  cols <- lapply(seq_along(markers), function(k) {
    .omip016FigureColumn(figureData, markers[[k]], markerMap, layers,
      strips = k == length(markers))
  })
  grid <- cowplot::plot_grid(plotlist = cols, nrow = 1, align = "h", axis = "tb")
  if (is.null(legend)) return(grid)
  cowplot::plot_grid(grid, legend, ncol = 1, rel_heights = c(1, 0.06))
}

# Figure 1E redrawn: pseudocolour CD4 T cells with the authors' rectangles
# and their frequencies (% of CD4 T cells) in red, for every tube.
.omip016AuthorFigure <- function(figureData, markerMap) {
  ref <- figureData$reference
  ref$label <- .omip016FigureLabel(ref$frequency)
  layers <- function(marker) {
    r <- ref[ref$marker == marker, ]
    list(
      ggplot2::geom_point(ggplot2::aes(colour = .data$density), size = 0.12, shape = 16),
      ggplot2::geom_rect(
        data = r, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
        inherit.aes = FALSE, fill = NA, colour = "red", linewidth = 0.4
      ),
      ggplot2::geom_text(
        data = r, ggplot2::aes(x = Inf, y = Inf, label = .data$label),
        inherit.aes = FALSE, colour = "red", hjust = 1.1, vjust = 1.4, size = 2.8
      ),
      ggplot2::scale_colour_gradientn(
        colours = c("#0000FF", "#00BFFF", "#00FF00", "#FFFF00", "#FF0000"),
        trans = "log10", guide = "none"
      )
    )
  }
  .omip016FigureGrid(figureData, markerMap, layers)
}

.omip016FigureGateScales <- function(methods, methodColours) {
  if (!all(methods %in% names(methodColours))) stop("Missing method colours.")
  list(
    ggplot2::scale_colour_manual(values = c(methodColours[methods], Manual = "black"),
      name = NULL, drop = FALSE),
    ggplot2::scale_linetype_manual(
      values = c(stats::setNames(rep("solid", length(methods)), methods), Manual = "22"),
      name = NULL, drop = FALSE
    )
  )
}

.omip016FigureLegend <- function(methods, methodColours) {
  gate <- factor(c(methods, "Manual"), levels = c(methods, "Manual"))
  d <- data.frame(gate = rep(gate, each = 2L), x = rep(0:1, length(gate)))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$x, colour = .data$gate, linetype = .data$gate)) +
    ggplot2::geom_line() +
    .omip016FigureGateScales(methods, methodColours) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom")
  cowplot::get_plot_component(p, "guide-box-bottom", return_all = FALSE)
}

# Per-panel label blocks: the manual gate first, then each method in its colour.
.omip016FigureLabels <- function(author, gates, methods, methodColours) {
  dplyr::bind_rows(
    data.frame(author[c("tube", "marker", "label")], line = 0L, colour = "black"),
    data.frame(gates[c("tube", "marker", "label")],
      line = match(as.character(gates$method), methods),
      colour = unname(methodColours[as.character(gates$method)])
    )
  )
}

# Bivariate panels for the stimulated tubes: grey density shading, the
# authors' rectangle (cutoff dashed black, other sides grey) and each
# method's gate as a vertical line across the rectangle's y interval.
# Labels are stimulated-tube % of CD4 T cells, without background subtraction.
.omip016MethodsFigure <- function(figureData, markerMap, methodColours) {
  gates <- figureData$methods
  author <- figureData$reference
  methods <- levels(gates$method)
  shade <- scales::rescale(log10(figureData$cells$density))
  figureData$cells$grey <- grDevices::grey(0.85 - 0.6 * shade)
  labels <- .omip016FigureLabels(author, gates, methods, methodColours)
  lines <- dplyr::bind_rows(
    data.frame(gates[is.finite(gates$x), c("tube", "marker", "x", "ymin", "ymax")],
      gate = as.character(gates$method[is.finite(gates$x)])),
    data.frame(author[c("tube", "marker", "ymin", "ymax")], x = author$xmin, gate = "Manual")
  )
  lines$gate <- factor(lines$gate, levels = c(methods, "Manual"))
  lines <- lines[order(lines$gate == "Manual"), ]
  layers <- function(marker) {
    a <- author[author$marker == marker, ]
    cells <- figureData$cells[figureData$cells$marker == marker, ]
    l <- lines[lines$marker == marker, ]
    lab <- labels[labels$marker == marker, ]
    list(
      ggplot2::geom_point(colour = cells$grey, size = 0.12, shape = 16),
      ggplot2::geom_segment(
        data = dplyr::bind_rows(
          dplyr::transmute(a, .data$tube, x = .data$xmin, xend = .data$xmax, y = .data$ymin, yend = .data$ymin),
          dplyr::transmute(a, .data$tube, x = .data$xmin, xend = .data$xmax, y = .data$ymax, yend = .data$ymax),
          dplyr::transmute(a, .data$tube, x = .data$xmax, xend = .data$xmax, y = .data$ymin, yend = .data$ymax)
        ),
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
        inherit.aes = FALSE, colour = "grey30", linewidth = 0.35
      ),
      ggplot2::geom_segment(
        data = l,
        ggplot2::aes(x = .data$x, xend = .data$x, y = .data$ymin, yend = .data$ymax,
          colour = .data$gate, linetype = .data$gate),
        inherit.aes = FALSE, linewidth = 0.6
      ),
      .omip016FigureGateScales(methods, methodColours),
      ggplot2::geom_text(
        data = lab,
        ggplot2::aes(x = -Inf, y = Inf, label = .data$label, vjust = 1.3 + 1.25 * .data$line),
        colour = lab$colour, inherit.aes = FALSE, hjust = -0.04, size = 2.3
      )
    )
  }
  .omip016FigureGrid(figureData, markerMap, layers, .omip016FigureLegend(methods, methodColours))
}

# Univariate panels: each marker in the stimulated tubes' CD4 T cells, with
# each method's gate and the authors' cutoff (dashed), labelled as above.
.omip016UnivariateFigure <- function(figureData, methodColours, bins = 80L) {
  gates <- figureData$methods
  author <- figureData$reference
  methods <- levels(gates$method)
  rows <- function(x) {
    x$column <- factor(.omip016FigureMarkers[as.character(x$marker)], levels = .omip016FigureMarkers)
    x
  }
  lines <- rows(dplyr::bind_rows(
    data.frame(gates[is.finite(gates$x), c("tube", "marker", "x")],
      gate = as.character(gates$method[is.finite(gates$x)])),
    data.frame(author[c("tube", "marker")], x = author$xmin, gate = "Manual")
  ))
  lines$gate <- factor(lines$gate, levels = c(methods, "Manual"))
  lines <- lines[order(lines$gate == "Manual"), ]
  labels <- rows(.omip016FigureLabels(author, gates, methods, methodColours))
  # Each marker's own values: the bivariate cells repeat a marker per y axis.
  cells <- rows(figureData$cells)
  raw <- c(0, 1e3, 1e4, 1e5)
  ggplot2::ggplot(cells, ggplot2::aes(.data$x)) +
    ggplot2::geom_histogram(bins = bins, fill = "grey60", colour = NA) +
    ggplot2::geom_blank(data = lines, ggplot2::aes(x = .data$x + 0.3)) +
    ggplot2::geom_vline(
      data = lines, ggplot2::aes(xintercept = .data$x, colour = .data$gate, linetype = .data$gate),
      linewidth = 0.6
    ) +
    .omip016FigureGateScales(methods, methodColours) +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = -Inf, y = Inf, label = .data$label, vjust = 1.3 + 1.25 * .data$line),
      colour = labels$colour, inherit.aes = FALSE, hjust = -0.04, size = 2.3
    ) +
    ggplot2::facet_grid(ggplot2::vars(.data$tube), ggplot2::vars(.data$column), scales = "free") +
    ggplot2::scale_x_continuous(breaks = .omip016Transform(raw),
      labels = c("0", "10³", "10⁴", "10⁵")) +
    ggplot2::scale_y_continuous(trans = "sqrt") +
    ggplot2::labs(x = NULL, y = "CD4 T cells (square-root scale)") +
    .omip016FigureTheme() +
    ggplot2::theme(legend.position = "bottom")
}
