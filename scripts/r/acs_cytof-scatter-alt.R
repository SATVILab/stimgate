# Alternative ACS frequency scatter plots for the thesis (Analysis 10): every
# population in one figure, one row per population (labelled once, on the
# right) and one column per cytokine. Source after analysis-plot-style.R and
# acs_cytof-plot_cyt.R.
#
# Two families, each with three line-of-best-fit variants (`fit`):
# - all methods together (`method = NULL`): colour is the method (shared
#   method colours), shape is the stimulation;
# - one method (`method = "stimgate"` etc.): colour and shape are the
#   stimulation (the ACS stimulation colours).
# `fit = "stim"` draws an ordinary least-squares line per stimulation (line
# type by stimulation), `fit = "panel"` one line through every point of the
# panel (per method when methods share a panel) and `fit = "none"` no lines.
#
# Each panel is square and has the same range on both axes (its own range, so
# panels still differ), so y = x is the diagonal in every panel. Both axes use
# a signed square root, which spreads the many points near zero, keeps
# negative background-subtracted frequencies, stops single large values from
# squashing the rest and leaves y = x straight. Lines are fitted on the
# original (percentage) scale and drawn through many points, so they curve as
# the axes do; they are cut off at the panel's range.
#
# Points are fillable shapes with an opaque outline and a translucent fill of
# the same colour: the outline keeps each point visible where points overlap,
# and the fill shows where they pile up. Methods are drawn Tailgate first and
# StimGate last (on top).

# Methods compared in the thesis, mapped to their shared colour/label keys.
.acsCytofScatterAltMethods <- function() {
  c(stimgate = "stimgate", fbeta_default = "fbeta", tailgate = "tailgate")
}

# Drawing order of the methods: later methods are drawn on top.
.acsCytofScatterAltDrawOrder <- function() c("tailgate", "fbeta_default", "stimgate")

.acsCytofScatterAltStims <- function() c("mtb", "p1", "p4", "ebv")

# Fillable versions of the stacked scatter shapes (triangle, square, diamond
# for the plus, circle).
.acsCytofScatterAltShapes <- function() c(mtb = 24, p1 = 22, p4 = 23, ebv = 21)

.acsCytofScatterAltLinetypes <- function() {
  c(mtb = "solid", p1 = "42", p4 = "22", ebv = "1343")
}

.acsCytofScatterAltPopulations <- function() {
  c("CD4 T cells", "CD8 T cells", "TCRgd T cells", "B cells", "NK cells")
}

.acsCytofScatterAltCytokines <- function() c("IFNg", "IL2", "TNF", "IL17", "IL22", "IL6")

# Plotted rows: the selected methods, with failed (or unknown-status)
# thresholds and non-finite frequencies removed, as in the stacked scatter.
# Rows are sorted so that methods are drawn in `.acsCytofScatterAltDrawOrder()`.
.acsCytofScatterAltData <- function(comparisonTbl, methods) {
  .acsCytofValidationValidateComparisonTable(comparisonTbl)
  rows <- comparisonTbl |>
    dplyr::filter(as.character(.data$method) %in% .env$methods)
  failed <- if ("thresholdFailed" %in% names(rows)) {
    unname(rows$thresholdFailed)
  } else {
    rep(FALSE, nrow(rows))
  }
  keep <- !is.na(failed) & !failed &
    is.finite(rows$freq_bs_man) & is.finite(rows$freq_bs_auto)
  rows <- rows[keep, , drop = FALSE]
  populations <- .acsCytofScatterAltPopulations()
  available <- unique(as.character(comparisonTbl$pop))
  cytokines <- .acsCytofScatterAltCytokines()
  stims <- .acsCytofScatterAltStims()
  out <- rows |>
    dplyr::transmute(
      method = factor(as.character(.data$method), levels = methods),
      pop = factor(
        as.character(.data$pop),
        levels = c(intersect(populations, available), setdiff(available, populations))
      ),
      cyt = factor(
        as.character(.data$cyt),
        levels = c(cytokines, setdiff(unique(as.character(comparisonTbl$cyt)), cytokines))
      ),
      stim = factor(
        as.character(.data$stim),
        levels = c(stims, setdiff(unique(as.character(comparisonTbl$stim)), stims))
      ),
      SampleID = .data$SampleID,
      freq_bs_man = .data$freq_bs_man,
      freq_bs_auto = .data$freq_bs_auto
    )
  out[order(match(as.character(out$method), .acsCytofScatterAltDrawOrder())), , drop = FALSE]
}

# Signed square root, sign(x) * sqrt(|x|), for frequencies that can be
# negative. Breaks are zero, the largest one-digit round value within the
# range and a quarter of the range's maximum rounded down (halfway along the
# square-root axis), plus a negative break when negative values reach that far.
.acsCytofSignedSqrtBreaks <- function(limits) {
  limits <- limits[is.finite(limits)]
  if (length(limits) < 2L || diff(range(limits)) <= 0) return(limits)
  roundDown <- function(x) {
    out <- numeric(length(x))
    pos <- x > 0
    f <- 10^floor(log10(x[pos]))
    out[pos] <- signif(floor(x[pos] / f + 1e-9) * f, 6)
    out
  }
  lo <- min(limits)
  hi <- max(limits)
  top <- roundDown(max(hi, 0))
  breaks <- c(0, roundDown(top / 4), top)
  if (lo < 0 && -lo >= max(hi, 0) / 4) breaks <- c(-roundDown(-lo), breaks)
  breaks <- sort(unique(breaks))
  breaks[breaks >= lo & breaks <= hi]
}

.acsCytofSignedSqrtTrans <- function() {
  scales::trans_new(
    "signed_sqrt",
    transform = function(x) sign(x) * sqrt(abs(x)),
    inverse = function(x) sign(x) * x^2,
    breaks = .acsCytofSignedSqrtBreaks,
    domain = c(-Inf, Inf)
  )
}

# The common x and y range of each population/cytokine panel: from the lower
# of zero and the smallest frequency on either axis to the largest.
.acsCytofScatterAltLimits <- function(plotTbl) {
  plotTbl |>
    dplyr::group_by(.data$pop, .data$cyt) |>
    dplyr::summarise(
      lo = min(0, .data$freq_bs_man, .data$freq_bs_auto),
      hi = max(0, .data$freq_bs_man, .data$freq_bs_auto),
      .groups = "drop"
    ) |>
    dplyr::mutate(hi = ifelse(.data$hi > .data$lo, .data$hi, .data$lo + 0.01)) |>
    as.data.frame()
}

# Least-squares line of automated on manual frequency (original scale) within
# each group, as `n` points over the group's manual range (evenly spaced on
# the signed square-root axis) so the line can be drawn on transformed axes. Groups with fewer than `minPoints` points or
# fewer than `minDistinct` distinct manual frequencies get no line. With
# `limits` (from `.acsCytofScatterAltLimits()`), points outside the panel's
# range are dropped, so a line stops at the panel edge.
.acsCytofScatterAltFits <- function(
  plotTbl, by, limits = NULL, n = 100L, minPoints = 4L, minDistinct = 3L
) {
  groups <- split(plotTbl, plotTbl[by], drop = TRUE)
  fits <- lapply(groups, function(g) {
    x <- g$freq_bs_man
    if (nrow(g) < minPoints || length(unique(x)) < minDistinct) return(NULL)
    coef <- stats::coef(stats::lm(freq_bs_auto ~ freq_bs_man, data = g))
    # Evenly spaced on the signed square-root axis, where the line is drawn.
    root <- seq(sign(min(x)) * sqrt(abs(min(x))), sign(max(x)) * sqrt(abs(max(x))), length.out = n)
    xs <- sign(root) * root^2
    out <- g[rep(1L, n), by, drop = FALSE]
    out$freq_bs_man <- xs
    out$freq_bs_auto <- unname(coef[[1L]] + coef[[2L]] * xs)
    out$intercept <- unname(coef[[1L]])
    out$slope <- unname(coef[[2L]])
    out
  })
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) {
    out <- plotTbl[0, by, drop = FALSE]
    out$freq_bs_man <- out$freq_bs_auto <- out$intercept <- out$slope <- numeric()
    return(out)
  }
  out <- do.call(rbind, fits)
  rownames(out) <- NULL
  if (!is.null(limits)) {
    lim <- merge(out, limits, by = c("pop", "cyt"), sort = FALSE)
    out <- lim[lim$freq_bs_auto >= lim$lo & lim$freq_bs_auto <= lim$hi, names(out), drop = FALSE]
    rownames(out) <- NULL
  }
  out
}

# Plotmath method labels for legends (F-beta as F with subscript beta).
.acsCytofScatterAltMethodLabels <- function(methods) {
  keys <- .acsCytofScatterAltMethods()[methods]
  labels <- c(stimgate = "'StimGate'", fbeta = "F[beta]", tailgate = "'Tailgate'")
  parse(text = unname(labels[keys]))
}

# One population's row: square panels, one per cytokine, each with its own
# common x/y range.
.acsCytofScatterAltRow <- function(plotTbl, fitTbl, limTbl, byMethod, fit, methods) {
  stimLabels <- .acsCytofValidationStimLabels()
  stimLevels <- levels(plotTbl$stim)
  stimName <- "Stimulation"
  fillAlpha <- 0.35
  corners <- rbind(
    data.frame(limTbl[c("pop", "cyt")], freq_bs_man = limTbl$lo, freq_bs_auto = limTbl$lo),
    data.frame(limTbl[c("pop", "cyt")], freq_bs_man = limTbl$hi, freq_bs_auto = limTbl$hi)
  )
  xy <- ggplot2::aes(x = .data$freq_bs_man, y = .data$freq_bs_auto)
  pointAes <- if (byMethod) {
    ggplot2::aes(
      x = .data$freq_bs_man, y = .data$freq_bs_auto,
      colour = .data$method, fill = .data$method, shape = .data$stim
    )
  } else {
    ggplot2::aes(
      x = .data$freq_bs_man, y = .data$freq_bs_auto,
      colour = .data$stim, fill = .data$stim, shape = .data$stim
    )
  }
  sqrtTrans <- .acsCytofSignedSqrtTrans()

  plot <- ggplot2::ggplot(plotTbl) +
    .analysis_theme(base_size = 8) +
    .analysis_thesis_theme() +
    ggplot2::geom_blank(data = corners, mapping = xy) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey65", linewidth = 0.4)

  if (fit != "none" && nrow(fitTbl)) {
    colourVar <- if (byMethod) "method" else if (fit == "stim") "stim"
    lineArgs <- list(x = quote(freq_bs_man), y = quote(freq_bs_auto), group = quote(.line))
    if (!is.null(colourVar)) lineArgs$colour <- as.name(colourVar)
    if (fit == "stim") lineArgs$linetype <- quote(stim)
    lineAes <- do.call(ggplot2::aes, lineArgs)
    plot <- plot + if (is.null(colourVar)) {
      # One method, one line per panel: black, apart from the stimulations and
      # from the lighter y = x line.
      ggplot2::geom_line(data = fitTbl, mapping = lineAes, linewidth = 0.7, colour = "black")
    } else {
      ggplot2::geom_line(data = fitTbl, mapping = lineAes, linewidth = 0.55)
    }
  }

  plot <- plot +
    ggplot2::geom_point(mapping = pointAes, size = 1.4, stroke = 0.45) +
    ggplot2::facet_wrap(ggplot2::vars(cyt), nrow = 1L, scales = "free", drop = FALSE) +
    ggplot2::scale_x_continuous(
      transform = sqrtTrans, labels = .analysis_label_number,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    ) +
    ggplot2::scale_y_continuous(
      transform = sqrtTrans, labels = .analysis_label_number,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    ) +
    ggplot2::scale_shape_manual(
      values = .acsCytofScatterAltShapes(), labels = stimLabels,
      limits = stimLevels, name = stimName
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      aspect.ratio = 1,
      strip.text = ggplot2::element_text(size = 7, margin = ggplot2::margin(1.5, 1, 1.5, 1)),
      axis.text = ggplot2::element_text(size = 6),
      axis.text.x = ggplot2::element_text(angle = 0),
      panel.spacing = ggplot2::unit(1.5, "mm"),
      legend.box = "vertical",
      legend.spacing.y = ggplot2::unit(0, "mm"),
      plot.margin = ggplot2::margin(1, 1, 1, 1)
    )

  if (fit == "stim") {
    plot <- plot +
      ggplot2::scale_linetype_manual(
        values = .acsCytofScatterAltLinetypes(), labels = stimLabels,
        limits = stimLevels, name = stimName
      ) +
      # Longer keys so the stimulation line types can be told apart.
      ggplot2::theme(legend.key.width = ggplot2::unit(9, "mm"))
  }
  if (byMethod) {
    keys <- .acsCytofScatterAltMethods()[methods]
    colours <- stats::setNames(unname(.analysis_method_colours[keys]), methods)
    methodLabels <- .acsCytofScatterAltMethodLabels(methods)
    plot +
      ggplot2::scale_colour_manual(
        values = colours, labels = methodLabels, limits = methods, name = "Method"
      ) +
      ggplot2::scale_fill_manual(
        values = ggplot2::alpha(colours, fillAlpha), labels = methodLabels,
        limits = methods, name = "Method"
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(order = 1, override.aes = list(shape = 21, size = 2.2)),
        fill = ggplot2::guide_legend(order = 1),
        shape = ggplot2::guide_legend(
          order = 2, nrow = 2,
          override.aes = list(colour = "grey20", fill = ggplot2::alpha("grey20", fillAlpha), size = 2.2)
        ),
        linetype = ggplot2::guide_legend(order = 2, nrow = 2)
      )
  } else {
    colours <- .acsCytofValidationStimColours()
    plot +
      ggplot2::scale_colour_manual(
        values = colours, labels = stimLabels, limits = stimLevels, name = stimName
      ) +
      ggplot2::scale_fill_manual(
        values = ggplot2::alpha(colours, fillAlpha), labels = stimLabels,
        limits = stimLevels, name = stimName
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(nrow = 2, override.aes = list(size = 2.2)),
        fill = ggplot2::guide_legend(nrow = 2),
        shape = ggplot2::guide_legend(nrow = 2),
        linetype = ggplot2::guide_legend(nrow = 2)
      )
  }
}

# Rows of populations (labelled once, on the right) and columns of cytokines.
# Returns the composite plot, the row plots, the data, panel ranges and lines,
# and the printed size in mm.
.acsCytofScatterAlt <- function(
  comparisonTbl, method = NULL, fit = c("stim", "panel", "none")
) {
  fit <- match.arg(fit)
  methods <- names(.acsCytofScatterAltMethods())
  if (!is.null(method)) {
    if (length(method) != 1L || !method %in% methods) {
      stop("method must be NULL or one of: ", paste(methods, collapse = ", "), ".")
    }
    methods <- method
  }
  plotTbl <- .acsCytofScatterAltData(comparisonTbl, methods)
  if (!nrow(plotTbl)) stop("No matched ACS rows are available to plot.")
  byMethod <- is.null(method)
  limTbl <- .acsCytofScatterAltLimits(plotTbl)
  by <- c(if (byMethod) "method", "pop", "cyt", if (fit == "stim") "stim")
  fitTbl <- if (fit == "none") NULL else .acsCytofScatterAltFits(plotTbl, by, limits = limTbl)
  if (!is.null(fitTbl)) fitTbl$.line <- interaction(fitTbl[by], drop = TRUE)

  populations <- levels(droplevels(plotTbl$pop))
  rows <- stats::setNames(lapply(populations, function(p) {
    .acsCytofScatterAltRow(
      plotTbl[plotTbl$pop == p, , drop = FALSE],
      if (is.null(fitTbl)) NULL else fitTbl[fitTbl$pop == p, , drop = FALSE],
      limTbl[limTbl$pop == p, , drop = FALSE],
      byMethod, fit, methods
    )
  }), populations)

  textColour <- "#333333"
  rowLabel <- function(p) {
    cowplot::ggdraw() +
      cowplot::draw_label(p, angle = -90, size = 7.5, colour = textColour)
  }
  rowGrobs <- lapply(populations, function(p) {
    cowplot::plot_grid(
      rows[[p]] + ggplot2::theme(legend.position = "none"), rowLabel(p),
      nrow = 1, rel_widths = c(1, 0.035)
    )
  })
  panels <- cowplot::plot_grid(plotlist = rowGrobs, ncol = 1)
  yTitle <- cowplot::ggdraw() + cowplot::draw_label(
    "Background-subtracted frequency (automated gating, %)",
    angle = 90, size = 8.5, colour = textColour
  )
  xTitle <- cowplot::ggdraw() + cowplot::draw_label(
    "Background-subtracted frequency (manual gating, %)",
    size = 8.5, colour = textColour
  )
  legend <- cowplot::get_plot_component(rows[[1L]], "guide-box-bottom", return_all = TRUE)
  if (is.list(legend) && !inherits(legend, "grob")) legend <- legend[[1L]]
  rowHeight <- 27
  legendHeight <- if (byMethod) 16 else 11
  body <- cowplot::plot_grid(yTitle, panels, nrow = 1, rel_widths = c(0.035, 1))
  plot <- cowplot::plot_grid(
    body, xTitle, legend, ncol = 1,
    rel_heights = c(length(populations) * rowHeight, 5, legendHeight)
  ) + ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  list(
    plot = plot, rows = rows, data = plotTbl, limits = limTbl, fits = fitTbl,
    width = .analysis_thesis_width,
    height = length(populations) * rowHeight + 5 + legendHeight
  )
}
