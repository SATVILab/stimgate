# Alternative ACS frequency scatter plots for the thesis (Analysis 10): every
# population in one figure, with populations as facet rows and cytokines as
# facet columns. Source after analysis-plot-style.R and acs_cytof-plot_cyt.R.
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
# Points are fillable shapes with an opaque outline and a translucent fill of
# the same colour: the outline keeps each point visible where points overlap,
# and the fill shows where they pile up near zero.

# Methods compared in the thesis, mapped to their shared colour/label keys.
.acsCytofScatterAltMethods <- function() {
  c(stimgate = "stimgate", fbeta_default = "fbeta", tailgate = "tailgate")
}

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
  rows |>
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
}

# Least-squares line of automated on manual frequency within each group,
# drawn over that group's manual range. Groups with fewer than two points or
# no spread in manual frequency have no line.
.acsCytofScatterAltFits <- function(plotTbl, by) {
  groups <- split(plotTbl, plotTbl[by], drop = TRUE)
  fits <- lapply(groups, function(g) {
    x <- g$freq_bs_man
    if (nrow(g) < 2L || diff(range(x)) <= 0) return(NULL)
    coef <- stats::coef(stats::lm(freq_bs_auto ~ freq_bs_man, data = g))
    xs <- range(x)
    out <- g[c(1L, 1L), by, drop = FALSE]
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
  out
}

# Plotmath method labels for legends (F-beta as F with subscript beta).
.acsCytofScatterAltMethodLabels <- function(methods) {
  keys <- .acsCytofScatterAltMethods()[methods]
  parse(text = ifelse(
    keys == "fbeta", "F[beta]",
    paste0("'", unname(.analysis_method_labels[keys]), "'")
  ))
}

# Populations as facet rows and cytokines as columns, each panel with its own
# x and y range (as in the stacked scatter). Returns the plot with its printed
# size in mm.
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
  stimLabels <- .acsCytofValidationStimLabels()
  stimLevels <- levels(plotTbl$stim)
  stimName <- "Stimulation"
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

  plot <- ggplot2::ggplot(plotTbl) +
    .analysis_theme(base_size = 8) +
    .analysis_thesis_theme() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey65", linewidth = 0.4)

  if (fit != "none") {
    by <- c(if (byMethod) "method", "pop", "cyt", if (fit == "stim") "stim")
    fitTbl <- .acsCytofScatterAltFits(plotTbl, by)
    fitTbl$.line <- interaction(fitTbl[by], drop = TRUE)
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
    ggplot2::geom_point(mapping = pointAes, size = 1.5, stroke = 0.45) +
    ggplot2::facet_wrap(
      ggplot2::vars(pop, cyt), ncol = 6L, scales = "free"
    ) +
    ggplot2::scale_x_continuous(
      labels = .analysis_label_number, n.breaks = 3,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    ) +
    ggplot2::scale_y_continuous(labels = .analysis_label_number, n.breaks = 3) +
    ggplot2::scale_shape_manual(
      values = .acsCytofScatterAltShapes(), labels = stimLabels,
      limits = stimLevels, name = stimName
    ) +
    ggplot2::labs(
      x = "Background-subtracted frequency (manual gating, %)",
      y = "Background-subtracted frequency (automated gating, %)"
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 6.5, margin = ggplot2::margin(1.5, 1, 1.5, 1)),
      axis.text = ggplot2::element_text(size = 6),
      panel.spacing = ggplot2::unit(1.5, "mm"),
      legend.box = "vertical",
      legend.spacing.y = ggplot2::unit(0, "mm")
    )

  if (fit == "stim") {
    plot <- plot + ggplot2::scale_linetype_manual(
      values = .acsCytofScatterAltLinetypes(), labels = stimLabels,
      limits = stimLevels, name = stimName
    )
  }
  fillAlpha <- 0.35
  if (byMethod) {
    keys <- .acsCytofScatterAltMethods()[methods]
    colours <- stats::setNames(unname(.analysis_method_colours[keys]), methods)
    methodLabels <- .acsCytofScatterAltMethodLabels(methods)
    plot <- plot +
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
    plot <- plot +
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
  if (fit == "stim") {
    # Longer keys so the stimulation line types can be told apart.
    plot <- plot + ggplot2::theme(legend.key.width = ggplot2::unit(9, "mm"))
  }
  nRow <- ceiling(length(unique(interaction(plotTbl$pop, plotTbl$cyt, drop = TRUE))) / 6)
  list(
    plot = plot, data = plotTbl,
    width = .analysis_thesis_width,
    height = nRow * 30 + if (byMethod) 28 else 22
  )
}
