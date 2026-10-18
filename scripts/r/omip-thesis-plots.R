# Saved-result OMIP thesis figures. Source analysis-plot-style.R,
# acs_cytof-plot_cyt.R and omip016-methods.R first.
.omip111ThesisRows <- function(results, classification) {
  keys <- c("strain", "mouse", "population", "marker", "method")
  if (!all(c(keys, "propBs", "manualBs", "referenceAccepted") %in% names(results)) ||
      !all(c(keys, "f1") %in% names(classification)) ||
      anyDuplicated(results[keys]) || anyDuplicated(classification[keys])) {
    stop("Missing fields or duplicate OMIP-111 thesis outcomes.")
  }
  if (nrow(dplyr::anti_join(results, classification, by = keys)) ||
      nrow(dplyr::anti_join(classification, results, by = keys))) {
    stop("OMIP-111 frequency and classification outcomes do not match.")
  }
  # Keep every available cell-count denominator for auditing F1, while leaving
  # the original saved frequency columns intact and avoiding join suffixes.
  auditFields <- c("f1", "tp", "fp", "fn", "tn", "n", "n_defined", "nCellStim",
                   "n_selected", "n_manual_positive", "n_manual_negative", "n_f1",
                   "precision", "sensitivity", "specificity", "fdp")
  auditFields <- intersect(auditFields, names(classification))
  overlapping <- intersect(auditFields, names(results))
  aligned <- dplyr::left_join(results[keys], classification[c(keys, overlapping)], by = keys)
  for (field in overlapping) {
    if (!isTRUE(all.equal(results[[field]], aligned[[field]], check.attributes = FALSE))) {
      stop("Saved frequency and classification audit fields disagree: ", field)
    }
  }
  rows <- dplyr::left_join(results,
    classification[c(keys, setdiff(auditFields, names(results)))], by = keys)
  rows$condition <- rows$strain
  rows$reference <- 100 * rows$manualBs
  rows$estimate <- 100 * rows$propBs
  rows$scatterEligible <- rows$referenceAccepted
  rows
}

.omipThesisSummary <- function(rows) {
  required <- c("method", "population", "marker", "condition", "reference",
                "estimate", "f1", "scatterEligible")
  if (!all(required %in% names(rows)) || !nrow(rows) ||
      anyNA(rows[c("method", "population", "marker", "condition")]) ||
      !is.logical(rows$scatterEligible) ||
      any(is.finite(rows$f1) & (rows$f1 < 0 | rows$f1 > 1))) {
    stop("Invalid normalized OMIP thesis rows.")
  }
  groups <- c("condition", "population", "marker")
  correlations <- .omipFrequencyCorrelationTable(
    rows, groups, "estimate", "reference", "population"
  )
  f1 <- rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(groups, "method")))) |>
    dplyr::summarise(
      n_f1_total = dplyr::n(), n_f1_finite = sum(is.finite(.data$f1)),
      n_f1_missing = sum(!is.finite(.data$f1)),
      n_scatter_eligible = sum(.data$scatterEligible %in% TRUE),
      f1 = if (any(is.finite(.data$f1))) mean(.data$f1[is.finite(.data$f1)]) else NA_real_,
      .groups = "drop"
    )
  list(strata = dplyr::left_join(
    correlations[correlations$scope == "stratum", ], f1,
    by = c(groups, "method")
  ), pooled = correlations[correlations$scope == "pooled", ])
}

# parseMethods: method labels are plotmath (e.g. "F[beta]") for the strips.
.omipThesisHeatmap <- function(summary, metric, palette, methodLabels, parseMethods = FALSE) {
  limits <- if (metric == "f1") c(0, 1) else c(-1, 1)
  methods <- names(methodLabels)
  rows <- summary
  rows$method <- factor(rows$method, levels = methods, labels = methodLabels)
  rows$marker <- factor(rows$marker, levels = unique(summary$marker))
  rows$population <- factor(rows$population, levels = rev(unique(summary$population)))
  rows$condition <- factor(rows$condition, levels = unique(summary$condition))
  colours <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  stops <- seq(limits[1], limits[2], length.out = 11)
  if (palette == "midpoint_0_5") {
    colours[6] <- "#FFFFFF"
    colours <- c(colours[1], colours)
    stops <- c(-1, seq(0, 1, length.out = 11))
  }
  fill <- ggplot2::scale_fill_gradientn(
    colours = colours, values = scales::rescale(stops, from = limits),
    limits = limits, na.value = "white",
    breaks = if (metric == "f1") seq(0, 1, 0.25) else seq(-1, 1, 0.5),
    name = switch(metric, pearson = "Pearson correlation", ccc = "Concordance correlation",
                  f1 = "Mean F1 score"),
    guide = ggplot2::guide_colourbar(title.position = "top",
      barwidth = ggplot2::unit(38, "mm"), barheight = ggplot2::unit(3, "mm"))
  )
  if (palette == "viridis" || metric == "f1") {
    fill$palette <- function(x) {
      .acsCytofF1Palette(pmax(scales::rescale(x, to = limits, from = c(0, 1)), 0))
    }
  }
  values <- rows[[metric]]
  mapped <- fill$palette(scales::rescale(values, from = limits))
  rgb <- grDevices::col2rgb(ifelse(is.na(mapped), "white", mapped)) / 255
  linear <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  luminance <- as.numeric(c(0.2126, 0.7152, 0.0722) %*% linear)
  rows$textColour <- ifelse(luminance < sqrt(0.0525) - 0.05, "white", "black")
  rows$label <- ifelse(is.finite(values), sprintf("%.2f", values), "NA")
  # Only observed combinations get tiles/text; unmeasured marker panels stay blank.
  separators <- data.frame(method = factor(methodLabels[-length(methods)], levels = methodLabels),
                           boundary = 0.5)
  ggplot2::ggplot(rows, ggplot2::aes(.data$marker, .data$population)) +
    .analysis_theme(grid = "none") +
    ggplot2::geom_tile(ggplot2::aes(fill = .data[[metric]]),
      colour = if (palette == "standard" && metric != "f1") NA else "grey80", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, colour = .data$textColour), size = 2.5) +
    ggplot2::scale_colour_identity() + fill +
    ggplot2::facet_grid(ggplot2::vars(method), ggplot2::vars(condition), drop = FALSE,
      labeller = if (parseMethods) ggplot2::labeller(method = ggplot2::label_parsed) else "label_value") +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::geom_hline(data = separators, ggplot2::aes(yintercept = .data$boundary),
      inherit.aes = FALSE, linewidth = 1.2, colour = "black") +
    .analysis_thesis_theme() +
    ggplot2::labs(x = "Cytokine", y = "Population") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      strip.text = ggplot2::element_text(size = 8), panel.grid = ggplot2::element_blank())
}

.omipThesisScatterPanel <- function(rows, marker, ranges, conditionColours, conditionShapes,
                                   referenceLabel) {
  d <- rows[rows$marker == marker & rows$scatterEligible %in% TRUE, ]
  # The marker is named by a facet strip (figures carry no titles); the blank
  # layer keeps the strip when the marker has no eligible rows.
  d$.marker_strip <- rep(marker, nrow(d))
  ggplot2::ggplot(d, ggplot2::aes(.data$reference, .data$estimate)) +
    ggplot2::geom_blank(data = data.frame(.marker_strip = marker), inherit.aes = FALSE) +
    .analysis_theme(grid = "xy") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$condition, shape = .data$condition),
      alpha = 0.5, size = 1.2, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = conditionColours, drop = FALSE) +
    ggplot2::scale_shape_manual(values = conditionShapes, drop = FALSE) +
    ggplot2::scale_x_continuous(limits = ranges[[marker]], n.breaks = 3,
      labels = .analysis_label_number, guide = ggplot2::guide_axis(check.overlap = TRUE)) +
    ggplot2::scale_y_continuous(limits = ranges[[marker]], n.breaks = 3,
      labels = .analysis_label_number) +
    .analysis_thesis_theme() +
    ggplot2::facet_wrap(ggplot2::vars(.marker_strip)) +
    ggplot2::labs(x = NULL, y = NULL, colour = NULL, shape = NULL) +
    ggplot2::theme(legend.position = "none", axis.title = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_text(size = 7), strip.text = ggplot2::element_text(size = 8))
}

# One row of cytokine panels per population: colour gives the method and
# shape the condition (stimulation or strain). Axes are square per cytokine.
.omipThesisScatterCombined <- function(rows, ranges, methodLabels, conditionShapes,
                                       conditionLabels) {
  d <- rows[rows$scatterEligible %in% TRUE, ]
  markers <- unique(d$marker)
  d$marker <- factor(d$marker, levels = markers)
  d$method <- factor(methodLabels[as.character(d$method)], levels = unname(methodLabels))
  d$condition <- factor(d$condition, levels = names(conditionShapes))
  keys <- match(unname(methodLabels), .analysis_method_labels)
  colours <- stats::setNames(ifelse(is.na(keys), "grey40",
    .analysis_method_colours[names(.analysis_method_labels)[keys]]), unname(methodLabels))
  limits <- do.call(rbind, lapply(markers, function(marker) data.frame(
    marker = factor(marker, levels = markers), reference = ranges[[marker]],
    estimate = ranges[[marker]])))
  ggplot2::ggplot(d, ggplot2::aes(.data$reference, .data$estimate)) +
    .analysis_theme(grid = "xy") +
    ggplot2::geom_blank(data = limits) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
    # Hollow, semi-transparent shapes keep methods that agree closely visible
    # on top of each other.
    ggplot2::geom_point(ggplot2::aes(colour = .data$method, shape = .data$condition),
      size = 1.8, stroke = 0.6, alpha = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(ggplot2::vars(.data$marker), nrow = 1, scales = "free") +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(1, 2, 0, 5, 3, 4)[seq_along(conditionShapes)],
      labels = conditionLabels, drop = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.04), n.breaks = 3,
      labels = .analysis_label_number, guide = ggplot2::guide_axis(check.overlap = TRUE)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.04), n.breaks = 3,
      labels = .analysis_label_number) +
    .analysis_thesis_theme() +
    ggplot2::labs(x = "Background-subtracted frequency (manual gating, %)",
      y = "Background-subtracted frequency\n(automated gating, %)",
      colour = NULL, shape = NULL) +
    ggplot2::theme(aspect.ratio = 1, axis.text = ggplot2::element_text(size = 7),
      axis.title = ggplot2::element_text(size = 8), strip.text = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8), legend.box = "horizontal",
      legend.box.spacing = ggplot2::unit(1, "mm"), legend.margin = ggplot2::margin(0, 8, 0, 8),
      plot.margin = ggplot2::margin(4, 4, 2, 4))
}

.omipThesisFigures <- function(rows, pathOutput, referenceLabel, methodLabels = NULL,
                               conditionLabels = NULL) {
  summaries <- .omipThesisSummary(rows)
  methods <- unique(as.character(rows$method))
  methods <- c(methods[grepl("[Ss]tim[Gg]ate", methods)],
    methods[grepl("[Ff].?beta", methods)], methods[grepl("[Tt]ailgate", methods)],
    methods[!grepl("[Ss]tim[Gg]ate|[Ff].?beta|[Tt]ailgate", methods)])
  if (is.null(methodLabels)) methodLabels <- stats::setNames(methods, methods)
  if (anyDuplicated(methods) || !setequal(names(methodLabels), methods) ||
      anyDuplicated(unname(methodLabels))) stop("Invalid OMIP thesis method labels.")
  methodLabels <- methodLabels[methods]
  markers <- unique(rows$marker)
  conditions <- unique(rows$condition)
  populations <- unique(rows$population)
  conditionColours <- stats::setNames(rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                                           "#E69F00", "#56B4E9", "#000000"),
                                         length.out = length(conditions)), conditions)
  conditionShapes <- stats::setNames(rep(c(16, 17, 15, 18, 3, 8),
                                        length.out = length(conditions)), conditions)
  if (is.null(conditionLabels)) conditionLabels <- stats::setNames(conditions, conditions)
  if (!all(conditions %in% names(conditionLabels))) stop("Invalid OMIP thesis condition labels.")
  ranges <- stats::setNames(lapply(markers, function(marker) {
    d <- rows[rows$marker == marker & rows$scatterEligible %in% TRUE, ]
    vals <- c(0, d$reference, d$estimate)
    r <- range(vals[is.finite(vals)])
    if (diff(r) == 0) r <- r + c(-0.5, 0.5)
    r + c(-1, 1) * diff(r) * 0.05
  }), markers)
  slug <- function(x) gsub("[^a-z0-9]+", "-", tolower(x))
  files <- character()
  .acsCytofReplaceDir(pathOutput, function(staging) {
    utils::write.csv(summaries$strata, file.path(staging, "stratum-scores-coverage.csv"), row.names = FALSE)
    utils::write.csv(summaries$pooled, file.path(staging, "pooled-correlations.csv"), row.names = FALSE)
    utils::write.csv(rows, file.path(staging, "individual-scores.csv"), row.names = FALSE)
    save <- function(p, name, height) {
      .analysis_save_thesis_fig(p, file.path(staging, name),
                               min(height, .analysis_fig_max_height * 10))
      files <<- c(files, paste0(name, c(".pdf", ".png")))
    }
    for (metric in c("pearson", "ccc", "f1")) {
      for (palette in if (metric == "f1") "viridis" else c("standard", "midpoint_0_5", "viridis")) {
        caption <- if (metric == "f1") {
          paste0("Reference: ", referenceLabel, ". Equal-weight mean stimulated-tube F1. No intervals.")
        } else {
          paste0("Reference: ", referenceLabel, ". Background-subtracted frequencies. No intervals.")
        }
        if (metric != "f1" && all(summaries$strata$n_total == 1L)) {
          caption <- paste(caption, "One donor per stratum: correlations undefined.")
        }
        caption <- paste(strwrap(caption, width = 100), collapse = "\n")
        plot <- .omipThesisHeatmap(summaries$strata, metric, palette, methodLabels) +
          ggplot2::labs(caption = caption) +
          ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6),
            plot.caption.position = "plot",
            plot.caption = ggplot2::element_text(size = 7, hjust = 0.5))
        save(plot,
          file.path("heatmaps", palette, paste0(metric, "-stacked")),
          45 + length(methods) * max(25, 12 * length(unique(rows$population))))
      }
    }
    panel <- function(d, marker) .omipThesisScatterPanel(
      d, marker, ranges, conditionColours, conditionShapes, referenceLabel)
    legendPlot <- panel(rows, markers[[1]]) + ggplot2::theme(legend.position = "bottom")
    legend <- cowplot::get_legend(legendPlot)
    stack <- function(blocks, labels, ncols, nrows) {
      bodies <- lapply(seq_along(blocks), function(i) cowplot::plot_grid(
        cowplot::ggdraw() + cowplot::draw_label(labels[[i]], size = 8, fontface = "bold"),
        cowplot::plot_grid(plotlist = blocks[[i]], ncol = ncols, nrow = nrows),
        ncol = 1, rel_heights = c(0.12, 1)))
      content <- cowplot::plot_grid(cowplot::plot_grid(plotlist = bodies, ncol = 1), legend,
                                    ncol = 1, rel_heights = c(1, 0.07))
      cowplot::ggdraw() +
        cowplot::draw_plot(content, x = 0.07, y = 0.06, width = 0.93, height = 0.94) +
        cowplot::draw_label("Background-subtracted frequency (manual gating, %)",
                           x = 0.54, y = 0.015, size = 8) +
        cowplot::draw_label("Background-subtracted frequency (automated gating, %)",
                           x = 0.015, y = 0.55, angle = 90, size = 8) +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
    }
    for (population in populations) {
      d <- rows[rows$population == population, ]
      ownMarkers <- unique(d$marker)
      save(.omipThesisScatterCombined(d, ranges, methodLabels, conditionShapes,
          conditionLabels[conditions]),
        file.path("scatter", paste0(slug(population), "-methods-one-row")),
        22 + .analysis_thesis_width / length(ownMarkers))
      for (layout in "compact") {
        ncols <- if (layout == "compact") min(3L, length(ownMarkers)) else length(ownMarkers)
        nrows <- ceiling(length(ownMarkers) / ncols)
        blocks <- lapply(methods, function(method) lapply(ownMarkers, function(marker) {
          panel(d[d$method == method, ], marker)
        }))
        save(stack(blocks, methodLabels, ncols, nrows),
          file.path("scatter", paste0(slug(population), "-methods-", layout)),
          15 + length(methods) * (12 + 38 * nrows))
      }
    }
    for (method in methods) {
      blocks <- lapply(populations, function(population) lapply(markers, function(marker) {
        d <- rows[rows$method == method & rows$population == population, ]
        if (!marker %in% d$marker) cowplot::ggdraw() else panel(d, marker)
      }))
      save(stack(blocks, paste(methodLabels[[method]], populations, sep = ": "), length(markers), 1),
        file.path("scatter", paste0(slug(method), "-populations-one-row")),
        15 + length(populations) * 50)
    }
    writeLines(c(referenceLabel,
      "Correlations use signed net frequencies in all existing report rows; finite pairs and intended coverage are recorded.",
      "F1 is the equal-weight mean of defined stimulated-tube F1 scores against supplied reference cell masks.",
      "Undefined observed scores are NA; unmeasured marker/population combinations are blank.",
      "Scatter figures retain the existing scatter eligibility cohort. Colours/shapes identify condition; marker-specific axes are common across methods/populations.",
      "Correlation legends retain [-1,1]; viridis colours for zero and negative correlations use the zero endpoint."),
      file.path(staging, "figure-notes.txt"))
  })
  invisible(list(files = file.path(pathOutput, files), summary = summaries))
}
