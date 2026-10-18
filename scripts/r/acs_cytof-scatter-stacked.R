# Saved, paired raw-frequency comparisons: retain the same strata as the
# existing scatterplots (the signal eligibility filter applies to correlations).
.acsCytofStackedScatterBlock <- function(comparisonTbl, method, population, ncol) {
  rows <- comparisonTbl |>
    dplyr::filter(.data$method == .env$method, .data$pop == .env$population)
  title <- paste(.acsCytofValidationMethodLabels()[[method]], population, sep = " — ")
  if (!nrow(rows)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::labs(title = title, subtitle = "No matched rows available"))
  }
  failed <- if ("thresholdFailed" %in% names(rows)) rows$thresholdFailed else rep(FALSE, nrow(rows))
  # Unknown failure status is missing, never a valid automated estimate.
  rows$freq_bs_auto[is.na(failed) | failed] <- NA_real_
  rows$freq_bs_man[!is.finite(rows$freq_bs_man)] <- NA_real_
  rows$freq_bs_auto[!is.finite(rows$freq_bs_auto)] <- NA_real_
  paired <- is.finite(rows$freq_bs_man) & is.finite(rows$freq_bs_auto)
  plot <- .acsCytofValidationPlotScatter(rows, method)
  for (i in seq_along(plot$layers)) {
    if (inherits(plot$layers[[i]]$geom, "GeomPoint")) plot$layers[[i]]$geom_params$na.rm <- TRUE
  }
  # Replace inherited scales directly rather than adding duplicate scales.
  plot$scales$scales <- Filter(function(s) !any(c("colour", "shape") %in% s$aesthetics), plot$scales$scales)
  plot +
    ggplot2::facet_wrap(ggplot2::vars(cyt), scales = "free", ncol = ncol, drop = FALSE) +
    ggplot2::scale_colour_manual(
      values = .acsCytofValidationStimColours(),
      labels = .acsCytofValidationStimLabels(), drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = c(mtb = 17, p1 = 15, p4 = 3, ebv = 16),
      labels = .acsCytofValidationStimLabels(), drop = FALSE
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        sum(paired), "/", nrow(rows), " matched pairs plotted; ",
        sum(!paired), " missing; ", sum(failed, na.rm = TRUE), " failed thresholds"
      )
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 9)
    )
}

# Select one population to stack methods, or one method to stack populations.
# Returns the composite and physical dimensions so six cytokines remain legible.
.acsCytofStackedScatter <- function(
  comparisonTbl, population = NULL, method = NULL, ncol = 6L
) {
  .acsCytofValidationValidateComparisonTable(comparisonTbl)
  if (is.null(population) == is.null(method)) {
    stop("Select exactly one population or one method to stack.")
  }
  if (!ncol %in% c(3L, 6L)) stop("Scatter columns must be 3 or 6.")
  methods <- c("stimgate", "fbeta_default", "tailgate")
  populations <- c("CD4 T cells", "CD8 T cells", "TCRgd T cells", "B cells", "NK cells")
  available <- unique(as.character(comparisonTbl$pop))
  populations <- c(intersect(populations, available), setdiff(available, populations))
  if (!is.null(population)) {
    population <- match.arg(population, populations)
    groups <- methods
    blocks <- lapply(methods, function(m) {
      .acsCytofStackedScatterBlock(comparisonTbl, m, population, ncol)
    })
  } else {
    if (length(method) != 1L || !method %in% methods) {
      stop("method must be one of: ", paste(methods, collapse = ", "), ".")
    }
    groups <- populations
    blocks <- lapply(populations, function(p) {
      .acsCytofStackedScatterBlock(comparisonTbl, method, p, ncol)
    })
  }
  legends <- lapply(blocks, cowplot::get_legend)
  legends <- Filter(Negate(is.null), legends)
  panels <- cowplot::plot_grid(
    plotlist = lapply(blocks, function(p) p + ggplot2::theme(legend.position = "none")),
    ncol = 1, align = "v", axis = "lr"
  )
  if (length(legends)) {
    panels <- cowplot::plot_grid(
      panels, legends[[1L]], ncol = 1,
      rel_heights = c(length(blocks) * if (ncol == 3L) 160 else 100, 15)
    )
  }
  list(
    plot = panels + ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    ), blocks = blocks, groups = groups,
    width = if (ncol == 3L) 230 else 420,
    height = length(blocks) * (if (ncol == 3L) 160 else 100) + 15
  )
}
