.acsCytofValidationStimColours <- function() {
  c(
    mtb = "#e34234",
    p1 = "#40826d",
    p4 = "#ff00ff",
    ebv = "#7700cf"
  )
}

.acsCytofValidationStimLabels <- function() {
  c(
    mtb = "Live Mtb",
    p1 = "Secreted Mtb proteins",
    ebv = "EBV and CMV",
    p4 = "Non-secreted Mtb proteins"
  )
}

.acsCytofValidationMethodLabels <- function() {
  c(
    stimgate = "StimGate",
    fbeta = "F-beta",
    tailgate = "Tailgate"
  )
}

.acsCytofValidationRealPopulations <- function() {
  c("CD4 T cells", "CD8 T cells", "TCRgd T cells")
}

.acsCytofValidationCcc <- function(x, y) {
  complete <- stats::complete.cases(x, y)
  x <- x[complete]
  y <- y[complete]
  if (length(x) < 2L || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }

  covariance <- stats::cov(x, y)
  (2 * covariance) /
    (stats::var(x) + stats::var(y) + (mean(x) - mean(y))^2)
}

.acsCytofValidationCorrelationTable <- function(comparisonTbl) {
  comparisonTbl |>
    dplyr::group_by(method, pop, cyt, stim) |>
    dplyr::summarise(
      n = sum(stats::complete.cases(.data$freq_bs_auto, .data$freq_bs_man)),
      pcc = if (n > 1L) {
        suppressWarnings(stats::cor(
          .data$freq_bs_auto,
          .data$freq_bs_man,
          use = "complete.obs"
        ))
      } else {
        NA_real_
      },
      ccc = .acsCytofValidationCcc(
        .data$freq_bs_auto,
        .data$freq_bs_man
      ),
      .groups = "drop"
    )
}

.acsCytofValidationPlotScatter <- function(comparisonTbl, method) {
  method <- match.arg(method, c("stimgate", "fbeta", "tailgate"))
  plotTbl <- comparisonTbl |>
    dplyr::filter(.data$method == .env$method) |>
    dplyr::mutate(
      cyt = factor(
        .data$cyt,
        levels = c("IFNg", "IL2", "TNF", "IL17", "IL22", "IL6")
      ),
      pop = factor(
        .data$pop,
        levels = c(
          "CD4 T cells",
          "CD8 T cells",
          "TCRgd T cells",
          "B cells",
          "NK cells"
        )
      )
    )

  ggplot2::ggplot(plotTbl) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white"),
      panel.background = ggplot2::element_rect(fill = "white")
    ) +
    cowplot::background_grid(major = "xy") +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$freq_bs_man,
        y = .data$freq_bs_auto,
        colour = .data$stim
      )
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(pop, cyt),
      scales = "free",
      ncol = dplyr::n_distinct(plotTbl$cyt)
    ) +
    ggplot2::scale_colour_manual(
      values = .acsCytofValidationStimColours(),
      labels = .acsCytofValidationStimLabels()
    ) +
    ggplot2::labs(
      title = unname(.acsCytofValidationMethodLabels()[[method]]),
      x = "Background-subtracted frequency\n(manual gating)",
      y = "Background-subtracted frequency\n(automated gating)",
      colour = NULL
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
}

.acsCytofValidationPlotCorrelation <- function(
  correlationTbl,
  method,
  metric = c("pcc", "ccc"),
  realPopulationsOnly = TRUE
) {
  metric <- match.arg(metric)
  method <- match.arg(method, c("stimgate", "fbeta", "tailgate"))
  plotTbl <- correlationTbl |>
    dplyr::filter(
      .data$method == .env$method,
      .data$stim != "p4"
    )
  if (isTRUE(realPopulationsOnly)) {
    plotTbl <- plotTbl |>
      dplyr::filter(.data$pop %in% .acsCytofValidationRealPopulations())
  }

  plotTbl <- plotTbl |>
    dplyr::mutate(
      cyt = factor(
        .data$cyt,
        levels = c("IFNg", "IL2", "TNF", "IL17", "IL22", "IL6")
      ),
      pop = factor(
        .data$pop,
        levels = c(
          "CD4 T cells",
          "CD8 T cells",
          "TCRgd T cells",
          "B cells",
          "NK cells"
        )
      ),
      stim = factor(
        .data$stim,
        levels = c("p1", "mtb", "ebv"),
        labels = .acsCytofValidationStimLabels()[c("p1", "mtb", "ebv")]
      )
    )

  metricLabel <- if (metric == "pcc") {
    "Pearson correlation"
  } else {
    "Concordance correlation"
  }
  colourValues <- rev(RColorBrewer::brewer.pal(11, "RdBu"))

  ggplot2::ggplot(
    plotTbl,
    ggplot2::aes(x = .data$cyt, y = .data$pop)
  ) +
    cowplot::theme_cowplot() +
    ggplot2::geom_raster(ggplot2::aes(fill = .data[[metric]])) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(.data[[metric]], 2)),
      size = 2.25
    ) +
    ggplot2::facet_wrap(ggplot2::vars(stim), ncol = 3, scales = "fixed") +
    ggplot2::scale_fill_gradientn(
      colours = colourValues,
      values = seq(0, 1, length.out = length(colourValues)),
      limits = c(-1, 1),
      na.value = "gray75",
      name = metricLabel
    ) +
    ggplot2::labs(
      title = unname(.acsCytofValidationMethodLabels()[[method]]),
      x = "Cytokine",
      y = "Population"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(size = 9.5),
      legend.title = ggplot2::element_text(size = 10)
    )
}

.acsCytofValidationSavePlots <- function(comparisonTbl, pathDirSave) {
  dir.create(pathDirSave, recursive = TRUE, showWarnings = FALSE)
  correlationTbl <- .acsCytofValidationCorrelationTable(comparisonTbl)
  utils::write.csv(
    correlationTbl,
    file.path(pathDirSave, "manual-comparison-correlations.csv"),
    row.names = FALSE
  )
  saveRDS(
    correlationTbl,
    file.path(pathDirSave, "manual-comparison-correlations.rds")
  )

  methods <- intersect(
    c("stimgate", "fbeta", "tailgate"),
    as.character(unique(comparisonTbl$method))
  )
  for (method in methods) {
    scatter <- .acsCytofValidationPlotScatter(comparisonTbl, method)
    ggplot2::ggsave(
      file.path(pathDirSave, paste0("scatter-", method, ".pdf")),
      plot = scatter,
      height = 25,
      width = 40,
      units = "cm"
    )

    for (realOnly in c(TRUE, FALSE)) {
      populationSuffix <- if (realOnly) "real-pops" else "all-pops"
      for (metric in c("pcc", "ccc")) {
        correlationPlot <- .acsCytofValidationPlotCorrelation(
          correlationTbl = correlationTbl,
          method = method,
          metric = metric,
          realPopulationsOnly = realOnly
        )
        ggplot2::ggsave(
          file.path(
            pathDirSave,
            paste0("heatmap-", metric, "-", method, "-", populationSuffix, ".pdf")
          ),
          plot = correlationPlot,
          height = 12.5,
          width = 30,
          units = "cm"
        )
      }
    }
  }

  invisible(correlationTbl)
}
