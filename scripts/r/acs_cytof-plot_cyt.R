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

.acsCytofValidationMethods <- function() {
  c("stimgate", "fbeta", "tailgate")
}

.acsCytofValidationValidateComparisonTable <- function(
  comparisonTbl,
  requiredMethods = NULL
) {
  if (!is.data.frame(comparisonTbl) || nrow(comparisonTbl) == 0L) {
    stop("ACS validation comparison table must be a non-empty data frame.")
  }

  requiredCols <- c(
    "method",
    "pop",
    "cyt",
    "stim",
    "SampleID",
    "freq_bs_auto",
    "freq_bs_man",
    "freq_stim_man",
    "freq_uns_man"
  )
  missingCols <- setdiff(requiredCols, names(comparisonTbl))
  if (length(missingCols) > 0L) {
    stop(
      "ACS validation comparison table is missing required column(s): ",
      paste(missingCols, collapse = ", "),
      "."
    )
  }

  frequencyCols <- c(
    "freq_bs_auto",
    "freq_bs_man",
    "freq_stim_man",
    "freq_uns_man"
  )
  nonNumeric <- frequencyCols[
    !vapply(comparisonTbl[frequencyCols], is.numeric, logical(1))
  ]
  if (length(nonNumeric) > 0L) {
    stop(
      "ACS validation frequency column(s) must be numeric: ",
      paste(nonNumeric, collapse = ", "),
      "."
    )
  }

  keyCols <- c("method", "pop", "cyt", "stim", "SampleID")
  if (any(vapply(comparisonTbl[keyCols], anyNA, logical(1)))) {
    stop("ACS validation comparison keys must not contain missing values.")
  }

  duplicateKeys <- comparisonTbl |>
    dplyr::count(dplyr::across(dplyr::all_of(keyCols)), name = "n") |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicateKeys) > 0L) {
    stop(
      "ACS validation comparison table contains duplicate ",
      "method/population/cytokine/stimulation/sample rows."
    )
  }

  if (!is.null(requiredMethods)) {
    missingMethods <- setdiff(
      as.character(requiredMethods),
      unique(as.character(comparisonTbl$method))
    )
    if (length(missingMethods) > 0L) {
      stop(
        "ACS validation comparison table is missing required method(s): ",
        paste(missingMethods, collapse = ", "),
        "."
      )
    }
  }

  invisible(TRUE)
}

# For paired, non-repeated two-method data this is algebraically equivalent
# to the U-statistic CCC used by cccrm::cccUst() in the ACS manuscript code.
.acsCytofValidationCcc <- function(x, y) {
  finite <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[finite])
  y <- as.numeric(y[finite])
  n <- length(x)
  if (n < 2L) {
    return(NA_real_)
  }

  meanX <- mean(x)
  meanY <- mean(y)
  varianceX <- mean((x - meanX)^2)
  varianceY <- mean((y - meanY)^2)
  covariance <- mean((x - meanX) * (y - meanY))
  denominator <- varianceX + varianceY + (meanX - meanY)^2

  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }

  (2 * covariance) / denominator
}

.acsCytofValidationCorrelationTable <- function(comparisonTbl) {
  .acsCytofValidationValidateComparisonTable(comparisonTbl)

  comparisonTbl |>
    dplyr::group_by(method, pop, cyt, stim) |>
    dplyr::filter(
      stats::quantile(.data$freq_stim_man, 0.75, na.rm = TRUE) >
        3 * max(0.01, stats::median(.data$freq_uns_man, na.rm = TRUE)),
      stats::quantile(.data$freq_bs_man, 0.75, na.rm = TRUE) > 0.02
    ) |>
    dplyr::summarise(
      n = sum(is.finite(.data$freq_bs_auto) & is.finite(.data$freq_bs_man)),
      pcc = {
        finite <- is.finite(.data$freq_bs_auto) &
          is.finite(.data$freq_bs_man)
        if (sum(finite) > 1L) {
          suppressWarnings(stats::cor(
            .data$freq_bs_auto[finite],
            .data$freq_bs_man[finite]
          ))
        } else {
          NA_real_
        }
      },
      ccc = .acsCytofValidationCcc(
        .data$freq_bs_auto,
        .data$freq_bs_man
      ),
      .groups = "drop"
    )
}

.acsCytofValidationPlotScatter <- function(comparisonTbl, method) {
  .acsCytofValidationValidateComparisonTable(comparisonTbl)
  method <- match.arg(method, .acsCytofValidationMethods())
  if (!method %in% as.character(unique(comparisonTbl$method))) {
    stop("No ACS validation rows are available for method: ", method, ".")
  }
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
  method <- match.arg(method, .acsCytofValidationMethods())
  requiredCols <- c("method", "pop", "cyt", "stim", metric)
  missingCols <- setdiff(requiredCols, names(correlationTbl))
  if (length(missingCols) > 0L) {
    stop(
      "ACS validation correlation table is missing required column(s): ",
      paste(missingCols, collapse = ", "),
      "."
    )
  }
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
      strip.background = ggplot2::element_rect(
        fill = "white",
        colour = "black"
      ),
      strip.text = ggplot2::element_text(size = 9.5),
      legend.title = ggplot2::element_text(size = 10)
    )
}

.acsCytofValidationReplaceDirectory <- function(stagedDir, targetDir) {
  parentDir <- dirname(targetDir)
  dir.create(parentDir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(targetDir) && !dir.exists(targetDir)) {
    stop("ACS validation output path exists but is not a directory: ", targetDir)
  }

  backupDir <- tempfile(
    pattern = paste0(".", basename(targetDir), "-previous-"),
    tmpdir = parentDir
  )
  hadTarget <- dir.exists(targetDir)

  if (hadTarget && !file.rename(targetDir, backupDir)) {
    stop("Could not move the previous ACS validation output directory aside.")
  }

  if (!file.rename(stagedDir, targetDir)) {
    restored <- TRUE
    if (hadTarget && dir.exists(backupDir)) {
      restored <- file.rename(backupDir, targetDir)
    }
    stop(
      "Could not promote the staged ACS validation output directory.",
      if (!isTRUE(restored)) {
        " Restoring the previous output directory also failed."
      } else {
        ""
      }
    )
  }

  if (dir.exists(backupDir)) {
    unlink(backupDir, recursive = TRUE, force = TRUE)
  }

  invisible(TRUE)
}

.acsCytofValidationSavePlots <- function(comparisonTbl, pathDirSave) {
  .acsCytofValidationValidateComparisonTable(comparisonTbl)

  parentDir <- dirname(pathDirSave)
  dir.create(parentDir, recursive = TRUE, showWarnings = FALSE)
  stagedDir <- tempfile(
    pattern = paste0(".", basename(pathDirSave), "-next-"),
    tmpdir = parentDir
  )
  dir.create(stagedDir, recursive = TRUE, showWarnings = FALSE)
  on.exit(
    if (dir.exists(stagedDir)) {
      unlink(stagedDir, recursive = TRUE, force = TRUE)
    },
    add = TRUE
  )

  pathDirSaveHeatmap <- file.path(stagedDir, "heatmaps")
  pathDirSaveScatter <- file.path(stagedDir, "scatter-plots")
  dir.create(pathDirSaveHeatmap, recursive = TRUE, showWarnings = FALSE)
  dir.create(pathDirSaveScatter, recursive = TRUE, showWarnings = FALSE)

  correlationTbl <- .acsCytofValidationCorrelationTable(comparisonTbl)
  utils::write.csv(
    correlationTbl,
    file.path(stagedDir, "manual-comparison-correlations.csv"),
    row.names = FALSE
  )
  .write_rds_atomic(
    correlationTbl,
    file.path(stagedDir, "manual-comparison-correlations.rds")
  )

  methods <- intersect(
    .acsCytofValidationMethods(),
    as.character(unique(comparisonTbl$method))
  )
  if (length(methods) == 0L) {
    stop("No supported ACS validation methods are available to plot.")
  }

  for (method in methods) {
    scatter <- .acsCytofValidationPlotScatter(comparisonTbl, method)
    ggplot2::ggsave(
      file.path(pathDirSaveScatter, paste0(method, ".pdf")),
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
            pathDirSaveHeatmap,
            paste0(
              metric,
              "-",
              method,
              "-",
              populationSuffix,
              ".pdf"
            )
          ),
          plot = correlationPlot,
          height = 12.5,
          width = 30,
          units = "cm"
        )
      }
    }
  }

  .acsCytofValidationReplaceDirectory(stagedDir, pathDirSave)

  invisible(correlationTbl)
}
