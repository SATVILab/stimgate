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
    tailgate = "Tailgate",
    fbeta_default = "F-beta (default settings)",
    tailgate_default = "Tailgate (default settings)"
  )
}

# Method sets for every method-comparison figure: all methods, and without
# Tailgate. Names are the figure subfolders; `label` is the heading.
.acsCytofMethodSets <- function() {
  list(
    all_methods = list(
      label = "All methods",
      methods = c("stimgate", "tailgate", "fbeta")
    ),
    no_tailgate = list(
      label = "Without Tailgate",
      methods = c("stimgate", "fbeta")
    )
  )
}

.acsCytofValidationRealPopulations <- function() {
  c("CD4 T cells", "CD8 T cells", "TCRgd T cells")
}

.acsCytofValidationMethods <- function() {
  c("stimgate", "fbeta", "tailgate")
}

# Tailgate and F-beta at their published defaults (all cells, automatic
# Tailgate tolerance, no bias), reported in the appendix of Analysis 10.
.acsCytofValidationAppendixMethods <- function() {
  c("fbeta_default", "tailgate_default")
}

# The main results use StimGate and the tuned comparators only.
.acsCytofMainComparisonTable <- function(comparisonTbl) {
  out <- comparisonTbl |>
    dplyr::filter(.data$method %in% .acsCytofValidationMethods())
  if (is.factor(out$method)) out$method <- droplevels(out$method)
  for (nm in c("manifest", "exclusions", "summary")) {
    attr(out, nm) <- attr(comparisonTbl, nm)
  }
  out
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

  if (!"thresholdFailed" %in% names(comparisonTbl)) comparisonTbl$thresholdFailed <- FALSE
  correlationTbl <- comparisonTbl |>
    dplyr::mutate(freq_bs_auto = dplyr::if_else(.data$thresholdFailed, NA_real_, .data$freq_bs_auto)) |>
    dplyr::group_by(method, pop, cyt, stim) |>
    dplyr::filter(
      stats::quantile(.data$freq_stim_man, 0.75, na.rm = TRUE) >
        3 * max(0.01, stats::median(.data$freq_uns_man, na.rm = TRUE)),
      stats::quantile(.data$freq_bs_man, 0.75, na.rm = TRUE) > 0.02
    ) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_failed = sum(.data$thresholdFailed),
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
  keyCols <- c("method", "pop", "cyt", "stim")
  attr(correlationTbl, "excludedStrata") <- comparisonTbl |>
    dplyr::distinct(dplyr::across(dplyr::all_of(keyCols))) |>
    dplyr::anti_join(correlationTbl, by = keyCols)
  correlationTbl
}

.acsCytofValidationPlotScatter <- function(comparisonTbl, method) {
  .acsCytofValidationValidateComparisonTable(comparisonTbl)
  method <- match.arg(
    method,
    c(.acsCytofValidationMethods(), .acsCytofValidationAppendixMethods())
  )
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
    .analysis_theme() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 7)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$freq_bs_man,
        y = .data$freq_bs_auto,
        colour = .data$stim,
        shape = .data$stim
      ),
      alpha = 0.5,
      size = 0.9
    ) +
    ggplot2::scale_x_continuous(
      labels = .analysis_label_number,
      n.breaks = 3,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    ) +
    ggplot2::scale_y_continuous(labels = .analysis_label_number, n.breaks = 3) +
    ggplot2::facet_wrap(
      ggplot2::vars(cyt),
      scales = "free",
      ncol = 3
    ) +
    ggplot2::scale_colour_manual(
      values = .acsCytofValidationStimColours(),
      labels = .acsCytofValidationStimLabels()
    ) +
    ggplot2::scale_shape_manual(
      values = c(mtb = 17, p1 = 15, p4 = 3, ebv = 16),
      labels = .acsCytofValidationStimLabels()
    ) +
    ggplot2::labs(
      x = "Background-subtracted frequency\n(manual gating, %)",
      y = "Background-subtracted frequency\n(automated gating, %)",
      colour = NULL,
      shape = NULL
    )
}

.acsCytofValidationPlotCorrelation <- function(
  correlationTbl,
  method,
  metric = c("pcc", "ccc"),
  realPopulationsOnly = TRUE
) {
  metric <- match.arg(metric)
  method <- match.arg(
    method,
    c(.acsCytofValidationMethods(), .acsCytofValidationAppendixMethods())
  )
  requiredCols <- c("method", "pop", "cyt", "stim", metric)
  missingCols <- setdiff(requiredCols, names(correlationTbl))
  if (length(missingCols) > 0L) {
    stop(
      "ACS validation correlation table is missing required column(s): ",
      paste(missingCols, collapse = ", "),
      "."
    )
  }
  plotTbl <- correlationTbl |> dplyr::mutate(excluded = FALSE)
  excludedStrata <- attr(correlationTbl, "excludedStrata")
  if (!is.null(excludedStrata)) {
    plotTbl <- dplyr::bind_rows(
      plotTbl, dplyr::mutate(excludedStrata, excluded = TRUE)
    )
  }
  plotTbl <- plotTbl |>
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

  plotTbl$cellLabel <- ifelse(
    plotTbl$excluded, "excl.",
    ifelse(
      is.finite(plotTbl[[metric]]),
      format(round(plotTbl[[metric]], 2), trim = TRUE), "NA"
    )
  )
  plotTbl$textColour <- ifelse(
    !plotTbl$excluded & is.finite(plotTbl[[metric]]) & abs(plotTbl[[metric]]) >= 0.6,
    "white", "black"
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
    .analysis_theme(grid = "none") +
    ggplot2::geom_tile(ggplot2::aes(fill = .data[[metric]])) +
    ggplot2::geom_tile(
      data = dplyr::filter(plotTbl, .data$excluded), fill = "grey85"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$cellLabel, colour = .data$textColour),
      size = 3
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::facet_wrap(ggplot2::vars(stim), ncol = 3, scales = "fixed") +
    ggplot2::scale_fill_gradientn(
      colours = colourValues,
      values = seq(0, 1, length.out = length(colourValues)),
      limits = c(-1, 1),
      na.value = "white",
      name = metricLabel
    ) +
    ggplot2::labs(
      x = "Cytokine",
      y = "Population",
      caption = paste0(
        "Grey / excl.: excluded by signal eligibility rules.\n",
        "White / NA: eligible but correlation unavailable."
      )
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
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
    c(.acsCytofValidationMethods(), .acsCytofValidationAppendixMethods()),
    as.character(unique(comparisonTbl$method))
  )
  if (length(methods) == 0L) {
    stop("No supported ACS validation methods are available to plot.")
  }

  # Every validation figure shows one method, so each is saved once.
  pathDirSaveHeatmap <- file.path(stagedDir, "heatmaps")
  pathDirSaveScatter <- file.path(stagedDir, "scatter-plots")
  dir.create(pathDirSaveHeatmap, recursive = TRUE, showWarnings = FALSE)
  dir.create(pathDirSaveScatter, recursive = TRUE, showWarnings = FALSE)

  for (method in methods) {
    methodRows <- comparisonTbl |> dplyr::filter(.data$method == .env$method)
    for (population in unique(as.character(methodRows$pop))) {
      rows <- methodRows |> dplyr::filter(.data$pop == .env$population)
      populationDir <- file.path(
        pathDirSaveScatter, gsub("[^A-Za-z0-9_.-]+", "_", population)
      )
      dir.create(populationDir, recursive = TRUE, showWarnings = FALSE)
      scatter <- .acsCytofValidationPlotScatter(rows, method)
      .analysis_save_fig(
        scatter,
        file.path(populationDir, paste0(method, ".pdf")),
        height = 12
      )
    }

    for (realOnly in c(TRUE, FALSE)) {
      populationSuffix <- if (realOnly) "real-pops" else "all-pops"
      for (metric in c("pcc", "ccc")) {
        correlationPlot <- .acsCytofValidationPlotCorrelation(
          correlationTbl = correlationTbl,
          method = method,
          metric = metric,
          realPopulationsOnly = realOnly
        )
        .analysis_save_fig(
          correlationPlot,
          file.path(
            pathDirSaveHeatmap,
            paste0(metric, "-", method, "-", populationSuffix, ".pdf")
          ),
          height = if (realOnly) 9 else 12
        )
      }
    }
  }

  .acsCytofValidationReplaceDirectory(stagedDir, pathDirSave)

  invisible(correlationTbl)
}
