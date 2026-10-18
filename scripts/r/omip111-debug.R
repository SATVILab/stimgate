# Per-mouse diagnostics using the same recorder and panels as Analysis 13.
.omip111DebugLines <- function(rows, references, cofactor) {
  sg <- rows[rows$method == "StimGate", ]
  cyt <- sg$gateCyt
  if (isTRUE(all.equal(cyt, sg$threshold))) cyt <- NA_real_
  tibble::tibble(
    line = c("final gate", "conditional gate", "F-beta", "Tailgate", "author stim", "author uns"),
    label = c(
      "final StimGate gate", "StimGate: another cytokine positive",
      "F-beta gate", "Tailgate gate", "author cytokine cutoff: stimulated",
      "author cytokine cutoff: control"
    ),
    x = c(
      sg$threshold, cyt, rows$threshold[rows$method == "F-beta"],
      rows$threshold[rows$method == "Tailgate"],
      asinh(references$cytokineLowerRaw[match(c(sg$sampleStim, sg$sampleUns), references$sample)] / cofactor)
    ),
    colour = c(
      .analysis_method_colours[["stimgate"]], .analysis_method_colours[["stimgate"]],
      .analysis_method_colours[["fbeta"]],
      .analysis_method_colours[["tailgate"]], "#7B3294", "#7B3294"
    ),
    linetype = c("solid", "dashed", "solid", "solid", "longdash", "dotted"),
    linewidth = c(1.2, 0.8, 0.9, 0.9, 0.9, 0.9)
  ) |>
    dplyr::filter(is.finite(.data$x))
}

.omip111DebugFigure <- function(rec, rows, references, cofactor,
                                fitLabel = "CD44+ parent of the cytokine",
                                referenceLabel = "author gate (CD44+ parent)") {
  sg <- rows[rows$method == "StimGate", ]
  lines <- .omip111DebugLines(rows, references, cofactor)
  info <- .simDebugLocInfo(rec)
  info$estimate <- info$estimate |>
    dplyr::filter(!.data$name %in% c("response frequency (true)", "relative error", "false discovery", "sensitivity"))
  info$simulation <- tibble::tibble(
    name = c("strain / mouse", "population / marker", "fitted on", "stimulated tube", "control tube", "expression scale"),
    value = c(
      paste(sg$strain, sg$mouse), paste(sg$population, sg$marker), fitLabel,
      sg$sampleStim, sg$sampleUns, paste0("asinh(raw / ", cofactor, ")")
    )
  )
  info$final <- tibble::tibble(
    name = c(
      "saved gates reproduced", "StimGate source", "StimGate reason",
      paste0(referenceLabel, ": stim / control / net (%)"),
      paste0(rows$method, ": stim / control / net (%)"),
      paste0(rows$method, ": error (percentage points)")
    ),
    value = c(
      "yes (base and conditional gates)", sg$locSource, sg$locReason,
      sprintf("%.3f / %.3f / %.3f", 100 * sg$manualStim, 100 * sg$manualUns, 100 * sg$manualBs),
      sprintf("%.3f / %.3f / %.3f", 100 * rows$propStim, 100 * rows$propUns, 100 * rows$propBs),
      sprintf("%+.3f", rows$errorPp)
    )
  )
  headings <- c(
    simulation = "Sample", gating = "Gating settings",
    estimate = "Initial local-FDR gate", final = "Saved final method results"
  )
  if (!is.data.frame(rec$densTblRaw) && !is.data.frame(rec$dataMod)) {
    return(cowplot::plot_grid(plotlist = lapply(names(headings), function(nm) {
      .simDebugLocInfoPanel(info[[nm]], headings[[nm]], nRow = 30L)
    }), nrow = 1L))
  }
  .simDebugLocPlotGrid(.simDebugLocPlots(rec,
    xlim = .acsCytofDebugXlim(rec, lines), extraLines = lines
  ), info, headings)
}

# Print every diagnostic page in a collapsed block under strain, population and
# marker headings. An optional `note` column extends each block's summary.
.omip111DiagnosticPages <- function(index, pathDebug) {
  if (is.null(index$fit)) index$fit <- "cd44"
  fitTitles <- c(cd44 = "Fitted on CD44+ cells", nonnaive = "Fitted on all non-naïve cells")
  for (fit in intersect(names(fitTitles), unique(index$fit))) {
    if (length(unique(index$fit)) > 1L) cat("\n\n### ", fitTitles[[fit]], "\n\n", sep = "")
    .omip111DiagnosticFitPages(index[index$fit == fit, ], pathDebug,
      if (fit == "cd44") "CD44-high" else "non-naïve",
      if (length(unique(index$fit)) > 1L) "####" else "###")
  }
  invisible(NULL)
}

.omip111DiagnosticFitPages <- function(index, pathDebug, parentLabel, level) {
  sub <- paste0(level, "#")
  for (strain in c("C57", "BALB")) {
    for (population in c("CD4", "CD8")) {
      cat("\n\n", level, " ", strain, " ", parentLabel, " ", population, "\n\n", sep = "")
      subset <- index[index$strain == strain & index$population == population, ]
      for (marker in unique(subset$marker)) {
        cat("\n\n", sub, " ", marker, "\n\n", sep = "")
        pages <- subset[subset$marker == marker, ]
        for (i in seq_len(nrow(pages))) {
          note <- if (is.null(pages$note)) "" else pages$note[[i]]
          cat("\n<details><summary>", pages$mouse[[i]], " — ", marker,
            " (PDF page ", pages$page[[i]], ")", note, "</summary>\n\n",
            sep = ""
          )
          cat("![Gate diagnostics for ", pages$mouse[[i]], " ", population,
            " ", marker, "](", knitr::image_uri(file.path(pathDebug, pages$png[[i]])),
            ")\n",
            sep = ""
          )
          cat("\n\n</details>\n\n")
        }
      }
    }
  }
  invisible(NULL)
}

# The rerun must reproduce the saved gating. Gates may differ in the last bits
# when OpenBLAS uses a different thread count, so they must agree to `tol`,
# and every final combination count must be identical: the same cells are
# then positive under the base and conditional gates.
.omip111DebugCheckParity <- function(gates, saved, stats, savedStats, label, tol = 1e-10) {
  final <- function(x) x[x$gateName == "loc_minClust", , drop = FALSE]
  gates <- final(gates)
  saved <- final(saved)
  nSaved <- nrow(saved)
  key <- function(x) paste(x$ind, x$marker)
  saved <- saved[match(key(gates), key(saved)), , drop = FALSE]
  close <- function(a, b) {
    identical(is.na(a), is.na(b)) &&
      all(abs(a - b)[!is.na(a)] <= tol * pmax(1, abs(b[!is.na(b)])))
  }
  if (nrow(gates) == 0L || nrow(gates) != nSaved || anyNA(saved$ind) ||
    !close(gates$gate, saved$gate) || !close(gates$gateCyt, saved$gateCyt)) {
    stop("Diagnostic rerun differs from saved Analysis 14 gates: ", label)
  }
  countKey <- c("ind", "cytCombn")
  stats <- final(stats)[c(countKey, "countStim", "countUns")]
  savedStats <- final(savedStats)[c(countKey, "countStim", "countUns")]
  sortCounts <- function(x) x[do.call(order, unname(as.list(x[countKey]))), , drop = FALSE]
  stats <- sortCounts(stats)
  savedStats <- sortCounts(savedStats)
  rownames(stats) <- rownames(savedStats) <- NULL
  if (!isTRUE(all.equal(stats, savedStats, check.attributes = FALSE, tolerance = 0))) {
    stop("Diagnostic rerun changes saved Analysis 14 combination counts: ", label)
  }
  invisible(TRUE)
}

# Record and draw the pages of one StimGate rerun. `rowsFor(sid, marker)`
# returns the saved method rows of a stimulated tube; pages go under `fit`.
.omip111DebugPages <- function(pathTmp, records, fit, strain, population, markers, sampleMap,
                               rowsFor, prepared, settings, fitLabel, referenceLabel) {
  index <- list()
  pairs <- sampleMap[sampleMap$condition == "stim", ]
  for (marker in markers) {
    relDir <- file.path(fit, strain, population)
    dir.create(file.path(pathTmp, relDir), recursive = TRUE, showWarnings = FALSE)
    relPdf <- file.path(relDir, paste0(marker, ".pdf"))
    grDevices::pdf(file.path(pathTmp, relPdf), width = 36 / 2.54, height = 48 / 2.54)
    tryCatch(
      {
        for (i in seq_len(nrow(pairs))) {
          rows <- rowsFor(pairs$sample[[i]], marker)
          stopifnot(nrow(rows) == 3L)
          rec <- readRDS(file.path(records, paste0(pairs$ind[[i]], "_", marker, ".rds")))
          refs <- prepared$references[prepared$references$population == population &
            prepared$references$marker == marker, ]
          fig <- .omip111DebugFigure(rec, rows, refs, settings$preprocessing$cofactor,
            fitLabel = fitLabel, referenceLabel = referenceLabel)
          print(fig)
          relPng <- file.path(relDir, paste0(pairs$mouse[[i]], "-", marker, ".png"))
          ggplot2::ggsave(file.path(pathTmp, relPng), fig,
            width = 36 / 2.54, height = 48 / 2.54, dpi = 110, limitsize = FALSE,
            bg = "white"
          )
          index[[length(index) + 1L]] <- data.frame(fit, strain, population, marker,
            mouse = pairs$mouse[[i]], png = relPng, pdf = relPdf, page = i, matchesSavedGates = TRUE
          )
        }
      },
      finally = grDevices::dev.off()
    )
  }
  dplyr::bind_rows(index)
}

.omip111DebugGateStim <- function(settings, project, data, batches, markers, sampleMap, records) {
  dir.create(records, recursive = TRUE)
  .analysis_with_seed(settings$seed, .simDebugLoc(
    stimgate::gateStim(
      pathProject = project, .data = data, batchList = batches,
      marker = markers, biasUns = NULL, calcCytPosGates = isTRUE(settings$control$calcCytPosGates),
      control = do.call(
        stimgate::stimControl,
        settings$control[setdiff(names(settings$control), "calcCytPosGates")]
      )
    ),
    sample = NULL, tubeInfo = function(ind) list(sample = sampleMap$mouse[[as.integer(ind)]]),
    onRecord = function(rec) {
      saveRDS(.acsCytofDebugSlim(rec), file.path(records, paste0(rec$ind, "_", rec$chnl, ".rds")))
      NULL
    }
  ))
}

# Pages for both fits: the CD44+ fit (one StimGate run per cytokine, scored on
# the CD44+ parent) and the non-naive fit (one joint run per strain and
# population, scored on every non-naive cell, combination D). Each rerun must
# reproduce the saved gates and counts exactly.
.omip111DebugRun <- function(prepared, comparison, settings, pathDebug, pathResults,
                             comparisonNonNaive) {
  if (!"package:ggplot2" %in% search()) {
    suppressPackageStartupMessages(attachNamespace("ggplot2"))
  }
  nonNaiveRows <- .omip111FitEval(prepared, comparison, comparisonNonNaive,
    "nonnaive", "nonnaive", "bivariate")
  nonNaiveRows$gateCyt <- comparisonNonNaive$gateCyt[match(
    paste(nonNaiveRows$sampleStim, nonNaiveRows$population, nonNaiveRows$marker, nonNaiveRows$method),
    paste(comparisonNonNaive$sampleStim, comparisonNonNaive$population, comparisonNonNaive$marker,
      comparisonNonNaive$method)
  )]
  nonNaiveRows$errorPp <- 100 * (nonNaiveRows$propBs - nonNaiveRows$manualBs)
  .acsCytofReplaceDir(pathDebug, function(pathTmp) {
    index <- list()
    for (strain in c("C57", "BALB")) {
      sampleMap <- prepared$samples[prepared$samples$strain == strain, ]
      sampleMap$ind <- seq_len(nrow(sampleMap))
      batches <- lapply(split(sampleMap, sampleMap$mouse), function(pair) {
        pair$ind[match(c("uns", "stim"), pair$condition)]
      })
      for (population in c("CD4", "CD8")) {
        markers <- names(.omip111Markers(population))
        # CD44+ fit: one StimGate run per cytokine, as in .omip111Run().
        for (marker in markers) {
          data <- stats::setNames(lapply(sampleMap$sample, function(sid) {
            .omip111MarkerInput(prepared$input[[sid]][[population]], marker)$expression
          }), sampleMap$sample)
          project <- file.path(pathTmp, "work", strain, population, marker, "stimgate")
          records <- file.path(pathTmp, "work", strain, population, marker, "records")
          message("Recording gate diagnostics (CD44+ fit): ", strain, " ", population, " ", marker)
          .omip111DebugGateStim(settings, project, data, batches, marker, sampleMap, records)
          pathSaved <- file.path(pathResults, strain, population, marker, "stimgate")
          .omip111DebugCheckParity(
            stimgate::getStimGates(project), stimgate::getStimGates(pathSaved),
            stimgate::getStimStats(project), stimgate::getStimStats(pathSaved),
            paste(strain, population, marker)
          )
          index[[length(index) + 1L]] <- .omip111DebugPages(pathTmp, records, "cd44", strain,
            population, marker, sampleMap, function(sid, m) {
              comparison[comparison$sampleStim == sid & comparison$population == population &
                comparison$marker == m, ]
            }, prepared, settings,
            fitLabel = "CD44+ parent of the cytokine",
            referenceLabel = "author gate (CD44+ parent)"
          )
          unlink(c(project, records), recursive = TRUE)
        }
        # Non-naive fit: one joint StimGate run, as in .omip111RunNonNaive().
        data <- stats::setNames(lapply(sampleMap$sample, function(sid) {
          input <- prepared$input[[sid]][[population]]
          input$expression[input$nonNaive, markers, drop = FALSE]
        }), sampleMap$sample)
        project <- file.path(pathTmp, "work", strain, population, "non-naive", "stimgate")
        records <- file.path(pathTmp, "work", strain, population, "non-naive", "records")
        message("Recording gate diagnostics (non-naive fit): ", strain, " ", population)
        .omip111DebugGateStim(settings, project, data, batches, markers, sampleMap, records)
        pathSaved <- file.path(pathResults, strain, population, "non-naive", "stimgate")
        gates <- stimgate::getStimGates(project)
        .omip111DebugCheckParity(
          gates, stimgate::getStimGates(pathSaved),
          stimgate::getStimStats(project), stimgate::getStimStats(pathSaved),
          paste(strain, population, "non-naive")
        )
        gates <- gates[gates$gateName == "loc_minClust", ]
        index[[length(index) + 1L]] <- .omip111DebugPages(pathTmp, records, "nonnaive", strain,
          population, markers, sampleMap, function(sid, m) {
            rows <- nonNaiveRows[nonNaiveRows$sampleStim == sid & nonNaiveRows$population == population &
              nonNaiveRows$marker == m, ]
            g <- gates[as.integer(as.character(gates$ind)) == match(sid, sampleMap$sample) &
              gates$marker == m, ]
            rows$locSource <- as.character(g$locSource)
            rows$locReason <- as.character(g$locReason)
            rows
          }, prepared, settings,
          fitLabel = "all non-naive cells",
          referenceLabel = "author gate (non-naive denominator)"
        )
        unlink(c(project, records), recursive = TRUE)
      }
    }
    unlink(file.path(pathTmp, "work"), recursive = TRUE)
    index <- dplyr::bind_rows(index)
    stopifnot(nrow(index) == 160L)
    saveRDS(list(index = index, settings = settings, context = prepared$manifest), file.path(pathTmp, "diagnostics.rds"))
    utils::write.csv(index, file.path(pathTmp, "index.csv"), row.names = FALSE)
    writeLines("complete", file.path(pathTmp, "COMPLETE"))
  })
}
