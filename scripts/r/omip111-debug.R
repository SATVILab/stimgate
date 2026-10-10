# Per-mouse diagnostics using the same recorder and panels as Analysis 13.
.omip111DebugLines <- function(rows, references, cofactor) {
  sg <- rows[rows$method == "StimGate", ]
  cyt <- sg$gateCyt
  if (isTRUE(all.equal(cyt, sg$threshold))) cyt <- NA_real_
  tibble::tibble(
    line = c("final gate", "conditional gate", "F-beta", "Tailgate", "author stim", "author uns"),
    label = c(
      "final StimGate gate", "StimGate: another cytokine positive",
      "F-beta gate", "Tailgate gate", "projected author cutoff: stimulated",
      "projected author cutoff: control"
    ),
    x = c(
      sg$threshold, cyt, rows$threshold[rows$method == "F-beta"],
      rows$threshold[rows$method == "Tailgate"],
      asinh(references$cytokineLowerRaw[match(c(sg$sampleStim, sg$sampleUns), references$sample)] / cofactor)
    ),
    colour = c(
      "#D55E00", "#D55E00", .analysis_method_colours[["fbeta"]],
      .analysis_method_colours[["tailgate"]], "#7B3294", "#7B3294"
    ),
    linetype = c("solid", "dashed", "solid", "solid", "longdash", "dotted"),
    linewidth = c(1.2, 0.8, 0.9, 0.9, 0.9, 0.9)
  ) |>
    dplyr::filter(is.finite(.data$x))
}

.omip111DebugFigure <- function(rec, rows, references, cofactor) {
  sg <- rows[rows$method == "StimGate", ]
  lines <- .omip111DebugLines(rows, references, cofactor)
  info <- .simDebugLocInfo(rec)
  info$estimate <- info$estimate |>
    dplyr::filter(!.data$name %in% c("response frequency (true)", "relative error", "false discovery", "sensitivity"))
  info$simulation <- tibble::tibble(
    name = c("strain / mouse", "population / marker", "stimulated tube", "control tube", "expression scale"),
    value = c(
      paste(sg$strain, sg$mouse), paste(sg$population, sg$marker),
      sg$sampleStim, sg$sampleUns, paste0("asinh(raw / ", cofactor, ")")
    )
  )
  info$final <- tibble::tibble(
    name = c(
      "saved gates reproduced", "StimGate source", "StimGate reason",
      "projected reference: stim / control / net (%)",
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

.omip111DebugRun <- function(prepared, comparison, settings, pathDebug, pathResults) {
  if (!"package:ggplot2" %in% search()) {
    suppressPackageStartupMessages(attachNamespace("ggplot2"))
  }
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
        data <- stats::setNames(lapply(sampleMap$sample, function(sid) {
          prepared$input[[sid]][[population]]$expression
        }), sampleMap$sample)
        project <- file.path(pathTmp, strain, population, "stimgate")
        records <- file.path(pathTmp, strain, population, "records")
        dir.create(records, recursive = TRUE)
        message("Recording gate diagnostics: ", strain, " ", population)
        .analysis_with_seed(settings$seed, .simDebugLoc(
          stimgate::gateStim(
            pathProject = project, .data = data, batchList = batches,
            marker = markers, biasUns = NULL, control = do.call(stimgate::stimControl, settings$control)
          ),
          sample = NULL, tubeInfo = function(ind) list(sample = sampleMap$mouse[[as.integer(ind)]]),
          onRecord = function(rec) {
            saveRDS(.acsCytofDebugSlim(rec), file.path(records, paste0(rec$ind, "_", rec$chnl, ".rds")))
            NULL
          }
        ))
        pathSaved <- file.path(pathResults, strain, population, "stimgate")
        .omip111DebugCheckParity(
          stimgate::getStimGates(project), stimgate::getStimGates(pathSaved),
          stimgate::getStimStats(project), stimgate::getStimStats(pathSaved),
          paste(strain, population)
        )
        for (marker in markers) {
          relPdf <- file.path(strain, population, paste0(marker, ".pdf"))
          grDevices::pdf(file.path(pathTmp, relPdf), width = 36 / 2.54, height = 48 / 2.54)
          tryCatch(
            {
              pairs <- sampleMap[sampleMap$condition == "stim", ]
              for (i in seq_len(nrow(pairs))) {
                rows <- comparison[comparison$sampleStim == pairs$sample[[i]] &
                  comparison$population == population & comparison$marker == marker, ]
                stopifnot(nrow(rows) == 3L)
                rec <- readRDS(file.path(records, paste0(pairs$ind[[i]], "_", marker, ".rds")))
                refs <- prepared$references[prepared$references$population == population & prepared$references$marker == marker, ]
                fig <- .omip111DebugFigure(rec, rows, refs, settings$preprocessing$cofactor)
                print(fig)
                relPng <- file.path(strain, population, paste0(pairs$mouse[[i]], "-", marker, ".png"))
                ggplot2::ggsave(file.path(pathTmp, relPng), fig,
                  width = 36 / 2.54, height = 48 / 2.54, dpi = 110, limitsize = FALSE,
                  bg = "white"
                )
                index[[length(index) + 1L]] <- data.frame(strain, population, marker,
                  mouse = pairs$mouse[[i]], png = relPng, pdf = relPdf, page = i, matchesSavedGates = TRUE
                )
              }
            },
            finally = grDevices::dev.off()
          )
        }
        unlink(c(project, records), recursive = TRUE)
      }
    }
    index <- dplyr::bind_rows(index)
    stopifnot(nrow(index) == 80L)
    saveRDS(list(index = index, settings = settings, context = prepared$manifest), file.path(pathTmp, "diagnostics.rds"))
    utils::write.csv(index, file.path(pathTmp, "index.csv"), row.names = FALSE)
    writeLines("complete", file.path(pathTmp, "COMPLETE"))
  })
}
