# Analysis 13: StimGate's local-FDR gating of every ACS CyTOF stimulated tube,
# one diagnostic page per tube and cytokine.
#
# Source after analysis-runtime.R, analysis-plot-style.R, acs_cytof-helper.R,
# acs_cytof-gate.R, acs_cytof-methods.R and sim-debug-loc.R.
#
# Each population is re-gated with Analysis 9's settings and random-number
# stream, inside `.simDebugLoc()`, into this analysis's own cache folder.
# Analysis 9's GatingSets, StimGate projects, comparator results and manual
# comparison are only read. Every recorded tube is written to disk as soon as
# it is gated; once the final (clustered, cytokine-positive-refined) gates are
# known, the pages are drawn one at a time and the records deleted.

# Evaluate Analysis 9's "scientific-settings" chunk, so both analyses use the
# same populations, StimGate settings and seed without copying them.
.acsCytofDebugSettings <- function(qmdPath) {
  lines <- readLines(qmdPath, warn = FALSE)
  start <- grep("^#\\| label: scientific-settings$", lines)
  if (length(start) != 1L) {
    stop("Analysis 9 has no single 'scientific-settings' chunk: ", qmdPath)
  }
  end <- start + which(lines[-seq_len(start)] == "```")[1L]
  code <- lines[(start + 1L):(end - 1L)]
  env <- new.env(parent = globalenv())
  # Stage switches the chunk reads; they do not affect the settings below.
  env$run_preprocessing <- env$run_stimgate <- env$run_comparators <- FALSE
  eval(parse(text = code[!startsWith(code, "#|")]), envir = env)
  list(
    popVec = env$pop_vec,
    bwMtd = env$stimgate_bw_mtd,
    bwScope = env$stimgate_bw_scope,
    locThresholdMethod = env$stimgate_loc_threshold_method,
    biasUns = env$bias_uns_vec_by_pop,
    seed = env$analysis_seed
  )
}

# "all" (or empty) selects everything (NULL); "none" selects nothing.
.acsCytofDebugList <- function(x) {
  x <- unlist(strsplit(paste(x, collapse = ","), "[,[:space:]]+"))
  x <- x[nzchar(x)]
  if (!length(x) || identical(tolower(x), "all")) {
    return(NULL)
  }
  if (identical(tolower(x), "none")) character(0L) else x
}

# SampleIDs whose pages are saved. "manual" keeps the donors with manual gating
# in Analysis 9's comparison; without that comparison every donor is kept.
.acsCytofDebugResolveSamples <- function(samples, comparison) {
  if (!identical(samples, "manual")) {
    return(samples)
  }
  if (is.null(comparison) || nrow(comparison) == 0L) {
    message("No manual comparison is available: saving pages for every sample.")
    return(NULL)
  }
  sort(unique(as.character(comparison$SampleID)))
}

# This analysis's own folders for one population: never inside Analysis 9's.
.acsCytofDebugPaths <- function(pop, pathDebugBase, outputGroup = NULL) {
  dir <- do.call(file.path, as.list(c(pathDebugBase, outputGroup, pop)))
  list(
    dir = dir,
    stimgate = file.path(dir, "stimgate"),
    records = file.path(dir, "records"),
    pages = file.path(dir, "pages"),
    html = file.path(dir, "html"),
    summary = file.path(dir, "summary.rds")
  )
}

# Tube description for `.simDebugLoc(tubeInfo = )`: donor, stimulus and
# population of a GatingSet index, from the ACS preprocessing sample map.
.acsCytofDebugTubeInfo <- function(sampleMap, pop) {
  function(ind) {
    i <- match(as.character(ind), as.character(sampleMap$ind))
    if (is.na(i)) {
      stop("GatingSet index ", ind, " is not in the ACS sample map.")
    }
    list(
      sample = as.character(sampleMap$SampleID[[i]]),
      stim = as.character(sampleMap$stim[[i]]),
      pop = pop
    )
  }
}

# Keep only the gated channel of the recorded expression, so a record stays
# small; the unstimulated cells with biasUns added are not plotted.
.acsCytofDebugSlim <- function(rec) {
  slim <- function(ex) {
    if (is.null(ex)) {
      return(NULL)
    }
    out <- ex[rec$chnl]
    for (a in c("ind", "chnlCut")) attr(out, a) <- attr(ex, a)
    out
  }
  nms <- c("exTblStimNoMin", "exTblStimOrig", "exTblUnsOrig")
  rec$inputs[nms] <- lapply(rec$inputs[nms], slim)
  rec$inputs$exTblUnsBias <- NULL
  rec$result <- NULL
  rec
}

# Stimulated minus unstimulated proportion above `gate` (strict `x > gate`).
.acsCytofDebugPropBs <- function(xStim, xUns, gate) {
  if (length(gate) != 1L || !is.finite(gate)) {
    return(NA_real_)
  }
  mean(xStim > gate) - mean(xUns > gate)
}

# The highest stimulated-tube value at which stimulated minus unstimulated
# proportion above it reaches the manual background-subtracted proportion.
# Derived from the manual frequency only; it is not a manual gate.
.acsCytofDebugManualMatch <- function(xStim, xUns, target) {
  if (length(target) != 1L || !is.finite(target) || target <= 0) {
    return(NA_real_)
  }
  g <- sort(unique(xStim[is.finite(xStim)]), decreasing = TRUE)
  above <- function(x) 1 - findInterval(g, sort(x)) / length(x)
  hit <- which(above(xStim) - above(xUns) >= target)
  if (length(hit)) g[[hit[[1L]]]] else NA_real_
}

# Analysis 9's StimGate rows against manual gating (read only), or NULL.
.acsCytofDebugReadComparison <- function(path, locThresholdMethod) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tbl <- readRDS(path)
  .acsCytofValidateComparisonManifest(tbl, locThresholdMethod)
  tbl |>
    dplyr::filter(as.character(.data$method) == "stimgate") |>
    dplyr::transmute(
      pop = as.character(.data$popCode),
      SampleID = as.character(.data$SampleID),
      stim = as.character(.data$stim),
      cyt = as.character(.data$cyt),
      freqStimManual = .data$freq_stim_man,
      freqUnsManual = .data$freq_uns_man,
      freqBsManual = .data$freq_bs_man,
      freqBsAnalysis9 = .data$freq_bs_auto,
      diffAnalysis9 = .data$diff,
      absDiffAnalysis9 = .data$abs_diff
    )
}

#' Choose the samples shown in the HTML report
#'
#' Within each population, stimulus and cytokine (or, with `byGroup = FALSE`,
#' across all of them together): the `nLargest` samples where StimGate's
#' Analysis 9 frequency differs most from the manual one, the `nClosest` that
#' agree best (among samples with a positive manual frequency when the
#' comparison has one), and `nRandom` of the rest. Without a manual
#' comparison, all are drawn at random.
#'
#' @param candidates data.frame `pop`, `stim`, `cyt`, `SampleID`.
#' @param comparison data.frame or NULL `.acsCytofDebugReadComparison()`.
#' @param nLargest,nRandom,nClosest integer Samples per group.
#' @param seed integer Seed for the random choices (RNG state is restored).
#' @param samples character or NULL Show these SampleIDs instead.
#' @param byGroup logical Choose within each population, stimulus and
#'   cytokine (TRUE) or across all candidates (FALSE).
#' @return tibble of `candidates` columns plus `htmlReason`.
.acsCytofDebugSelectHtml <- function(
  candidates,
  comparison = NULL,
  nLargest = 2L,
  nRandom = 1L,
  nClosest = 1L,
  seed = 20261008L,
  samples = NULL,
  byGroup = TRUE
) {
  keys <- c("pop", "stim", "cyt", "SampleID")
  candidates <- dplyr::distinct(tibble::as_tibble(candidates[keys]))
  if (!nrow(candidates)) {
    return(dplyr::mutate(candidates, htmlReason = character(0L)))
  }
  if (!is.null(samples)) {
    return(candidates |>
      dplyr::filter(.data$SampleID %in% samples) |>
      dplyr::mutate(htmlReason = "requested"))
  }
  joined <- if (is.null(comparison)) {
    NULL
  } else {
    dplyr::left_join(
      candidates,
      comparison[intersect(c(keys, "absDiffAnalysis9", "freqBsManual"), names(comparison))],
      by = keys
    )
  }
  absDiff <- joined$absDiffAnalysis9 %||% rep(NA_real_, nrow(candidates))
  # Close agreement where both are zero says little about the gate.
  closeOk <- if (is.null(joined$freqBsManual)) {
    rep(TRUE, nrow(candidates))
  } else {
    joined$freqBsManual %in% NA | joined$freqBsManual > 0
  }
  pick <- function(x, n) x[sample.int(length(x), min(n, length(x)))]
  groups <- if (isTRUE(byGroup)) {
    split(seq_len(nrow(candidates)), candidates[c("pop", "stim", "cyt")],
      drop = TRUE
    )
  } else {
    list(seq_len(nrow(candidates)))
  }
  out <- .analysis_with_seed(seed, lapply(groups, function(i) {
    ok <- i[is.finite(absDiff[i])]
    if (!length(ok)) {
      chosen <- pick(i, nLargest + nRandom + nClosest)
      return(dplyr::mutate(
        candidates[chosen, ],
        htmlReason = "random (no manual comparison)"
      ))
    }
    largest <- utils::head(ok[order(-absDiff[ok])], nLargest)
    rest <- setdiff(ok, largest)
    rest <- rest[closeOk[rest]]
    closest <- utils::head(rest[order(absDiff[rest])], nClosest)
    random <- pick(setdiff(i, c(largest, closest)), nRandom)
    dplyr::bind_rows(
      dplyr::mutate(candidates[largest, ], htmlReason = "largest error"),
      dplyr::mutate(candidates[random, ], htmlReason = "random"),
      dplyr::mutate(candidates[closest, ], htmlReason = "close agreement")
    )
  }))
  dplyr::bind_rows(out)
}

# Gates drawn on each page beside StimGate's initial local-FDR lines.
.acsCytofDebugLines <- function(row) {
  num <- function(x) if (length(x) == 1L) as.numeric(x) else NA_real_
  gateCyt <- num(row$gateCyt)
  if (isTRUE(all.equal(gateCyt, num(row$gate)))) gateCyt <- NA_real_
  tibble::tibble(
    line = c(
      "final gate", "final gate, other cytokine positive", "Tailgate gate",
      "F-beta gate", "manual-matched value"
    ),
    label = c(
      "final StimGate gate (loc_minClust)",
      "final gate if another cytokine is positive",
      "Tailgate gate (Analysis 9)", "F-beta gate (Analysis 9)",
      "value matching manual frequency (not a manual gate)"
    ),
    x = c(
      num(row$gate), gateCyt, num(row$tailgateGate), num(row$fbetaGate),
      num(row$manualMatchedX)
    ),
    colour = c(
      "#D55E00", "#D55E00", .analysis_method_colours[["tailgate"]],
      .analysis_method_colours[["fbeta"]], "#7B3294"
    ),
    linetype = c("solid", "dashed", "solid", "solid", "longdash"),
    linewidth = c(1.2, 0.8, 0.9, 0.9, 0.9)
  ) |>
    dplyr::filter(is.finite(.data$x))
}

# Common x range: the plotted data, widened only for gates near the data, so
# a distant fallback gate does not squash the plots (its value is in the key).
.acsCytofDebugXlim <- function(rec, lines) {
  x <- c(
    rec$densTblRaw$xStim,
    if (is.data.frame(rec$dataMod)) rec$dataMod[[rec$chnl]],
    rec$antimodeDensity$x
  )
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NULL)
  }
  r <- range(x)
  pad <- 0.25 * diff(r)
  g <- c(.simDebugLocLines(rec)$x, lines$x)
  range(r, g[is.finite(g) & g >= r[[1L]] - pad & g <= r[[2L]] + pad])
}

# Text blocks for one page.
.acsCytofDebugInfo <- function(rec, row) {
  num <- function(x) {
    if (length(x) == 1L && is.finite(x)) .simDebugFmt(x) else NA_character_
  }
  pct <- function(x) {
    if (length(x) == 1L && is.finite(x)) .simDebugPct(x) else NA_character_
  }
  pts <- function(x) {
    if (length(x) == 1L && is.finite(x)) {
      paste0(format(signif(x, 3L)), "%")
    } else {
      NA_character_
    }
  }
  chr <- function(x) {
    if (length(x) == 1L && !is.na(x)) as.character(x) else NA_character_
  }
  fallback <- function(gate, used) {
    g <- num(gate)
    if (!is.na(g) && isTRUE(used)) paste(g, "(fallback)") else g
  }
  info <- .simDebugLocInfo(rec)
  info$estimate <- info$estimate |>
    dplyr::filter(!.data$name %in% c(
      "response frequency (true)", "relative error", "false discovery",
      "sensitivity"
    ))
  info$simulation <- tibble::tibble(
    name = c(
      "population", "stimulation", "cytokine", "channel", "sample (donor)",
      "stimulated tube (GatingSet index)", "unstimulated tube"
    ),
    value = c(
      chr(row$pop), chr(row$stim), chr(row$cyt), chr(row$chnl),
      chr(row$SampleID), chr(row$ind), chr(row$indUns)
    )
  )
  match9 <- if (is.na(row$matchesAnalysis9)) {
    "not available"
  } else if (isTRUE(row$matchesAnalysis9)) "yes" else "NO"
  info$final <- tibble::tibble(
    name = c(
      "final gate (loc_minClust)", "final gate source", "final gate reason",
      "gate if another cytokine positive", "same gates as Analysis 9",
      "stim - unstim above final gate",
      "StimGate frequency scored in Analysis 9", "manual stim",
      "manual unstim", "manual stim - unstim",
      "Analysis 9 StimGate - manual",
      "value matching manual frequency", "Tailgate gate",
      "stim - unstim above Tailgate gate", "F-beta gate",
      "stim - unstim above F-beta gate"
    ),
    value = c(
      num(row$gate), chr(row$locSource), chr(row$locReason),
      num(row$gateCyt), match9, pct(row$propBsFinal),
      pts(row$freqBsAnalysis9), pts(row$freqStimManual),
      pts(row$freqUnsManual), pts(row$freqBsManual),
      if (length(row$diffAnalysis9) == 1L && is.finite(row$diffAnalysis9)) {
        sprintf("%+.3g percentage points", row$diffAnalysis9)
      } else {
        NA_character_
      },
      num(row$manualMatchedX),
      fallback(row$tailgateGate, row$tailgateFallback),
      pct(row$propBsTailgate),
      fallback(row$fbetaGate, row$fbetaFallback),
      pct(row$propBsFbeta)
    )
  ) |>
    dplyr::filter(!is.na(.data$value))
  info
}

# Block headings, in display order.
.acsCytofDebugHeadings <- c(
  simulation = "Sample",
  gating = "Gating settings",
  estimate = "Initial local-FDR gate",
  final = "Final gates and manual gating"
)

#' One page: the 2c diagnostic plots of a recorded tube and cytokine
#'
#' @param rec simDebugLoc Record from `.simDebugLoc()`.
#' @param row data.frame One summary row (final and comparison gates, manual
#'   frequencies).
#' @return ggplot object.
.acsCytofDebugFigure <- function(rec, row) {
  info <- .acsCytofDebugInfo(rec, row)
  headings <- .acsCytofDebugHeadings
  lines <- .acsCytofDebugLines(row)
  # Tubes that stopped early (almost no stimulated expression above the
  # minimum) have no densities to plot.
  if (!is.data.frame(rec$densTblRaw) && !is.data.frame(rec$dataMod)) {
    return(cowplot::plot_grid(
      plotlist = lapply(names(headings), function(nm) {
        .simDebugLocInfoPanel(info[[nm]], headings[[nm]], nRow = 30L)
      }),
      nrow = 1L
    ))
  }
  plots <- .simDebugLocPlots(
    rec,
    xlim = .acsCytofDebugXlim(rec, lines),
    extraLines = lines
  )
  .simDebugLocPlotGrid(plots, info, headings = headings)
}

# Proportions at the final and comparison gates, and the manual-matched value.
.acsCytofDebugGateProps <- function(rec, row) {
  xStim <- rec$inputs$exTblStimOrig[[rec$chnl]]
  xUns <- rec$inputs$exTblUnsOrig[[rec$chnl]]
  list(
    propBsFinal = .acsCytofDebugPropBs(xStim, xUns, row$gate),
    propBsTailgate = .acsCytofDebugPropBs(xStim, xUns, row$tailgateGate),
    propBsFbeta = .acsCytofDebugPropBs(xStim, xUns, row$fbetaGate),
    manualMatchedX = .acsCytofDebugManualMatch(
      xStim, xUns, row$freqBsManual / 100
    )
  )
}

# Draw one multi-page PDF per stimulus and cytokine (one page per sample,
# ordered by SampleID), deleting each record once drawn. Records chosen for
# the HTML report are kept, with their summary row, under `pathHtml`.
.acsCytofDebugWritePages <- function(summary, pathPages, pathHtml) {
  # Workers start without ggplot2 attached; attach it only now, after gating.
  if (!"package:ggplot2" %in% search()) {
    suppressPackageStartupMessages(attachNamespace("ggplot2"))
  }
  summary$page <- NA_integer_
  summary$pdf <- NA_character_
  summary$pageError <- NA_character_
  for (nm in c("propBsFinal", "propBsTailgate", "propBsFbeta", "manualMatchedX")) {
    summary[[nm]] <- NA_real_
  }
  todo <- which(!is.na(summary$record))
  groups <- split(todo, paste(summary$stim[todo], summary$cyt[todo]))
  for (g in groups) {
    g <- g[order(summary$SampleID[g], as.integer(summary$ind[g]))]
    rel <- file.path(summary$stim[[g[[1L]]]], paste0(summary$cyt[[g[[1L]]]], ".pdf"))
    path <- file.path(pathPages, rel)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(path, width = 36 / 2.54, height = 48 / 2.54)
    tryCatch(
      for (k in seq_along(g)) {
        i <- g[[k]]
        rec <- readRDS(summary$record[[i]])
        props <- .acsCytofDebugGateProps(rec, summary[i, ])
        for (nm in names(props)) summary[[nm]][[i]] <- props[[nm]]
        fig <- tryCatch(.acsCytofDebugFigure(rec, summary[i, ]), error = function(e) e)
        if (inherits(fig, "error")) {
          # Keep the page so page numbers stay in SampleID order.
          summary$pageError[[i]] <- conditionMessage(fig)
          fig <- ggplot() +
            annotate("text", x = 0, y = 0, label = paste(
              summary$SampleID[[i]], "figure failed:", conditionMessage(fig)
            )) +
            theme_void()
        }
        print(fig)
        summary$page[[i]] <- k
        summary$pdf[[i]] <- rel
        if (!is.na(summary$htmlReason[[i]])) {
          dir.create(pathHtml, recursive = TRUE, showWarnings = FALSE)
          saveRDS(
            list(rec = rec, row = summary[i, ]),
            file.path(pathHtml, paste0(summary$ind[[i]], "_", summary$chnl[[i]], ".rds"))
          )
        }
        unlink(summary$record[[i]])
      },
      finally = grDevices::dev.off()
    )
  }
  summary$record <- NULL
  summary
}

# Analysis 9's comparator thresholds for one population (read only).
.acsCytofDebugComparatorGates <- function(paths9, method) {
  path <- paths9[[method]]
  obj <- tryCatch(
    .acsCytofReadComparatorCache(path, method),
    error = function(e) NULL
  )
  if (is.null(obj)) {
    return(NULL)
  }
  obj$thresholds |>
    dplyr::transmute(
      ind = as.character(.data$ind), chnl = .data$chnl,
      "{method}Gate" := .data$threshold,
      "{method}Fallback" := .data$thresholdFallbackUsed
    )
}

#' Re-gate one ACS population while recording every tube, then draw its pages
#'
#' @param pop character Population code.
#' @param paths9 list Analysis 9's `.acsCytofPopulationPaths()`; read only.
#' @param pathsDebug list `.acsCytofDebugPaths()`; built in a temporary
#'   sibling and swapped in only on success.
#' @param settings list `.acsCytofDebugSettings()`.
#' @param select list Optional `stim`, `cyt` and `sample` filters for the
#'   pages (NULL: all). Every tube is still gated, and pages for `htmlKeys`
#'   are always drawn.
#' @param comparison data.frame or NULL `.acsCytofDebugReadComparison()`.
#' @param htmlKeys data.frame or NULL `.acsCytofDebugSelectHtml()`.
#' @param manifest list Run provenance stored with the summary.
#' @param rootDir character or NULL Checkout to load.
#' @return list(pop, success = TRUE, summary path).
.acsCytofDebugRunPopulation <- function(
  pop,
  paths9,
  pathsDebug,
  settings,
  select = list(),
  comparison = NULL,
  htmlKeys = NULL,
  manifest = list(),
  rootDir = NULL
) {
  if (!dir.exists(paths9$gs)) {
    stop("Analysis 9's cached GatingSet for '", pop, "' is missing: ", paths9$gs)
  }
  # The same order of steps as Analysis 9's population runner.
  .acsCytofEnsureCurrentCheckout(rootDir)
  # Read-only backend: the cached GatingSet cannot be changed.
  gs <- flowWorkspace::load_gs(paths9$gs, backend_readonly = TRUE)
  preprocessing <- .acsCytofReadPreprocessing(paths9$gs, gs)
  sampleMap <- tibble::as_tibble(preprocessing$sampleMap) |>
    dplyr::mutate(ind = as.character(.data$ind))
  batchList <- .acsCytofBatchList(sampleMap)
  channelMap <- .acsCytofChannelMap()
  htmlId <- if (is.null(htmlKeys)) {
    character(0L)
  } else {
    paste(htmlKeys$SampleID, htmlKeys$stim, htmlKeys$cyt)
  }
  keep <- function(tube, chnl) {
    cyt <- channelMap[[chnl]]
    paste(tube$sample, tube$stim, cyt) %in% htmlId || (
      (is.null(select$stim) || tube$stim %in% select$stim) &&
        (is.null(select$cyt) || cyt %in% select$cyt) &&
        (is.null(select$sample) || tube$sample %in% select$sample))
  }
  indUns <- unlist(lapply(batchList, function(b) {
    stats::setNames(rep(as.character(b[[1L]]), length(b) - 1L), b[-1L])
  }))

  .acsCytofReplaceDir(pathsDebug$dir, function(pathTmp) {
    pathRecords <- file.path(pathTmp, "records")
    pathProject <- file.path(pathTmp, "stimgate")
    dir.create(pathRecords, showWarnings = FALSE)
    onRecord <- function(rec) {
      file <- NA_character_
      if (keep(rec$tube, rec$chnl)) {
        file <- file.path(pathRecords, paste0(rec$ind, "_", rec$chnl, ".rds"))
        saveRDS(.acsCytofDebugSlim(rec), file)
      }
      .simDebugLocSummary(rec) |>
        dplyr::transmute(
          ind = as.character(.data$ind),
          chnl = rec$chnl,
          initialGate = .data$threshold,
          initialLocSource = .data$locSource,
          initialLocReason = .data$locReason,
          filterReason = .data$filterReason,
          nStim = .data$nStim,
          nUns = .data$nUns,
          record = file
        )
    }
    recorded <- .simDebugLoc(
      .acsCytofGateStim(
        gs = gs,
        pathProject = pathProject,
        batchList = batchList,
        biasUns = settings$biasUns[[pop]],
        bwMtd = settings$bwMtd,
        bwScope = settings$bwScope,
        locThresholdMethod = settings$locThresholdMethod
      ),
      sample = NULL,
      tubeInfo = .acsCytofDebugTubeInfo(sampleMap, pop),
      onRecord = onRecord
    )
    recorded <- dplyr::bind_rows(unclass(recorded))

    finalGates <- function(path) {
      stimgate::getStimGates(path) |>
        dplyr::filter(.data$gateName == "loc_minClust") |>
        dplyr::mutate(ind = as.character(.data$ind))
    }
    summary <- finalGates(pathProject) |>
      dplyr::transmute(
        ind = .data$ind, chnl = .data$chnl, gate = .data$gate,
        gateCyt = .data$gateCyt, locGenerated = .data$locGenerated,
        locSource = .data$locSource, locReason = .data$locReason
      ) |>
      dplyr::left_join(
        dplyr::select(sampleMap, "ind", "SampleID", "stim"),
        by = "ind"
      ) |>
      dplyr::mutate(
        pop = .env$pop,
        cyt = unname(channelMap[.data$chnl]),
        indUns = unname(indUns[.data$ind]),
        .before = 1L
      ) |>
      dplyr::left_join(recorded, by = c("ind", "chnl"))
    if (is.null(summary$record)) summary$record <- NA_character_
    summary$recorded <- !is.na(summary$initialGate) | !is.na(summary$record)

    # Gates must equal Analysis 9's saved gates: same data, settings and
    # random-number stream.
    gates9 <- if (dir.exists(paths9$stimgate)) {
      finalGates(paths9$stimgate) |>
        dplyr::select("ind", "chnl", gateAnalysis9 = "gate", gateCytAnalysis9 = "gateCyt")
    }
    if (is.null(gates9)) {
      summary$gateAnalysis9 <- summary$gateCytAnalysis9 <- NA_real_
    } else {
      summary <- dplyr::left_join(summary, gates9, by = c("ind", "chnl"))
    }
    same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
    summary$matchesAnalysis9 <- ifelse(
      is.na(summary$gateAnalysis9) & is.null(gates9), NA,
      same(summary$gate, summary$gateAnalysis9) &
        same(summary$gateCyt, summary$gateCytAnalysis9)
    )

    for (method in c("tailgate", "fbeta")) {
      gatesMethod <- .acsCytofDebugComparatorGates(paths9, method)
      if (is.null(gatesMethod)) {
        summary[[paste0(method, "Gate")]] <- NA_real_
        summary[[paste0(method, "Fallback")]] <- NA
      } else {
        summary <- dplyr::left_join(summary, gatesMethod, by = c("ind", "chnl"))
      }
    }
    manualCols <- c(
      "freqStimManual", "freqUnsManual", "freqBsManual", "freqBsAnalysis9",
      "diffAnalysis9", "absDiffAnalysis9"
    )
    if (is.null(comparison)) {
      summary[manualCols] <- NA_real_
    } else {
      summary <- dplyr::left_join(
        summary,
        comparison[c("pop", "SampleID", "stim", "cyt", manualCols)],
        by = c("pop", "SampleID", "stim", "cyt")
      )
    }
    summary$htmlReason <- NA_character_
    if (!is.null(htmlKeys) && nrow(htmlKeys)) {
      summary <- dplyr::rows_update(
        summary,
        htmlKeys[c("pop", "SampleID", "stim", "cyt", "htmlReason")],
        by = c("pop", "SampleID", "stim", "cyt"),
        unmatched = "ignore"
      )
    }

    summary <- .acsCytofDebugWritePages(
      summary,
      pathPages = file.path(pathTmp, "pages"),
      pathHtml = file.path(pathTmp, "html")
    )
    # The re-gated project and records are large and no longer needed.
    unlink(c(pathProject, pathRecords), recursive = TRUE)
    saveRDS(
      list(manifest = c(manifest, list(pop = pop)), table = summary),
      file.path(pathTmp, "summary.rds")
    )
  })
  list(pop = pop, success = TRUE, error = NULL, summary = pathsDebug$summary)
}
