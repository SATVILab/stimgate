# Step into, and plot, StimGate's local-FDR gating of one simulated sample.
#
# Source after analysis-plot-style.R, with the package loaded by
# devtools::load_all().
#
# `.simDebugLoc()` evaluates a QMD's own "rerun one simulation" call (for
# example `.simBandwidthRunRow(row, ...)` or `.simCompareRunScenario(row, ...)`)
# unchanged, so the data are exactly those of the stored run. It only watches
# from outside, through `trace()`: it records the local-FDR intermediates of
# one stimulated sample and can open the browser in chosen package functions
# for that sample only. No random numbers are drawn by the tracing.

# Package functions whose inputs or outputs are recorded.
.simDebugLocCaptureFns <- c(
  ".getCpUnsLocCondition",
  ".getCpUnsLocGetDensRaw",
  ".getCpUnsLocGetProbFit",
  ".getCpUnsLocGetCp",
  ".getCpUnsLocFilterAfterSmoothing",
  ".getCpUnsLocAntimodeDensity"
)

#' Rerun a simulation and record the local-FDR gating of one or more samples
#'
#' `code` runs in full, unchanged, so every earlier sample and dataset is
#' gated exactly as in the stored run before the target is reached.
#'
#' @param code expression The rerun call, evaluated once after tracing is set
#'   up (e.g. `.simBandwidthRunRow(row_target, ...)`).
#' @param sample integer Sample (donor) number within the simulated dataset.
#'   Its first stimulated tube is gated: GatingSet index
#'   `(sample - 1) * nCondition + 2`. Default: 1, the first sample.
#' @param ind character or NULL GatingSet index of the stimulated tube;
#'   overrides `sample`.
#' @param dataset integer Which time that tube is gated in `code`, i.e. the
#'   simulated dataset (`iter`) when one call gates several datasets.
#'   Default: 1.
#' @param subsequent logical Also record (and browse) every tube gated after
#'   the target, across later datasets too. Default: FALSE.
#' @param nCondition integer Tubes per sample (first is unstimulated).
#' @param browse character Package functions in which to open the browser for
#'   recorded tubes only (e.g. `".getCpUnsLocFilterAfterSmoothing"`). Type `n`
#'   to step, `c` to continue to the next recorded tube, `Q` to quit.
#' @param stopAfter logical Stop the rerun once the last recorded tube is
#'   gated (ignored when `subsequent = TRUE`). Set FALSE to also return the
#'   rerun output (`$result`). Default: TRUE.
#' @return `simDebugLoc` list for one tube, or, when `subsequent = TRUE`, a
#'   `simDebugLocList` named `dataset<k>_ind<i>`, in gating order, with the
#'   rerun output as attribute `"result"`.
.simDebugLoc <- function(
  code,
  sample = 1L,
  ind = NULL,
  dataset = 1L,
  subsequent = FALSE,
  nCondition = 2L,
  browse = character(0L),
  stopAfter = TRUE
) {
  ns <- asNamespace("stimgate")
  ind <- as.character(ind %||% ((as.integer(sample) - 1L) * nCondition + 2L))
  missingFns <- browse[!vapply(
    browse, exists, logical(1L),
    envir = ns, inherits = FALSE
  )]
  if (length(missingFns) > 0L) {
    stop("Not stimgate functions: ", paste(missingFns, collapse = ", "))
  }

  state <- new.env(parent = emptyenv())
  state$seen <- integer(0L)
  state$started <- FALSE
  state$active <- FALSE
  state$tubes <- list()
  doneCondition <- structure(
    class = c("simDebugLocDone", "condition"),
    list(message = "Target sample gated.", call = NULL)
  )

  entryHooks <- list(
    .getCpUnsLocCondition = function(frame) {
      if (!identical(frame$stage, "init")) {
        return(invisible())
      }
      indCurr <- attr(frame$exTblStimNoMin, "ind")
      prev <- if (indCurr %in% names(state$seen)) state$seen[[indCurr]] else 0L
      state$seen[[indCurr]] <- prev + 1L
      isTarget <- identical(indCurr, ind) && state$seen[[indCurr]] == dataset
      if (state$started && !isTRUE(subsequent)) {
        # A plain condition, not an error, so analysis tryCatch(error = )
        # handlers do not turn the early stop into an error row.
        if (isTRUE(stopAfter)) signalCondition(doneCondition)
        return(invisible())
      }
      if (!isTarget && !state$started) {
        return(invisible())
      }
      state$started <- TRUE
      state$active <- TRUE
      state$capture <- list(
        ind = indCurr,
        sample = (as.integer(indCurr) - 2L) %/% nCondition + 1L,
        dataset = state$seen[[indCurr]],
        inputs = mget(
          c(
            "exTblStimNoMin", "exTblUnsBias", "exTblStimOrig", "exTblUnsOrig",
            "chnlSettings", "bias"
          ),
          envir = frame
        ),
        experiment = state$experiment
      )
    },
    .getCpUnsLocGetCp = function(frame) {
      if (state$active) state$capture$dataMod <- frame$dataMod
    }
  )
  exitHooks <- list(
    .getCpUnsLocCondition = function(value) {
      if (state$active) {
        state$capture$cp <- value
        state$active <- FALSE
        cap <- state$capture
        nm <- paste0("dataset", cap$dataset, "_ind", cap$ind)
        state$tubes[[nm]] <- .simDebugLocAssemble(cap)
      }
    },
    .getCpUnsLocGetDensRaw = function(value) {
      if (state$active) {
        state$capture$densTblRaw <- c(state$capture$densTblRaw, list(value))
      }
    },
    .getCpUnsLocGetProbFit = function(value) {
      if (state$active) {
        state$capture$probFit <- c(state$capture$probFit, list(value))
      }
    },
    .getCpUnsLocFilterAfterSmoothing = function(value) {
      if (state$active) state$capture$filter <- value
    },
    .getCpUnsLocAntimodeDensity = function(value) {
      if (state$active) state$capture$antimodeDensity <- value
    }
  )

  traced <- character(0L)
  on.exit(
    for (fn in traced) {
      try(suppressMessages(untrace(fn, where = ns)), silent = TRUE)
    },
    add = TRUE
  )
  for (fn in union(.simDebugLocCaptureFns, browse)) {
    args <- list(what = fn, where = ns, print = FALSE)
    tracer <- list()
    if (!is.null(entryHooks[[fn]])) {
      tracer <- c(tracer, list(bquote(.(entryHooks[[fn]])(environment()))))
    }
    if (fn %in% browse) {
      tracer <- c(tracer, list(bquote(if (isTRUE(.(state)$active)) browser())))
    }
    if (length(tracer) > 0L) {
      args$tracer <- as.call(c(as.name("{"), tracer))
    }
    if (!is.null(exitHooks[[fn]])) {
      args$exit <- bquote(.(exitHooks[[fn]])(returnValue()))
    }
    suppressMessages(do.call(trace, args, quote = TRUE))
    traced <- c(traced, fn)
  }
  # Keep the latest simulated experiment for the true population labels.
  if (requireNamespace("simcyto", quietly = TRUE)) {
    simNs <- asNamespace("simcyto")
    suppressMessages(trace(
      "simCytExperiment",
      exit = bquote(assign("experiment", returnValue(), envir = .(state))),
      where = simNs,
      print = FALSE
    ))
    on.exit(
      try(
        suppressMessages(untrace("simCytExperiment", where = simNs)),
        silent = TRUE
      ),
      add = TRUE
    )
  }

  result <- tryCatch(code, simDebugLocDone = function(e) NULL)

  if (length(state$tubes) == 0L) {
    warning(
      "Tube ", ind, " was not gated by local FDR in dataset ", dataset,
      "; it may have had too few cells, or the call gates fewer datasets."
    )
  }
  if (isTRUE(subsequent)) {
    return(structure(
      state$tubes,
      class = c("simDebugLocList", "list"),
      result = result
    ))
  }
  out <- if (length(state$tubes) > 0L) {
    state$tubes[[1L]]
  } else {
    .simDebugLocAssemble(list(ind = ind, sample = sample, dataset = dataset))
  }
  out$result <- result
  out
}

#' Assemble the record of one gated tube
#'
#' @param cap list Raw captures for one tube.
#' @return list of class `simDebugLoc`.
.simDebugLocAssemble <- function(cap) {
  structure(
    list(
      ind = cap$ind,
      sample = cap$sample,
      dataset = cap$dataset,
      found = !is.null(cap$cp),
      inputs = cap$inputs,
      densTblRaw = cap$densTblRaw[[1L]],
      probFit = cap$probFit,
      dataMod = cap$dataMod,
      filter = cap$filter,
      antimodeDensity = cap$antimodeDensity,
      cp = cap$cp,
      truth = .simDebugLocTruth(cap$inputs, cap$experiment)
    ),
    class = c("simDebugLoc", "list")
  )
}

#' True population labels of the target stimulated and unstimulated tubes
#'
#' Labels come from the simulated experiment. They are matched to the gated
#' (single-precision) expression by rank, and only when the simulated and
#' gated expression agree, so that later data changes cannot misalign them.
#'
#' @param inputs list Captured `.getCpUnsLocCondition()` inputs.
#' @param experiment list or NULL `simcyto::simCytExperiment()` output.
#' @return tibble (`condition`, `x`, `label`) or NULL.
.simDebugLocTruth <- function(inputs, experiment) {
  if (is.null(inputs) || is.null(experiment)) {
    return(NULL)
  }
  chnl <- attr(inputs$exTblStimOrig, "chnlCut")
  one <- function(ex, condition) {
    i <- as.integer(attr(ex, "ind"))
    fr <- experiment$flowFrameList[[i]]
    labels <- experiment$labelsList[[i]]
    if (is.null(fr) || !chnl %in% colnames(flowCore::exprs(fr))) {
      return(NULL)
    }
    xSim <- flowCore::exprs(fr)[, chnl]
    xGated <- sort(ex[[chnl]])
    if (
      length(xSim) != length(xGated) || length(labels) != length(xSim) ||
        !isTRUE(all.equal(sort(xSim), xGated, tolerance = 1e-6))
    ) {
      return(NULL)
    }
    tibble::tibble(
      condition = condition,
      x = xGated,
      label = labels[order(xSim)]
    )
  }
  stim <- one(inputs$exTblStimOrig, "stim")
  uns <- one(inputs$exTblUnsOrig, "unstim")
  if (is.null(stim) || is.null(uns)) {
    message("Simulated and gated expression differ; true labels omitted.")
    return(NULL)
  }
  dplyr::bind_rows(stim, uns)
}

#' Named threshold positions for one debugged sample
#'
#' @param dbg simDebugLoc Output of `.simDebugLoc()`.
#' @return tibble (`line`, `x`) of finite positions.
.simDebugLocLines <- function(dbg) {
  final <- dbg$filter$info$final %||% list()
  x <- c(
    minProbXPos = attr(dbg$dataMod, "minProbXPos") %||% NA_real_,
    xClearInit = final$xClearInit %||% NA_real_,
    xDom = final$xDom %||% NA_real_,
    xQual = final$xQual %||% NA_real_,
    xAntimode = final$xAntimode %||% NA_real_,
    xSum = final$xSum %||% NA_real_,
    threshold = dbg$cp$cp %||% NA_real_
  )
  x <- suppressWarnings(as.numeric(x)) |> stats::setNames(names(x))
  tibble::tibble(
    line = factor(names(x), levels = names(x)),
    x = unname(x)
  ) |>
    dplyr::filter(is.finite(.data$x))
}

#' Threshold and true-label summary for one debugged sample
#'
#' Positives use the strict `x > threshold` rule; the unstimulated count uses
#' raw (unbiased) unstimulated expression.
#'
#' @param dbg simDebugLoc Output of `.simDebugLoc()`.
#' @return tibble One row.
.simDebugLocSummary <- function(dbg) {
  cp <- dbg$cp$cp %||% NA_real_
  chnl <- attr(dbg$inputs$exTblStimOrig, "chnlCut")
  xStim <- dbg$inputs$exTblStimOrig[[chnl]]
  xUns <- dbg$inputs$exTblUnsOrig[[chnl]]
  out <- tibble::tibble(
    ind = dbg$ind,
    threshold = cp,
    locGenerated = dbg$cp$locGenerated %||% NA,
    locSource = dbg$cp$locSource %||% NA_character_,
    locReason = dbg$cp$locReason %||% NA_character_,
    filterReason = dbg$filter$info$reason %||% NA_character_,
    nStim = length(xStim),
    nUns = length(xUns),
    propStimEst = mean(xStim > cp),
    propUnsEst = mean(xUns > cp),
    propRespEst = mean(xStim > cp) - mean(xUns > cp)
  )
  if (is.null(dbg$truth)) {
    return(out)
  }
  stim <- dbg$truth[dbg$truth$condition == "stim", ]
  uns <- dbg$truth[dbg$truth$condition == "unstim", ]
  pos <- stim$x > cp
  gp <- stim$label == "gp"
  out |>
    dplyr::mutate(
      propRespTruth = mean(gp) - mean(uns$label == "gp"),
      nGpStim = sum(gp),
      nPosStim = sum(pos),
      nTruePos = sum(pos & gp),
      fdp = if (sum(pos) > 0L) sum(pos & !gp) / sum(pos) else NA_real_,
      sensitivity = if (sum(gp) > 0L) sum(pos & gp) / sum(gp) else NA_real_
    )
}

#' Diagnostic plots for one debugged sample
#'
#' @param dbg simDebugLoc Output of `.simDebugLoc()`.
#' @return named list of ggplot objects (NULL where data are unavailable):
#'   `density` (raw stim/unstim densities), `prob` (raw and smoothed response
#'   probability), `deriv` (probability derivative), `taut` (taut-string
#'   antimode density), `respCells` (expected responding cells per bin) and
#'   `truth` (stimulated expression by true label).
.simDebugLocPlots <- function(dbg) {
  if (!isTRUE(dbg$found)) {
    stop("The target sample was not gated; nothing to plot.")
  }
  lines <- .simDebugLocLines(dbg)
  vlines <- function() {
    geom_vline(
      data = lines,
      aes(xintercept = .data$x, linetype = .data$line, group = .data$line),
      colour = "grey30"
    )
  }
  scaleLines <- scale_linetype_manual(
    values = c(
      minProbXPos = "dotted", xClearInit = "dotdash", xDom = "longdash",
      xQual = "twodash", xAntimode = "dashed", xSum = "solid",
      threshold = "solid"
    ),
    drop = TRUE,
    name = NULL
  )
  style <- function(p) {
    p + vlines() + scaleLines + .analysis_theme() +
      theme(legend.position = "bottom")
  }
  chnl <- attr(dbg$inputs$exTblStimOrig, "chnlCut")
  out <- list()

  if (is.data.frame(dbg$densTblRaw)) {
    densTbl <- dbg$densTblRaw |>
      dplyr::mutate(condition = ifelse(.data$stim == "yes", "stim", "unstim"))
    out$density <- style(
      ggplot(densTbl, aes(.data$xStim, .data$dens, colour = .data$condition)) +
        geom_line() +
        labs(x = chnl, y = "Raw density", colour = NULL)
    )
  }

  dataMod <- dbg$dataMod
  if (is.data.frame(dataMod)) {
    probCols <- intersect(c("probSmooth", "pred"), names(dataMod))
    probTbl <- tibble::tibble(x = dataMod[[chnl]]) |>
      dplyr::bind_cols(tibble::as_tibble(dataMod[probCols])) |>
      tidyr::pivot_longer(dplyr::all_of(probCols), names_to = "type") |>
      dplyr::mutate(type = dplyr::recode(
        .data$type,
        probSmooth = "raw", pred = "smoothed"
      ))
    out$prob <- style(
      ggplot(probTbl, aes(.data$x, .data$value, colour = .data$type)) +
        geom_line() +
        expand_limits(y = c(0, 1)) +
        labs(x = chnl, y = "Response probability", colour = NULL)
    )

    probCol <- dbg$filter$info$probCol %||% "pred"
    derivTbl <- .getCpUnsLocDerivTbl(dataMod, probCol)
    if (is.data.frame(derivTbl)) {
      out$deriv <- style(
        ggplot(derivTbl, aes(.data$x, .data$deriv)) +
          geom_line() +
          labs(x = chnl, y = "Response probability derivative")
      )
    }

    binVec <- attr(dataMod, "binVec")
    xStim <- dbg$inputs$exTblStimNoMin[[chnl]]
    xRange <- range(dataMod[[chnl]], na.rm = TRUE)
    if (length(binVec) > 1L && length(xStim) > 0L) {
      inRange <- xStim >= xRange[1] & xStim <= xRange[2]
      bin <- findInterval(xStim[inRange], binVec, rightmost.closed = TRUE)
      respTbl <- purrr::map(probCols, function(col) {
        prob <- stats::approx(
          dataMod[[chnl]], dataMod[[col]],
          xout = xStim[inRange], ties = mean, rule = 2
        )$y
        tibble::tibble(bin = bin, prob = prob) |>
          dplyr::group_by(.data$bin) |>
          dplyr::summarise(nResp = sum(.data$prob), .groups = "drop") |>
          dplyr::mutate(
            x = (binVec[.data$bin] +
              binVec[pmin(.data$bin + 1L, length(binVec))]) / 2,
            type = col
          )
      }) |>
        dplyr::bind_rows() |>
        dplyr::mutate(type = dplyr::recode(
          .data$type,
          probSmooth = "raw", pred = "smoothed"
        ))
      out$respCells <- style(
        ggplot(respTbl, aes(.data$x, .data$nResp, colour = .data$type)) +
          geom_line() +
          labs(x = chnl, y = "Expected responding cells", colour = NULL)
      )
    }
  }

  if (!is.null(dbg$antimodeDensity)) {
    tautTbl <- tibble::tibble(
      x = dbg$antimodeDensity$x,
      y = dbg$antimodeDensity$y
    )
    out$taut <- style(
      ggplot(tautTbl, aes(.data$x, .data$y)) +
        geom_step() +
        labs(x = chnl, y = "Taut-string density")
    )
  }

  if (!is.null(dbg$truth)) {
    truthStim <- dbg$truth[dbg$truth$condition == "stim", ]
    out$truth <- style(
      ggplot(truthStim, aes(.data$x, fill = .data$label)) +
        geom_histogram(bins = 100L, position = "identity", alpha = 0.6) +
        scale_y_sqrt() +
        labs(x = chnl, y = "Cells (square-root scale)", fill = "True label")
    )
  }
  out
}

#' Settings and results shown beside the diagnostic plots
#'
#' @param dbg simDebugLoc Output of `.simDebugLoc()`.
#' @param row data.frame or NULL The simulation-grid row that was rerun.
#' @param tuningCols character Grid columns that set StimGate tuning rather
#'   than the simulated data; shown with the gating settings.
#' @return named list of tibbles (`name`, `value`): `gating`, `simulation`
#'   and `estimate`.
.simDebugLocInfo <- function(
  dbg,
  row = NULL,
  tuningCols = c("bw", "bias_uns", "bias_uns_setting")
) {
  fmt <- function(x) {
    if (is.null(x) || length(x) == 0L) {
      return(NA_character_)
    }
    x <- x[[1L]]
    if (is.numeric(x)) format(signif(x, 4L)) else as.character(x)
  }
  pct <- function(x) {
    if (is.finite(x)) paste0(format(signif(100 * x, 3L)), "%") else "NA"
  }
  row <- if (is.null(row)) list() else as.list(row)
  cs <- dbg$inputs$chnlSettings %||% list()
  bwUsed <- attr(dbg$dataMod, "locDensityBw") %||%
    attr(dbg$densTblRaw, "locDensityBw")
  gatingNames <- c(
    "bwScope", "bwMtd", "excMin", "cpMin", "minCell", "gateCombn",
    "clusterGates", "calcCytPosGates", "locProbCol", "locMinPeakProb",
    "locEnforceShapeThreshold"
  )
  gating <- tibble::tibble(
    name = c(
      paste0(intersect(tuningCols, names(row)), " (grid)"),
      "bandwidth used", "biasUns used", gatingNames
    ),
    value = c(
      unname(vapply(
        row[intersect(tuningCols, names(row))], fmt, character(1L)
      )),
      fmt(if (is.list(bwUsed)) NULL else bwUsed),
      fmt(dbg$inputs$bias),
      unname(vapply(cs[gatingNames], fmt, character(1L)))
    )
  )

  simCols <- setdiff(names(row), tuningCols)
  simulation <- tibble::tibble(
    name = c("sample", "dataset", "tube (GatingSet index)", simCols),
    value = c(
      fmt(dbg$sample), fmt(dbg$dataset), fmt(dbg$ind),
      unname(vapply(row[simCols], fmt, character(1L)))
    )
  )

  sm <- .simDebugLocSummary(dbg)
  relErr <- if (!is.null(sm$propRespTruth) && is.finite(sm$propRespTruth) &&
    sm$propRespTruth > 0) {
    (sm$propRespEst - sm$propRespTruth) / sm$propRespTruth
  } else {
    NA_real_
  }
  estimate <- tibble::tibble(
    name = c(
      "threshold", "threshold source", "threshold reason", "filter reason",
      "stim cells", "unstim cells", "stim above threshold",
      "unstim above threshold", "response frequency (estimated)",
      "response frequency (true)", "relative error", "true positives",
      "stim positives", "true responders (stim)", "false discovery",
      "sensitivity"
    ),
    value = c(
      fmt(sm$threshold), fmt(sm$locSource), fmt(sm$locReason),
      fmt(sm$filterReason), fmt(sm$nStim), fmt(sm$nUns),
      pct(sm$propStimEst), pct(sm$propUnsEst), pct(sm$propRespEst),
      pct(sm$propRespTruth %||% NA_real_),
      if (is.finite(relErr)) sprintf("%+.1f%%", 100 * relErr) else "NA",
      fmt(sm$nTruePos), fmt(sm$nPosStim), fmt(sm$nGpStim),
      pct(sm$fdp %||% NA_real_), pct(sm$sensitivity %||% NA_real_)
    )
  ) |>
    dplyr::filter(!is.na(.data$value))
  list(gating = gating, simulation = simulation, estimate = estimate)
}

#' Text block listing one table of settings or results
#'
#' @param tbl tibble (`name`, `value`).
#' @param heading character Block heading, shown as its first line.
#' @param nRow integer Rows to allow, so blocks side by side align.
#' @param width integer Characters per line before wrapping.
#' @return ggplot object.
.simDebugLocInfoPanel <- function(tbl, heading, nRow = NULL, width = 45L) {
  lines <- unlist(lapply(seq_len(nrow(tbl)), function(i) {
    wrapped <- strwrap(
      paste0(tbl$name[[i]], ": ", tbl$value[[i]]),
      width = width,
      exdent = 4L
    )
    if (length(wrapped) == 0L) "" else wrapped
  }))
  txt <- tibble::tibble(
    label = c(heading, lines),
    face = c("bold", rep("plain", length(lines)))
  )
  nRow <- max(nRow %||% 0L, nrow(txt))
  txt$y <- nRow - seq_len(nrow(txt)) + 1L
  ggplot(txt, aes(x = 0.03, y = .data$y, label = .data$label)) +
    geom_text(aes(fontface = .data$face), hjust = 0, size = 3.2) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, nRow + 0.5), expand = c(0, 0)) +
    theme_void()
}

#' Arrange the diagnostic plots in one figure
#'
#' @param plots list Output of `.simDebugLocPlots()`.
#' @param info list or NULL Output of `.simDebugLocInfo()`, shown as text
#'   blocks below the plots.
#' @return ggplot object.
.simDebugLocPlotGrid <- function(plots, info = NULL) {
  plots <- Filter(Negate(is.null), plots)
  grid <- cowplot::plot_grid(
    plotlist = plots,
    ncol = 2L,
    labels = LETTERS[seq_along(plots)],
    align = "hv"
  )
  if (is.null(info)) {
    return(grid)
  }
  headings <- c(
    gating = "Gating settings",
    simulation = "Simulation settings",
    estimate = "Threshold and response frequency"
  )
  info <- info[intersect(names(headings), names(info))]
  nLines <- vapply(info, function(tbl) {
    sum(vapply(
      paste0(tbl$name, ": ", tbl$value),
      function(x) length(strwrap(x, width = 45L, exdent = 4L)),
      integer(1L)
    )) + 1L
  }, integer(1L))
  panels <- purrr::imap(info, function(tbl, nm) {
    .simDebugLocInfoPanel(tbl, headings[[nm]], nRow = max(nLines))
  })
  infoRow <- cowplot::plot_grid(plotlist = panels, nrow = 1L)
  nPlotRows <- ceiling(length(plots) / 2)
  cowplot::plot_grid(
    grid, infoRow,
    ncol = 1L,
    rel_heights = c(10 * nPlotRows, 0.42 * max(nLines))
  )
}
