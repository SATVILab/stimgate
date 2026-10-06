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
      truth = .simDebugLocTruth(cap$inputs, cap$experiment),
      sim = .simDebugLocSim(cap$inputs, cap$experiment)
    ),
    class = c("simDebugLoc", "list")
  )
}

#' Simulated expression and labels of the target tubes
#'
#' These are the values the comparison methods receive in Analyses 7 and 8
#' (before StimGate's single-precision storage).
#'
#' @param inputs list Captured `.getCpUnsLocCondition()` inputs.
#' @param experiment list or NULL `simcyto::simCytExperiment()` output.
#' @return list (`stim`, `uns`, `labelsStim`, `labelsUns`) or NULL.
.simDebugLocSim <- function(inputs, experiment) {
  if (is.null(inputs) || is.null(experiment)) {
    return(NULL)
  }
  chnl <- attr(inputs$exTblStimOrig, "chnlCut")
  get <- function(ex) {
    i <- as.integer(attr(ex, "ind"))
    fr <- experiment$flowFrameList[[i]]
    if (is.null(fr) || !chnl %in% colnames(flowCore::exprs(fr))) {
      return(NULL)
    }
    list(x = as.numeric(flowCore::exprs(fr)[, chnl]), labels = experiment$labelsList[[i]])
  }
  stim <- get(inputs$exTblStimOrig)
  uns <- get(inputs$exTblUnsOrig)
  if (is.null(stim) || is.null(uns)) {
    return(NULL)
  }
  list(
    stim = stim$x, uns = uns$x,
    labelsStim = stim$labels, labelsUns = uns$labels
  )
}

# Format one value for the text blocks.
.simDebugFmt <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_character_)
  }
  x <- x[[1L]]
  if (is.numeric(x)) format(signif(x, 4L)) else as.character(x)
}

# Format one proportion as a percentage for the text blocks.
.simDebugPct <- function(x) {
  if (length(x) == 1L && is.finite(x)) {
    paste0(format(signif(100 * x, 3L)), "%")
  } else {
    "NA"
  }
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

# Colours, line types and widths of the diagnostic plots.
.simDebugLocColours <- c(
  stim = "#C2410C", unstim = "#1D4E89",
  raw = "#A3A3A3", smoothed = "#111111",
  gn = "#BDBDBD", gp = "#C2410C"
)
.simDebugLocLineStyle <- tibble::tibble(
  line = c(
    "minProbXPos", "xClearInit", "xDom", "xQual", "xAntimode", "xSum",
    "threshold"
  ),
  colour = c(
    "#8C8C8C", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2",
    "#000000"
  ),
  linetype = c(
    "dotted", "dotdash", "dashed", "longdash", "twodash", "solid", "solid"
  ),
  linewidth = c(0.6, 0.7, 0.7, 0.7, 0.7, 0.7, 1.1)
)

#' Diagnostic plots for one debugged sample
#'
#' All plots share one x-axis range: by default the widest range shown in
#' any of them.
#'
#' @param dbg simDebugLoc Output of `.simDebugLoc()`.
#' @param xlim numeric or NULL Common x-axis limits, e.g. to zoom in on the
#'   threshold. Default: the widest range of the plotted data.
#' @param extraLines tibble or NULL Further gates to draw (`line`, `x`,
#'   `colour`, `linetype`, `linewidth`), e.g. other methods' thresholds.
#' @return named list of ggplot objects (NULL where data are unavailable):
#'   `density` (raw stim/unstim densities), `prob` (raw and smoothed response
#'   probability), `deriv` (probability derivative), `taut` (taut-string
#'   antimode density), `respCells` (expected responding cells per bin) and
#'   `truth` (stimulated expression by true label).
.simDebugLocPlots <- function(dbg, xlim = NULL, extraLines = NULL) {
  if (!isTRUE(dbg$found)) {
    stop("The target sample was not gated; nothing to plot.")
  }
  lines <- .simDebugLocLines(dbg) |>
    dplyr::left_join(
      dplyr::mutate(
        .simDebugLocLineStyle,
        line = factor(.data$line, levels = levels(.simDebugLocLines(dbg)$line))
      ),
      by = "line"
    )
  if (is.data.frame(extraLines) && nrow(extraLines) > 0L) {
    lines <- dplyr::bind_rows(
      dplyr::mutate(lines, line = as.character(.data$line)),
      dplyr::mutate(extraLines, line = as.character(.data$line))
    ) |>
      dplyr::filter(is.finite(.data$x)) |>
      dplyr::mutate(line = factor(.data$line, levels = unique(.data$line)))
  }
  # Colours and widths are set per line, leaving the colour scale for the
  # plotted data; `.simDebugLocPlotGrid()` draws one key for all plots.
  vlines <- function() {
    geom_vline(
      data = lines,
      aes(xintercept = .data$x, linetype = .data$line, group = .data$line),
      colour = lines$colour,
      linewidth = lines$linewidth
    )
  }
  scaleLines <- scale_linetype_manual(
    values = stats::setNames(lines$linetype, lines$line),
    guide = "none"
  )
  style <- function(p) {
    p + vlines() + scaleLines + .analysis_theme() +
      cowplot::background_grid(major = "xy", minor = "x") +
      theme(legend.position = "bottom")
  }
  chnl <- attr(dbg$inputs$exTblStimOrig, "chnlCut")
  out <- list()

  if (is.data.frame(dbg$densTblRaw)) {
    densTbl <- dbg$densTblRaw |>
      dplyr::mutate(condition = ifelse(.data$stim == "yes", "stim", "unstim"))
    out$density <- style(
      ggplot(densTbl, aes(.data$xStim, .data$dens, colour = .data$condition)) +
        geom_line(linewidth = 0.7) +
        scale_colour_manual(values = .simDebugLocColours) +
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
        geom_line(linewidth = 0.7) +
        scale_colour_manual(values = .simDebugLocColours) +
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
          geom_line(linewidth = 0.7) +
          scale_colour_manual(values = .simDebugLocColours) +
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
        geom_histogram(bins = 100L, position = "identity", alpha = 0.8) +
        scale_fill_manual(values = .simDebugLocColours) +
        scale_y_sqrt() +
        labs(x = chnl, y = "Cells (square-root scale)", fill = "True label")
    )
  }

  if (is.null(xlim)) {
    xlim <- range(unlist(lapply(out, function(p) {
      x <- rlang::eval_tidy(p$mapping$x, p$data)
      x[is.finite(x)]
    })), lines$x)
  }
  out <- lapply(out, function(p) p + coord_cartesian(xlim = xlim))
  attr(out, "lines") <- lines
  out
}

#' Key to the threshold lines drawn on every diagnostic plot
#'
#' @param lines tibble Line positions and styles, as stored on the output of
#'   `.simDebugLocPlots()`.
#' @return ggplot object.
.simDebugLocLineKey <- function(lines) {
  key <- tibble::tibble(
    label = paste0(
      as.character(lines$line), " (",
      vapply(lines$x, function(x) format(signif(x, 4L)), character(1L)), ")"
    ),
    colour = lines$colour,
    linetype = lines$linetype,
    linewidth = lines$linewidth,
    x = (seq_len(nrow(lines)) - 1L) %% 4L,
    y = -((seq_len(nrow(lines)) - 1L) %/% 4L)
  )
  ggplot(key) +
    geom_segment(
      aes(
        x = .data$x + 0.02, xend = .data$x + 0.22, y = .data$y, yend = .data$y,
        colour = .data$colour, linetype = .data$linetype,
        linewidth = .data$linewidth
      )
    ) +
    geom_text(
      aes(x = .data$x + 0.25, y = .data$y, label = .data$label),
      hjust = 0, size = 3.4
    ) +
    scale_colour_identity() +
    scale_linetype_identity() +
    scale_linewidth_identity() +
    scale_x_continuous(limits = c(0, 4), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(key$y) - 0.5, 0.5), expand = c(0, 0)) +
    theme_void()
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
  fmt <- .simDebugFmt
  pct <- .simDebugPct
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
  lines <- attr(plots, "lines")
  plots <- Filter(Negate(is.null), plots)
  grid <- cowplot::plot_grid(
    plotlist = plots,
    ncol = 2L,
    labels = LETTERS[seq_along(plots)],
    align = "hv"
  )
  if (is.data.frame(lines) && nrow(lines) > 0L) {
    grid <- cowplot::plot_grid(
      grid, .simDebugLocLineKey(lines),
      ncol = 1L,
      rel_heights = c(10 * ceiling(length(plots) / 2), ceiling(nrow(lines) / 4))
    )
  }
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
