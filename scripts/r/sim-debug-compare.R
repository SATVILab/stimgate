# F-beta and Tailgate diagnostics for one debugged sample, shown beside
# StimGate's local-FDR diagnostics in one figure.
#
# Source after sim-compare-freq_bs.R and sim-debug-loc.R.
#
# Both comparison methods are run as in Analyses 7 and 8: on the simulated
# (double-precision) expression, with raw unstimulated cells (no biasUns),
# through `.simCompareFbetaThreshold()` and `.simCompareTailgateThreshold()`.
# The plotted intermediates are recomputed from the same inputs, and the
# Tailgate rule is checked against the cut-point cytoUtils returned.

# Comparison-method settings used by Analyses 7 and 8.
.simDebugFbetaDefaults <- list(
  beta = 0.8, theta = 2, width = 10, numBins = NULL
)
.simDebugTailgateDefaults <- list(
  x = "stim", adjust = 1, bandwidth = NULL, numPeaks = 1, refPeak = 1,
  tol = 1e-2, autoTol = TRUE
)

# Colours of the comparison gates where they are drawn on StimGate's plots.
.simDebugCompareLineStyle <- tibble::tibble(
  line = c("F-beta gate", "Tailgate gate"),
  colour = c("#7B3294", "#8C510A"),
  linetype = c("solid", "solid"),
  linewidth = c(0.9, 0.9)
)

#' Run F-beta and Tailgate on a debugged sample
#'
#' @param dbg simDebugLoc Output of `.simDebugLoc()`.
#' @param fbeta list Settings overriding `.simDebugFbetaDefaults`.
#' @param tailgate list Settings overriding `.simDebugTailgateDefaults`.
#' @param pathFbeta character or NULL Path to `fbeta.py`.
#' @return list with `fbeta` and `tailgate` results; a method that fails
#'   returns `list(error = <message>)`.
.simDebugCompare <- function(
  dbg,
  fbeta = list(),
  tailgate = list(),
  pathFbeta = NULL
) {
  sim <- dbg$sim
  if (is.null(sim)) {
    stop("The simulated expression of this sample was not recorded.")
  }
  run <- function(code) {
    tryCatch(code, error = function(e) list(error = conditionMessage(e)))
  }
  list(
    fbeta = run(.simDebugFbeta(
      sim, utils::modifyList(.simDebugFbetaDefaults, fbeta), pathFbeta
    )),
    tailgate = run(.simDebugTailgate(
      sim, utils::modifyList(.simDebugTailgateDefaults, tailgate)
    ))
  )
}

#' F-beta threshold with its smoothed histograms and scores
#'
#' @param sim list `dbg$sim`.
#' @param settings list F-beta settings.
#' @param pathFbeta character or NULL Path to `fbeta.py`.
#' @return list (`threshold`, `maxFscore`, `curves`, `settings`).
.simDebugFbeta <- function(sim, settings, pathFbeta = NULL) {
  res <- .simCompareFbetaThreshold(
    xUns = sim$uns,
    xStim = sim$stim,
    pathFbeta = pathFbeta,
    beta = settings$beta,
    theta = settings$theta,
    width = settings$width,
    numBins = settings$numBins
  )
  f <- res$fbeta
  settings$numBinsUsed <- settings$numBins %||%
    as.integer(sqrt(max(length(sim$uns), length(sim$stim))))
  list(
    threshold = res$threshold,
    maxFscore = res$thresholdMetric,
    curves = tibble::tibble(
      x = as.numeric(f$pdfx),
      unstim = as.numeric(f$pdfneg),
      stim = as.numeric(f$pdfpos),
      fscore = as.numeric(f$fscores),
      precision = as.numeric(f$precision),
      recall = as.numeric(f$recall)
    ),
    settings = settings
  )
}

#' Tailgate cut-point with its density-derivative rule
#'
#' With the first-derivative method, cytoUtils finds the reference density
#' peak, then the steepest point of the density's right shoulder (the first
#' valley of the derivative right of the peak), then the first point to its
#' right where |derivative| falls below the tolerance. With `autoTol` the
#' tolerance is 1% of the largest |derivative|.
#'
#' @param sim list `dbg$sim`.
#' @param settings list Tailgate settings.
#' @return list (`threshold`, `refPeak`, `shoulder`, `tolUsed`, `derivMax`,
#'   `bandwidth`, `reproduced`, `density`, `curves`, `settings`).
.simDebugTailgate <- function(sim, settings) {
  x <- switch(settings$x,
    stim = sim$stim,
    unstim = sim$uns,
    combined = c(sim$uns, sim$stim),
    stop("Unknown Tailgate x: ", settings$x)
  )
  x <- x[is.finite(x)]
  res <- .simCompareTailgateThreshold(
    x = x,
    adjust = settings$adjust,
    bandwidth = settings$bandwidth,
    numPeaks = settings$numPeaks,
    refPeak = settings$refPeak,
    method = "firstDeriv",
    tol = settings$tol,
    side = "right",
    strict = FALSE,
    autoTol = settings$autoTol
  )
  bandwidth <- settings$bandwidth %||%
    suppressWarnings(ks::hpi(x, deriv.order = 1L))
  peaks <- sort(openCyto:::.find_peaks(
    x,
    num_peaks = settings$numPeaks, adjust = settings$adjust
  )[, "x"])
  refPeak <- peaks[min(settings$refPeak, length(peaks))]
  derivOut <- cytoUtils:::.deriv_density(
    x = x, adjust = settings$adjust, deriv = 1, bandwidth = bandwidth
  )
  derivMax <- max(abs(derivOut$y))
  tolUsed <- if (isTRUE(settings$autoTol)) 0.01 * derivMax else settings$tol
  valleys <- openCyto:::.find_valleys(
    x = derivOut$x, y = derivOut$y, adjust = settings$adjust
  )
  shoulder <- sort(valleys[valleys > refPeak])[1]
  check <- derivOut$x[derivOut$x > shoulder & abs(derivOut$y) < tolUsed][1]
  dens <- stats::density(x, adjust = settings$adjust)
  list(
    threshold = res$threshold,
    refPeak = refPeak,
    shoulder = shoulder,
    tolUsed = tolUsed,
    derivMax = derivMax,
    bandwidth = bandwidth,
    reproduced = isTRUE(all.equal(check, res$threshold)),
    density = tibble::tibble(x = dens$x, y = dens$y),
    curves = tibble::tibble(
      x = derivOut$x,
      deriv = derivOut$y,
      ratio = abs(derivOut$y) / derivMax
    ),
    settings = settings
  )
}

#' Estimated against true response frequency for one gate
#'
#' Uses the simulated expression with the strict `x > threshold` rule and raw
#' unstimulated cells, as the comparison analyses do.
#'
#' @param threshold numeric Gate.
#' @param sim list `dbg$sim`.
#' @return tibble (`name`, `value`).
.simDebugGateRows <- function(threshold, sim) {
  pct <- .simDebugPct
  posStim <- sim$stim > threshold
  posUns <- sim$uns > threshold
  gpStim <- sim$labelsStim == "gp"
  gpUns <- sim$labelsUns == "gp"
  est <- mean(posStim) - mean(posUns)
  truth <- mean(gpStim) - mean(gpUns)
  relErr <- if (is.finite(truth) && truth > 0) (est - truth) / truth else NA
  tibble::tibble(
    name = c(
      "stim above gate", "unstim above gate",
      "response frequency (estimated)", "response frequency (true)",
      "relative error", "true positives", "stim positives",
      "true responders (stim)", "false discovery", "sensitivity"
    ),
    value = c(
      pct(mean(posStim)), pct(mean(posUns)), pct(est), pct(truth),
      if (is.finite(relErr)) sprintf("%+.1f%%", 100 * relErr) else "NA",
      .simDebugFmt(sum(posStim & gpStim)), .simDebugFmt(sum(posStim)),
      .simDebugFmt(sum(gpStim)),
      pct(if (sum(posStim) > 0L) sum(posStim & !gpStim) / sum(posStim) else NA),
      pct(if (sum(gpStim) > 0L) sum(posStim & gpStim) / sum(gpStim) else NA)
    )
  )
}

#' Settings and results text for the comparison methods
#'
#' @param compare list Output of `.simDebugCompare()`.
#' @param sim list `dbg$sim`.
#' @return list of tibbles (`name`, `value`): `fbetaSettings`, `fbetaResult`,
#'   `tailgateSettings`, `tailgateResult`.
.simDebugCompareInfo <- function(compare, sim) {
  fmt <- .simDebugFmt
  errorRows <- function(res) {
    tibble::tibble(name = "error", value = res$error)
  }
  fb <- compare$fbeta
  tg <- compare$tailgate
  out <- list()
  if (!is.null(fb$error)) {
    out$fbetaSettings <- errorRows(fb)
  } else {
    atGate <- fb$curves[which.min(abs(fb$curves$x - fb$threshold)), ]
    out$fbetaSettings <- tibble::tibble(
      name = c(
        "beta", "theta (positive/negative density ratio)",
        "moving-average width (bins)", "histogram bins",
        "unstimulated cells", "histogram range"
      ),
      value = c(
        fmt(fb$settings$beta), fmt(fb$settings$theta), fmt(fb$settings$width),
        fmt(fb$settings$numBinsUsed), "raw (no biasUns)", "both tubes"
      )
    )
    out$fbetaResult <- dplyr::bind_rows(
      tibble::tibble(
        name = c(
          "gate", "F-score at gate", "precision at gate", "recall at gate"
        ),
        value = c(
          fmt(fb$threshold), fmt(fb$maxFscore), fmt(atGate$precision),
          fmt(atGate$recall)
        )
      ),
      .simDebugGateRows(fb$threshold, sim)
    )
  }
  if (!is.null(tg$error)) {
    out$tailgateSettings <- errorRows(tg)
  } else {
    out$tailgateSettings <- tibble::tibble(
      name = c(
        "cells", "method", "side", "adjust", "bandwidth (hpi, 1st derivative)",
        "peaks searched", "reference peak", "automatic tolerance",
        "tolerance used (|slope|)", "tolerance / largest |slope|"
      ),
      value = c(
        tg$settings$x, "first derivative", "right", fmt(tg$settings$adjust),
        fmt(tg$bandwidth), fmt(tg$settings$numPeaks), fmt(tg$settings$refPeak),
        fmt(tg$settings$autoTol), fmt(tg$tolUsed), fmt(tg$tolUsed / tg$derivMax)
      )
    )
    out$tailgateResult <- dplyr::bind_rows(
      tibble::tibble(
        name = c(
          "reference peak", "steepest right-shoulder point", "gate",
          "rule reproduces gate"
        ),
        value = c(
          fmt(tg$refPeak), fmt(tg$shoulder), fmt(tg$threshold),
          fmt(tg$reproduced)
        )
      ),
      .simDebugGateRows(tg$threshold, sim)
    )
  }
  out
}

#' Comparison gates to draw on StimGate's plots
#'
#' @param compare list Output of `.simDebugCompare()`.
#' @return tibble of line positions and styles.
.simDebugCompareLines <- function(compare) {
  x <- c(
    "F-beta gate" = compare$fbeta$threshold %||% NA_real_,
    "Tailgate gate" = compare$tailgate$threshold %||% NA_real_
  )
  .simDebugCompareLineStyle |>
    dplyr::mutate(x = unname(x[.data$line])) |>
    dplyr::filter(is.finite(.data$x))
}

#' Plots of the comparison methods
#'
#' @param compare list Output of `.simDebugCompare()`.
#' @param chnl character Channel name for the x-axis label.
#' @return named list of ggplot objects: `fbetaPdf`, `fbetaScore`,
#'   `tailgateDensity` and `tailgateRatio` (NULL for a failed method).
.simDebugComparePlots <- function(compare, chnl) {
  cols <- .simDebugLocColours
  base <- function(p, lines) {
    p +
      geom_vline(
        data = lines,
        aes(xintercept = .data$x, group = .data$line),
        colour = lines$colour, linetype = lines$linetype,
        linewidth = lines$linewidth
      ) +
      .analysis_theme() +
      cowplot::background_grid(major = "xy", minor = "x") +
      theme(legend.position = "bottom")
  }
  out <- list()
  fb <- compare$fbeta
  if (is.null(fb$error)) {
    gate <- tibble::tibble(
      line = "F-beta gate", x = fb$threshold, colour = "#7B3294",
      linetype = "solid", linewidth = 0.9
    ) |>
      dplyr::filter(is.finite(.data$x))
    pdfTbl <- fb$curves |>
      tidyr::pivot_longer(c("stim", "unstim"), names_to = "condition") |>
      dplyr::filter(is.finite(.data$value))
    out$fbetaPdf <- base(
      ggplot(pdfTbl, aes(.data$x, .data$value, colour = .data$condition)) +
        geom_line(linewidth = 0.7) +
        scale_colour_manual(values = cols) +
        scale_y_sqrt() +
        labs(
          x = chnl, y = "Smoothed density (square-root scale)",
          colour = NULL
        ),
      gate
    )
    scoreTbl <- fb$curves |>
      tidyr::pivot_longer(
        c("fscore", "precision", "recall"),
        names_to = "score"
      ) |>
      dplyr::mutate(score = dplyr::recode(.data$score, fscore = "F-score")) |>
      dplyr::filter(is.finite(.data$value))
    out$fbetaScore <- base(
      ggplot(scoreTbl, aes(.data$x, .data$value, colour = .data$score)) +
        geom_line(linewidth = 0.7) +
        scale_colour_manual(values = c(
          "F-score" = "#7B3294", precision = "#1B9E77", recall = "#E6AB02"
        )) +
        expand_limits(y = c(0, 1)) +
        labs(x = chnl, y = "Score if gated here", colour = NULL),
      gate
    )
  }
  tg <- compare$tailgate
  if (is.null(tg$error)) {
    tgLines <- tibble::tibble(
      line = c("reference peak", "steepest right shoulder", "Tailgate gate"),
      x = c(tg$refPeak, tg$shoulder, tg$threshold),
      colour = c("#8C8C8C", "#E69F00", "#8C510A"),
      linetype = c("dotted", "dashed", "solid"),
      linewidth = c(0.7, 0.8, 0.9)
    ) |>
      dplyr::filter(is.finite(.data$x))
    out$tailgateDensity <- base(
      ggplot(tg$density, aes(.data$x, .data$y)) +
        geom_line(linewidth = 0.7, colour = cols[["stim"]]) +
        labs(x = chnl, y = paste0("Density (", tg$settings$x, " cells)")),
      tgLines
    )
    ratioTbl <- tg$curves |>
      dplyr::mutate(ratio = pmax(.data$ratio, 1e-6))
    # The rule searches right of the steepest shoulder point only.
    out$tailgateRatio <- base(
      ggplot(ratioTbl, aes(.data$x, .data$ratio)) +
        annotate(
          "rect",
          xmin = tg$shoulder, xmax = Inf,
          ymin = min(ratioTbl$ratio), ymax = Inf,
          fill = "#F6E8C3", alpha = 0.5
        ) +
        geom_line(linewidth = 0.7) +
        geom_hline(
          yintercept = tg$tolUsed / tg$derivMax,
          colour = "#D95F02", linetype = "longdash", linewidth = 0.8
        ) +
        scale_y_log10() +
        labs(x = chnl, y = "|slope| / largest |slope| (log scale)"),
      tgLines
    )
  }
  attr(out, "fbetaLines") <- if (is.null(fb$error)) {
    tibble::tibble(
      line = "F-beta gate", x = fb$threshold, colour = "#7B3294",
      linetype = "solid", linewidth = 0.9
    )
  }
  attr(out, "tailgateLines") <- if (is.null(tg$error)) {
    tibble::tibble(
      line = c(
        "reference peak", "steepest right shoulder", "Tailgate gate",
        "tolerance"
      ),
      x = c(tg$refPeak, tg$shoulder, tg$threshold, tg$tolUsed / tg$derivMax),
      colour = c("#8C8C8C", "#E69F00", "#8C510A", "#D95F02"),
      linetype = c("dotted", "dashed", "solid", "longdash"),
      linewidth = c(0.7, 0.8, 0.9, 0.8)
    )
  }
  out
}

#' One figure with the simulation settings and each method's diagnostics
#'
#' Sections, top to bottom: simulation settings; StimGate's plots (with the
#' comparison gates added), line key, settings and result; then F-beta and
#' Tailgate, each with its plots, settings and result. All plots share one
#' x-axis range.
#'
#' @param dbg simDebugLoc Output of `.simDebugLoc()`.
#' @param row data.frame or NULL The simulation-grid row that was rerun.
#' @param compare list or NULL Output of `.simDebugCompare()`; NULL shows
#'   StimGate only.
#' @param xlim numeric or NULL Common x-axis limits.
#' @return ggplot object with attribute `height_cm` (for a 30 cm width).
.simDebugFigure <- function(dbg, row = NULL, compare = NULL, xlim = NULL) {
  chnl <- attr(dbg$inputs$exTblStimOrig, "chnlCut")
  info <- .simDebugLocInfo(dbg, row)
  sgPlots <- .simDebugLocPlots(
    dbg,
    extraLines = if (!is.null(compare)) .simDebugCompareLines(compare)
  )
  cmpPlots <- if (!is.null(compare)) {
    .simDebugComparePlots(compare, chnl)
  } else {
    list()
  }
  # The Tailgate tolerance entry is a ratio on the y axis, not an x position.
  tgLines <- attr(cmpPlots, "tailgateLines")
  tgX <- if (is.data.frame(tgLines)) {
    tgLines$x[!startsWith(tgLines$line, "tolerance")]
  }
  if (is.null(xlim)) {
    xlim <- range(
      sgPlots[[1L]]$coordinates$limits$x,
      unlist(lapply(cmpPlots, function(p) {
        x <- rlang::eval_tidy(p$mapping$x, p$data)
        x[is.finite(x)]
      })),
      attr(cmpPlots, "fbetaLines")$x, tgX,
      na.rm = TRUE
    )
  }
  sgLines <- attr(sgPlots, "lines")
  cmpLines <- list(
    fbeta = attr(cmpPlots, "fbetaLines"),
    tailgate = attr(cmpPlots, "tailgateLines")
  )
  # Replace each plot's own range with the shared one.
  setXlim <- function(p) {
    p$coordinates <- coord_cartesian(xlim = xlim)
    p
  }
  sgPlots <- lapply(sgPlots, setXlim)
  cmpPlots <- lapply(cmpPlots, setXlim)

  nTextLines <- function(tbl) {
    sum(vapply(
      paste0(tbl$name, ": ", tbl$value),
      function(x) length(strwrap(x, width = 45L, exdent = 4L)),
      integer(1L)
    )) + 1L
  }
  header <- function(text) {
    ggplot() +
      annotate(
        "text",
        x = 0, y = 0, label = text, hjust = 0, size = 6, fontface = "bold"
      ) +
      scale_x_continuous(limits = c(-0.01, 1), expand = c(0, 0)) +
      theme_void()
  }
  textRow <- function(tbls, headings) {
    n <- max(vapply(tbls, nTextLines, integer(1L)))
    panels <- purrr::map2(tbls, headings, function(tbl, h) {
      .simDebugLocInfoPanel(tbl, h, nRow = n)
    })
    list(
      plot = cowplot::plot_grid(plotlist = panels, nrow = 1L),
      height = 0.45 * n
    )
  }
  plotRows <- function(plots, labels) {
    list(
      plot = cowplot::plot_grid(
        plotlist = plots, ncol = 2L, labels = labels, align = "hv"
      ),
      height = 10 * ceiling(length(plots) / 2)
    )
  }
  pieces <- list()
  add <- function(plot, height) {
    pieces[[length(pieces) + 1L]] <<- list(plot = plot, height = height)
  }

  # Simulation settings, in two columns.
  sim <- info$simulation
  half <- ceiling(nrow(sim) / 2)
  add(header("Simulation"), 1)
  simText <- textRow(
    list(sim[seq_len(half), ], sim[-seq_len(half), ]),
    c("Settings", "")
  )
  add(simText$plot, simText$height)

  # StimGate
  sgPlots <- Filter(Negate(is.null), sgPlots)
  nSg <- length(sgPlots)
  add(header("StimGate"), 1)
  sgRows <- plotRows(sgPlots, LETTERS[seq_len(nSg)])
  add(sgRows$plot, sgRows$height)
  add(.simDebugLocLineKey(sgLines), 0.8 * ceiling(nrow(sgLines) / 4))
  sgText <- textRow(
    list(info$gating, info$estimate),
    c("Settings", "Threshold and response frequency")
  )
  add(sgText$plot, sgText$height)

  if (!is.null(compare)) {
    cmpInfo <- .simDebugCompareInfo(compare, dbg$sim)
    used <- nSg
    for (method in c("fbeta", "tailgate")) {
      name <- c(fbeta = "F-beta", tailgate = "Tailgate")[[method]]
      plots <- Filter(
        Negate(is.null),
        cmpPlots[startsWith(names(cmpPlots), method)]
      )
      add(header(name), 1)
      if (length(plots) > 0L) {
        rows <- plotRows(plots, LETTERS[used + seq_along(plots)])
        used <- used + length(plots)
        add(rows$plot, rows$height)
        lines <- cmpLines[[method]]
        if (is.data.frame(lines) && any(is.finite(lines$x))) {
          lines <- lines[is.finite(lines$x), , drop = FALSE]
          add(.simDebugLocLineKey(lines), 0.8 * ceiling(nrow(lines) / 4))
        }
      }
      tbls <- Filter(Negate(is.null), list(
        cmpInfo[[paste0(method, "Settings")]],
        cmpInfo[[paste0(method, "Result")]]
      ))
      headings <- c("Settings", "Gate and response frequency")
      txt <- textRow(tbls, headings[seq_along(tbls)])
      add(txt$plot, txt$height)
    }
  }

  heights <- vapply(pieces, `[[`, numeric(1L), "height")
  fig <- cowplot::plot_grid(
    plotlist = lapply(pieces, `[[`, "plot"),
    ncol = 1L,
    rel_heights = heights
  )
  attr(fig, "height_cm") <- sum(heights)
  fig
}
