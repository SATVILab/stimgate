# Analysis 16 (16-explore-acs-cytof-cytpos.qmd): how each cytokine is
# distributed among ACS CyTOF cells that are positive for another cytokine,
# to inform the cytokine-positive refinement rule. Reads Analysis 9's cached
# GatingSets and StimGate gates; never re-gates or writes to Analysis 9.

# Expression of the six response cytokines for one stimulated tube and its
# control tube, with Analysis 9's final base and cytokine-positive gates.
.acsCytposTubeData <- function(gs, sampleMap, gates, sampleId, stim, channelMap) {
  row <- function(s) {
    out <- sampleMap[sampleMap$SampleID == sampleId & sampleMap$stim == s, ]
    if (nrow(out) != 1L) {
      stop("Expected one ", s, " tube for ", sampleId, "; found ", nrow(out), ".")
    }
    out
  }
  stimRow <- row(stim)
  unsRow <- row("uns")
  expr <- function(ind) {
    ex <- flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[as.integer(ind)]], "root"))
    ex <- ex[, names(channelMap), drop = FALSE]
    colnames(ex) <- unname(channelMap[colnames(ex)])
    as.data.frame(ex)
  }
  g <- gates[gates$gateName == "loc_minClust" & as.character(gates$ind) == stimRow$ind, ]
  g <- g[match(names(channelMap), g$chnl), ]
  list(
    stim = expr(stimRow$ind),
    uns = expr(unsRow$ind),
    gate = stats::setNames(g$gate, channelMap),
    gateCyt = stats::setNames(g$gateCyt, channelMap)
  )
}

# Cells above the base gate (strict `x > gate`) of at least one cytokine other
# than `cyt`. Cytokines without a gate count as negative.
.acsCytposOtherPos <- function(ex, gate, cyt) {
  others <- setdiff(names(gate), cyt)
  pos <- vapply(others, function(m) {
    (ex[[m]] > gate[[m]]) %in% TRUE
  }, logical(nrow(ex)))
  if (is.null(dim(pos))) pos <- matrix(pos, nrow = nrow(ex))
  rowSums(pos) > 0L
}

# One donor's IFNg and TNF values in both tubes, for the hexagon plots.
.acsCytposHexTbl <- function(dat, donor, x = "IFNg", y = "TNF") {
  dplyr::bind_rows(
    data.frame(tube = "stimulated", x = dat$stim[[x]], y = dat$stim[[y]]),
    data.frame(tube = "unstimulated", x = dat$uns[[x]], y = dat$uns[[y]])
  ) |>
    dplyr::mutate(donor = donor)
}

# One donor's values of `cyt` among all cells and among cells positive for
# another cytokine, in both tubes. Exact zeros are dropped (as StimGate drops
# each tube's minimum); `zeroShare` records how many were dropped.
.acsCytposDensityTbl <- function(dat, donor, cyt) {
  one <- function(ex, tube) {
    other <- .acsCytposOtherPos(ex, dat$gate, cyt)
    x <- ex[[cyt]]
    xOther <- x[other]
    dplyr::bind_rows(
      data.frame(tube = rep(tube, length(x)), cells = "all cells", x = x),
      data.frame(
        tube = rep(tube, length(xOther)),
        cells = rep("positive for another cytokine", length(xOther)), x = xOther
      )
    )
  }
  dplyr::bind_rows(one(dat$stim, "stimulated"), one(dat$uns, "unstimulated")) |>
    dplyr::mutate(donor = donor) |>
    dplyr::group_by(.data$donor, .data$tube, .data$cells) |>
    dplyr::mutate(n = dplyr::n(), zeroShare = mean(.data$x <= 0)) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$x > 0)
}

# Facet label per donor with the number of stimulated cells positive for
# another cytokine (before dropping zeros).
.acsCytposDonorLabels <- function(tbl) {
  lab <- tbl |>
    dplyr::filter(.data$tube == "stimulated", .data$cells != "all cells")
  lab <- lab[!duplicated(lab$donor), c("donor", "n")]
  out <- stats::setNames(paste0(lab$donor, " (n other+ = ", lab$n, ")"), lab$donor)
  missing <- setdiff(unique(tbl$donor), names(out))
  if (length(missing)) {
    out <- c(out, stats::setNames(paste0(missing, " (n other+ = 0)"), missing))
  }
  out
}

.acsCytposGateTbl <- function(gates, cyt) {
  dplyr::bind_rows(lapply(names(gates), function(donor) {
    data.frame(
      donor = donor,
      gate = c("base gate", "cytokine-positive gate"),
      x = c(gates[[donor]]$gate[[cyt]], gates[[donor]]$gateCyt[[cyt]])
    )
  })) |>
    dplyr::filter(is.finite(.data$x)) |>
    dplyr::mutate(group = paste(.data$donor, .data$gate))
}


# Hexagon plots of IFNg against TNF for every donor's stimulated and control
# tube, with Analysis 9's base gates.
.acsCytposHexPlot <- function(tbl, gates, x = "IFNg", y = "TNF") {
  tbl$panel <- paste0(tbl$donor, "\n", tbl$tube)
  lines <- dplyr::bind_rows(lapply(names(gates), function(donor) {
    data.frame(donor = donor, gx = gates[[donor]]$gate[[x]], gy = gates[[donor]]$gate[[y]])
  }))
  lines <- dplyr::bind_rows(
    dplyr::mutate(lines, tube = "stimulated"),
    dplyr::mutate(lines, tube = "unstimulated")
  )
  lines$panel <- paste0(lines$donor, "\n", lines$tube)
  ggplot(tbl, aes(x = .data$x, y = .data$y)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(trans = "log10", name = "Cells") +
    geom_vline(data = lines[is.finite(lines$gx), ], aes(xintercept = .data$gx),
      colour = "#D55E00", linewidth = 0.4) +
    geom_hline(data = lines[is.finite(lines$gy), ], aes(yintercept = .data$gy),
      colour = "#D55E00", linewidth = 0.4) +
    facet_wrap(~panel, ncol = 6) +
    .analysis_theme() +
    labs(x = x, y = y)
}

# Kernel densities of each donor's four curves (tube x cell set). Each curve
# uses its own Sheather-Jones bandwidth (nrd0 if that fails), at least `bwMin`:
# it adapts to multimodal subsets without being as narrow as a bandwidth
# chosen for tens of thousands of cells.
.acsCytposDensityCurves <- function(tbl, n = 512L, bwMin = 0.05) {
  dplyr::bind_rows(lapply(split(tbl, tbl$donor), function(d) {
    from <- min(d$x)
    to <- max(d$x)
    dplyr::bind_rows(lapply(split(d, list(d$tube, d$cells), drop = TRUE), function(g) {
      if (nrow(g) < 5L) {
        return(NULL)
      }
      bw <- tryCatch(stats::bw.SJ(g$x), error = function(e) stats::bw.nrd0(g$x))
      bw <- max(bw, bwMin)
      dens <- stats::density(g$x, bw = bw, from = from, to = to, n = n)
      data.frame(donor = g$donor[[1]], tube = g$tube[[1]], cells = g$cells[[1]],
        x = dens$x, y = dens$y, bw = bw)
    }))
  }))
}

# Densities of `cyt` among all cells and among cells positive for another
# cytokine, for every donor's stimulated and control tube.
.acsCytposDensityPlot <- function(tbl, gates, cyt) {
  labels <- .acsCytposDonorLabels(tbl)
  dens <- .acsCytposDensityCurves(tbl)
  gateTbl <- .acsCytposGateTbl(gates, cyt)
  # Each line has its own group, so a cytokine-positive gate equal to the base
  # gate is still drawn on top of it.
  ggplot(dens, aes(x = .data$x, y = .data$y, colour = .data$tube, linetype = .data$cells)) +
    geom_line(linewidth = 0.6) +
    geom_vline(
      data = gateTbl[gateTbl$gate == "base gate", ],
      aes(xintercept = .data$x, group = .data$group),
      colour = "grey20", linetype = "dashed", inherit.aes = FALSE, linewidth = 0.5
    ) +
    geom_vline(
      data = gateTbl[gateTbl$gate != "base gate", ],
      aes(xintercept = .data$x, group = .data$group),
      colour = "#D55E00", inherit.aes = FALSE, linewidth = 0.5
    ) +
    scale_y_sqrt() +
    scale_colour_manual(values = c(stimulated = "#B2182B", unstimulated = "#2166AC"),
      name = "Tube") +
    scale_linetype_manual(values = c("all cells" = "dotted",
      "positive for another cytokine" = "solid"), name = "Cells") +
    facet_wrap(~donor, ncol = 4, scales = "free_y",
      labeller = ggplot2::as_labeller(labels)) +
    .analysis_theme() +
    theme(legend.position = "bottom", legend.box = "vertical") +
    labs(x = paste0(cyt, " (non-zero values)"), y = "Density (square-root scale)")
}
