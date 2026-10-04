#' Plot stimulation gates
#'
#' Plot expression densities with saved gates from [gateStim()]. With two
#' markers, also draw hexbin plots (requires the hexbin package).
#'
#' @param ind numeric vector Sample indices to plot.
#' @param .data GatingSet Data used for gating; supplies uncached expression.
#' @param pathProject character Project directory from [gateStim()].
#' @param marker character vector or NULL One or two marker labels to plot;
#'   supply either `marker` or `chnl`. Default: NULL.
#' @param chnl character vector or NULL One or two channels to plot. Default: NULL.
#' @param pop character or NULL Gated population; NULL selects the single saved
#'   population and errors if several exist. Default: NULL.
#' @param indLab character vector or NULL Sample labels, named by index or in
#'   `ind` order. Default: NULL (sample indices).
#' @param axisLab character vector or NULL Axis labels, named by marker/channel
#'   or in their order. Default: NULL (marker/channel names).
#' @param excMin logical Exclude minimum expression values and show densities
#'   scaled by the retained fraction alongside raw densities. Default: TRUE.
#' @param limitsExpand list or NULL Axis limits to expand to, e.g.
#'   `list(x = c(0, 5), y = c(0, 5))`. Default: NULL.
#' @param limitsEqual logical Give bivariate axes equal ranges. Default: FALSE.
#' @param grid logical Arrange plots in a grid. Default: TRUE.
#' @param gridNCol integer Grid columns. Default: 2.
#' @param showGate logical Draw gate lines. Default: TRUE.
#' @param minCell numeric Minimum retained cell count to plot a sample.
#'   Default: 10.
#' @inheritParams getStimExpr
#' @return A ggplot grid if `grid = TRUE`; otherwise a list of bivariate plots
#'   by sample and univariate plots by marker. NULL if no sample meets `minCell`.
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' pathProject <- gateStim(
#'   tempfile("stimgate_"), gs, exampleData$batchList,
#'   marker = exampleData$marker
#' )
#' plotStim(exampleData$batchList[[1]], gs, pathProject,
#'   marker = exampleData$marker[1]
#' )
#' @export
plotStim <- function(
  ind,
  .data,
  pathProject,
  marker = NULL,
  chnl = NULL,
  pop = NULL,
  indLab = NULL,
  axisLab = NULL,
  excMin = TRUE,
  limitsExpand = NULL,
  limitsEqual = FALSE,
  grid = TRUE,
  gridNCol = 2,
  showGate = TRUE,
  minCell = 10,
  bias = FALSE,
  combnExc = NULL,
  chnlGate = NULL,
  markerGate = NULL,
  gateTypeCytPos = "cyt",
  mult = FALSE
) {
  if (is.null(marker) && is.null(chnl)) {
    stop("Must specify one of marker or chnl")
  }
  pop <- as.character(pop %||% .gateGetPop(pathProject))
  if (length(pop) > 1L) {
    stop("Cannot plot gates for multiple populations")
  }
  if (length(pop) == 0L || !nzchar(pop)) {
    stop("No population found for plotting gates")
  }
  # getStimExpr arguments shared by every plot
  exArgs <- list(
    pathProject = pathProject,
    .data = .data,
    pop = pop,
    marker = marker,
    chnl = chnl,
    excMin = excMin,
    bias = bias,
    combnExc = combnExc,
    chnlGate = chnlGate,
    markerGate = markerGate,
    gateTypeCytPos = gateTypeCytPos,
    mult = mult
  )
  pList <- append(
    .plotGateBv(
      ind = ind,
      indLab = indLab,
      marker = marker,
      chnl = chnl,
      pop = pop,
      axisLab = axisLab,
      pathProject = pathProject,
      limitsExpand = limitsExpand,
      limitsEqual = limitsEqual,
      showGate = showGate,
      minCell = minCell,
      exArgs = exArgs
    ),
    .plotGateUv(
      ind = ind,
      indLab = indLab,
      marker = marker,
      chnl = chnl,
      pop = pop,
      excMin = excMin,
      axisLab = axisLab,
      showGate = showGate,
      pathProject = pathProject,
      minCell = minCell,
      exArgs = exArgs
    )
  )
  if (length(pList) == 0L) {
    return(NULL)
  }
  .plotGrid(plot = grid, pList = pList, nCol = gridNCol)
}

#' @keywords internal
.plotGateBv <- function(
  ind,
  indLab,
  marker,
  chnl,
  pop,
  axisLab,
  pathProject,
  limitsExpand,
  limitsEqual,
  showGate,
  minCell,
  exArgs
) {
  if (xor(is.null(marker), is.null(chnl)) && length(marker %||% chnl) == 1L) {
    return(NULL)
  }
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    stop("The 'hexbin' package is required for bivariate plots.")
  }
  varVec <- (chnl %||% marker)[1:2]
  pList <- lapply(seq_along(ind), function(i) {
    exTbl <- do.call(getStimExpr, c(exArgs, list(ind = ind[[i]])))
    if (nrow(exTbl) < minCell) {
      return(NULL)
    }
    plotTbl <- exTbl[, varVec, drop = FALSE]
    colnames(plotTbl) <- c("V1", "V2")
    p <- ggplot(plotTbl, aes(x = V1, y = V2)) + # nolint
      cowplot::theme_cowplot(14) +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")
      ) +
      geom_hex() +
      scale_fill_viridis_c(trans = "log10", name = "Count") +
      cowplot::background_grid(major = "xy") +
      coord_equal()
    p <- .axisLimits(p, limitsExpand = limitsExpand, limitsEqual = limitsEqual)
    p <- .plotAddAxisTitle(p, marker, chnl, axisLab) +
      ggtitle(.plotGetLab(ind[[i]], indLab, i))
    .plotAddGate(p, ind[[i]], marker, chnl, pop, pathProject, showGate)
  }) |>
    stats::setNames(.plotGetLab(ind, indLab))
  pList <- Filter(Negate(is.null), pList)
  if (length(pList) == 0L) {
    return(NULL)
  }
  pList
}

#' @keywords internal
.plotAddAxisTitle <- function(p, val1, val2, valLab) {
  val <- if (!is.null(val1)) {
    val1
  } else {
    val2
  }
  lab <- .plotGetLab(val, valLab)
  p <- p + labs(x = lab[[1]])
  if (length(lab) > 1L) {
    p <- p + labs(y = lab[[2]])
  }
  p
}

#' @keywords internal
.plotGetLab <- function(val, valLab, i = NULL) {
  if (is.null(valLab)) {
    return(val)
  }
  lab <- if (!is.null(names(valLab))) {
    valLab[val]
  } else {
    if (!is.null(i)) valLab[i] else valLab
  }
  lab |> stats::setNames(NULL)
}

#' @keywords internal
.plotAddGate <- function(p, ind, marker, chnl, pop, pathProject, showGate) {
  if (!showGate) {
    return(p)
  }
  chnl <- chnl %||% stimgateMetaReadMarkerLab(pathProject)[marker]
  # only read gates for gated channels, but keep each channel's axis position
  chnlGate <- chnl[chnl %in% .gateGetChnl(pathProject, pop)]
  if (length(chnlGate) == 0L) {
    return(p)
  }
  gateTbl <- getStimGates(pathProject, pop = pop, chnl = chnlGate) |>
    dplyr::group_by(gateName, chnl, marker, ind, batch) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  gateTbl <- gateTbl[gateTbl[["ind"]] %in% ind, ]
  for (i in seq_len(min(2L, length(chnl)))) {
    gates <- gateTbl[["gate"]][gateTbl[["chnl"]] == chnl[i]]
    if (length(gates) == 0L) next
    # Distinct groups retain overlapping lines at identical thresholds.
    p <- p + if (i == 1L) {
      list(
        geom_vline(
          xintercept = gates, group = seq_along(gates), color = "red", alpha = 0.5
        ),
        expand_limits(x = gates * 1.1)
      )
    } else {
      list(
        geom_hline(
          yintercept = gates, group = seq_along(gates), color = "red", alpha = 0.5
        ),
        expand_limits(y = gates * 1.1)
      )
    }
  }
  p
}

#' @keywords internal
.plotGateUv <- function(
  ind,
  indLab,
  marker,
  chnl,
  pop,
  excMin,
  axisLab,
  showGate,
  pathProject,
  minCell,
  exArgs
) {
  varLoop <- if (!is.null(marker)) marker else chnl
  pList <- lapply(varLoop, function(v) {
    exArgs[["marker"]] <- if (!is.null(marker)) v
    exArgs[["chnl"]] <- if (!is.null(chnl)) v
    .plotGateUvMarker(
      ind = ind,
      indLab = indLab,
      marker = exArgs[["marker"]],
      chnl = exArgs[["chnl"]],
      pop = pop,
      excMin = excMin,
      axisLab = axisLab,
      showGate = showGate,
      pathProject = pathProject,
      minCell = minCell,
      exArgs = exArgs
    )
  }) |>
    stats::setNames(.plotGetLab(varLoop, axisLab))
  pList <- Filter(Negate(is.null), pList)
  if (length(pList) == 0L) {
    return(NULL)
  }
  pList
}

#' @keywords internal
.plotGateUvMarker <- function(
  ind,
  indLab,
  marker,
  chnl,
  pop,
  excMin,
  axisLab,
  showGate,
  pathProject,
  minCell,
  exArgs
) {
  if (length(ind) == 0L) {
    return(NULL)
  }
  chnlBw <- chnl %||% tryCatch(
    stimgateMetaReadMarkerLab(pathProject)[marker],
    error = function(e) NULL
  )
  pathBwProject <- if (!is.null(chnlBw)) {
    file.path(
      pathProject,
      "intermediateData",
      "init",
      chnlBw,
      "ind",
      ind[[1]],
      "bwCpUnsLoc.rds"
    )
  } else {
    ""
  }
  bw <- if (nzchar(pathBwProject) && file.exists(pathBwProject)) {
    tryCatch(readRDS(pathBwProject), error = function(e) "nrd0")
  } else {
    "nrd0"
  }
  .var <- if (!is.null(marker)) marker else chnl
  plotTblList <- lapply(seq_along(ind), function(i) {
    exTbl <- do.call(getStimExpr, c(exArgs, list(ind = ind[[i]])))
    if (nrow(exTbl) < minCell) {
      return(NULL)
    }
    densObj <- stats::density(exTbl[[.var]], na.rm = TRUE, bw = bw)
    plotTbl <- tibble::tibble(x = densObj$x, y = densObj$y, type = "raw")
    if (excMin) {
      # rescale the density to the proportion of cells above the minimum
      probGMin <- attr(exTbl, "probGMin")[[1]][[1]][[1]]
      plotTbl <- plotTbl |>
        dplyr::bind_rows(
          tibble::tibble(x = densObj$x, y = densObj$y * probGMin, type = "adj")
        )
    }
    plotTbl[, "ind"] <- as.character(ind[[i]])
    plotTbl[, "indLab"] <- .plotGetLab(as.character(ind[[i]]), indLab, i)
    plotTbl
  })
  plotTblList <- Filter(Negate(is.null), plotTblList)
  if (length(plotTblList) == 0L) {
    return(NULL)
  }
  aesVec <- c("x", "y", if (excMin) "alpha", if (length(ind) > 1L) "colour")
  p <- ggplot(
    Reduce(rbind, plotTblList),
    aes(x = x, y = y, alpha = type, colour = indLab)[aesVec] # nolint
  ) +
    if (excMin) scale_alpha_manual(values = c("raw" = 0.5, "adj" = 1))
  p <- p +
    geom_line() +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "x") +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white")
    )
  p <- .plotAddAxisTitle(p, marker, chnl, axisLab) +
    labs(y = "Density") +
    ggtitle(.plotGetLab(.var, axisLab))
  .plotAddGate(p, ind, marker, chnl, pop, pathProject, showGate)
}

#' @keywords internal
.plotGrid <- function(plot, pList, nCol) {
  if (!plot) {
    return(pList)
  }
  cowplot::plot_grid(
    plotlist = pList,
    ncol = nCol,
    align = "hv"
  ) +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white")
    )
}
