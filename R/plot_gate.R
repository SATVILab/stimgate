#' Plot stimulation gate
#'
#' Plot bivariate hex and univariate density plots for batches of samples, along
#' with their gates.
#'
#' @param ind numeric vector. Specifies indices in `.data` to plot.
#' @param .data GatingSet. Same GatingSet passed to `stimgate_gate`.
#' @param pathProject character.
#' Path to the project directory used for `stimgate_gate`.
#' @param marker character vector of length one or two. Specifies markers
#' to be plotted. If only one is passed, then only univariate plots are created.
#' @param chnl character vector of length one or two. Specifies channels
#' to be plotted. Ignored if `marker` is provided.
#' @param pop character. Specifies population within GatingSet that
#' gates were calculated on. If `NULL`, defaults to population specified
#' by folder name in `project_path/gates/pop_<pop>`, but throws
#' an error if more than one population is detected (i.e. more
#' than one directory in `gates/`). Default is `NULL`.
#' @param indLab named character vector.
#' Labels for `ind` used in plot.
#' Optional.
#' @param axisLab named character vector.
#' Labels for axis titles, applied to `marker` or `chnl`.
#' Optional.
#' @param excMin Logical.
#' If `TRUE`, excludes the minimum expression values when processing the data.
#' Default is `TRUE`.
#' @param limitsExpand list.
#' Expand the limits of the plot axes.
#' Default is `NULL`.
#' @param limitsEqual Logical.
#' If TRUE, forces equal lengths of the limits.
#' @param grid Logical.
#' If TRUE, arranges the resulting plots in a grid format
#' using `cowplot::plot_grid`.
#' Default is `TRUE`.
#' @param gridNCol Integer.
#' Number of columns in grid layout.
#' @param showGate Logical.
#' If `TRUE`, overlays gate lines on the plots.
#' Default is `TRUE`.
#' @param minCell integer.
#' Minimum number of cells to be plotted.
#' Will skip plots with fewer cells.
#' Default is 10.
#' @inheritParams getStimExpr
#'
#' @return A grid of plots if `grid` is TRUE, otherwise a list of ggplot objects.
#'
#' @import ggplot2
#'
#' @examples
#' # Create example data and run gating
#' exampleData <- get_example_data()
#' gs <- flowWorkspace::load_gs(exampleData$path_gs)
#' pathProject <- file.path(dirname(exampleData$path_gs), "stimgate")
#'
#' # Run gating
#' stimgate::stimgate_gate(
#'   .data = gs,
#'   pathProject = pathProject,
#'   popGate = "root",
#'   batchList = exampleData$batchList,
#'   marker = exampleData$marker
#' )
#'
#' # Create plots
#' plots <- stimgate_plot(
#'   ind = exampleData$batchList[[1]], # indices in `gs` to plot
#'   .data = gs, # GatingSet
#'   pathProject = pathProject,
#'   marker = exampleData$marker,
#'   grid = TRUE
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
  gateTypeSinglePos = "single",
  mult = FALSE,
  gateUnsMethod = "min"
) {
  if (is.null(marker) && is.null(chnl)) {
    stop("Must specify one of marker or chnl")
  }
  pop <- pop %% setdiff(.gateGetPop(pathProject), "")
  if (length(pop) > 1L) {
    stop("Cannot plot gates for multiple populations")
  }
  if (length(pop) == 0L || !nzchar(pop)) {
    stop("No population found for plotting gates")
  }
  pList <- .plotGate(
    ind = ind,
    indLab = indLab,
    .data = .data,
    marker = marker,
    chnl = chnl,
    pop = pop,
    axisLab = axisLab,
    pathProject = pathProject,
    excMin = excMin,
    limitsExpand = limitsExpand,
    limitsEqual = limitsEqual,
    showGate = showGate,
    minCell = minCell,
    bias = bias,
    combnExc = combnExc,
    chnlGate = chnlGate,
    markerGate = markerGate,
    gateTypeCytPos = gateTypeCytPos,
    gateTypeSinglePos = gateTypeSinglePos,
    mult = mult,
    gateUnsMethod = gateUnsMethod
  )
  if (length(pList) == 0L) {
    return(NULL)
  }
  .plotGrid(plot = grid, pList = pList, nCol = gridNCol)
}

#' @keywords internal
.plotGate <- function(
  marker,
  chnl,
  pop,
  ind,
  indLab,
  .data,
  axisLab,
  pathProject,
  excMin,
  limitsExpand,
  limitsEqual,
  showGate,
  minCell,
  bias,
  combnExc,
  chnlGate,
  markerGate,
  gateTypeCytPos,
  gateTypeSinglePos,
  mult,
  gateUnsMethod
) {
  # bv
  pListBv <- .plotGateBv(
    marker = marker,
    chnl = chnl,
    pop = pop,
    ind = ind,
    indLab = indLab,
    .data = .data,
    axisLab = axisLab,
    pathProject = pathProject,
    excMin = excMin,
    limitsExpand = limitsExpand,
    limitsEqual = limitsEqual,
    showGate = showGate,
    minCell = minCell,
    bias = bias,
    combnExc = combnExc,
    chnlGate = chnlGate,
    markerGate = markerGate,
    gateTypeCytPos = gateTypeCytPos,
    gateTypeSinglePos = gateTypeSinglePos,
    mult = mult,
    gateUnsMethod = gateUnsMethod
  )

  # uv
  pListUv <- .plotGateUv(
    ind = ind,
    indLab = indLab,
    .data = .data,
    marker = marker,
    chnl = chnl,
    pop = pop,
    excMin = excMin,
    axisLab = axisLab,
    showGate = showGate,
    pathProject = pathProject,
    minCell = minCell,
    bias = bias,
    combnExc = combnExc,
    chnlGate = chnlGate,
    markerGate = markerGate,
    gateTypeCytPos = gateTypeCytPos,
    gateTypeSinglePos = gateTypeSinglePos,
    mult = mult,
    gateUnsMethod = gateUnsMethod
  )

  pListBv |> append(pListUv)
}

#' @keywords internal
.plotGateBv <- function(
  marker,
  chnl,
  pop,
  ind,
  indLab,
  .data,
  axisLab,
  pathProject,
  excMin,
  limitsExpand,
  limitsEqual,
  showGate,
  minCell,
  bias,
  combnExc,
  chnlGate,
  markerGate,
  gateTypeCytPos,
  gateTypeSinglePos,
  mult,
  gateUnsMethod
) {
  oneChnl <- is.null(marker) && !is.null(chnl) && length(chnl) == 1L
  oneMarker <- is.null(chnl) && !is.null(marker) && length(marker) == 1L
  oneVar <- oneChnl || oneMarker
  if (oneVar) {
    return(NULL)
  }
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    if (interactive()) {
      promptAnswer <- readline(
        prompt = paste0(
          "The 'hexbin' package is required for bivariate plots. ",
          "Do you want to install it now? [y/n]: "
        )
      )
      if (tolower(promptAnswer) != "y") {
        stop("Cannot proceed without installing 'hexbin' package.")
      }
      utils::install.packages("hexbin")
    } else {
      stop("The 'hexbin' package is required but not installed.")
    }
  }
  pList <- lapply(seq_along(ind), function(i) {
    indCurr <- ind[[i]]
    exTbl <- .plotGetExTbl(
      ind = indCurr,
      .data = .data,
      pop = pop,
      marker = marker,
      chnl = chnl,
      excMin = excMin,
      pathProject = pathProject,
      bias = bias,
      combnExc = combnExc,
      chnlGate = chnlGate,
      markerGate = markerGate,
      gateTypeCytPos = gateTypeCytPos,
      gateTypeSinglePos = gateTypeSinglePos,
      mult = mult,
      gateUnsMethod = gateUnsMethod
    )
    if (nrow(exTbl) < minCell) {
      return(NULL)
    }
    p <- plotCyto(
      data = exTbl,
      marker = chnl %||% marker,
      excMin = FALSE,
      limitsExpand = limitsExpand,
      limitsEqual = limitsEqual
    ) +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")
      )
    p <- .plotAddAxisTitle(p, marker, chnl, axisLab)
    p <- .plotAddTitle(p, indCurr, i, indLab)
    p <- .plotAddGate(
      p = p,
      ind = indCurr,
      marker = marker,
      chnl = chnl,
      pop = pop,
      pathProject = pathProject,
      showGate = showGate
    )
    p
  }) |>
    stats::setNames(.plotGetLab(ind, indLab))
  pList <- pList[vapply(pList, Negate(is.null), logical(1))]
  if (length(pList) == 0L) {
    return(NULL)
  }
  pList
}

#' @keywords internal
.plotGetExTbl <- function(
  ind,
  .data,
  pop,
  marker,
  chnl,
  excMin,
  pathProject,
  bias,
  combnExc,
  chnlGate,
  markerGate,
  gateTypeCytPos,
  gateTypeSinglePos,
  mult,
  gateUnsMethod
) {
  lapply(ind, function(indCurr) {
    getStimExpr(
      pathProject = pathProject,
      .data = .data,
      pop = pop,
      ind = indCurr,
      chnl = chnl,
      marker = marker,
      bias = bias,
      excMin = excMin,
      combnExc = combnExc,
      chnlGate = chnlGate,
      markerGate = markerGate,
      gateTypeCytPos = gateTypeCytPos,
      gateTypeSinglePos = gateTypeSinglePos,
      mult = mult,
      gateUnsMethod = gateUnsMethod
    )
  }) |>
    Reduce(rbind, x = _)
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
.plotAddTitle <- function(p, ind, i, indLab) {
  p + ggtitle(.plotGetLab(ind, indLab, i))
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
  marker <- marker %||% stimgateMetaReadchnlLab(pathProject)[chnl]
  chnl <- chnl %||% stimgateMetaReadmarkerLab(pathProject)[marker]
  pop <- pop %||% .gateGetPop(pathProject)
  if (length(pop) > 1L) {
    stop("Cannot plot gates for multiple populations")
  }
  if (length(pop) == 0L || !nzchar(pop)) {
    stop("No population found for plotting gates")
  }
  gateTbl <- .plotGetGateTbl(ind, pop, marker, chnl, pathProject)
  chnlGate <- chnl[chnl %in% .gateGetChnl(pathProject, pop)]
  for (i in seq_along(chnlGate)) {
    gateVec <- gateTbl[["gate"]][gateTbl[["chnl"]] == chnlGate[i]]
    for (j in seq_along(gateVec)) {
      gate <- gateVec[[j]]
      if (i == 1) {
        p <- p +
          geom_vline(
            xintercept = gate,
            color = "red",
            alpha = 0.5
          ) +
          expand_limits(
            x = gate * 1.1
          )
      } else if (i == 2) {
        p <- p +
          geom_hline(
            yintercept = gate,
            color = "red",
            alpha = 0.5
          ) +
          expand_limits(y = gate * 1.1)
      }
    }
  }
  p
}

#' @keywords internal
.plotGetGateTbl <- function(ind, pop, marker, chnl, pathProject) {
  gateTbl <- getStimGates(
    pathProject = pathProject,
    pop = pop,
    chnl = chnl,
    marker = marker
  ) |>
    dplyr::group_by(gateName, chnl, marker, ind, batch) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  gateTbl <- gateTbl[gateTbl[["ind"]] %in% ind, ]
  indVec <- NULL
  if (!is.null(marker)) {
    for (i in seq_along(marker)) {
      indVec[[i]] <- which(gateTbl[["marker"]] == marker[i])
    }
  } else {
    for (i in seq_along(chnl)) {
      indVec[[i]] <- which(gateTbl[["chnl"]] == chnl[i])
    }
  }
  indVec <- indVec |> unlist()
  gateTbl[indVec, ]
}

#' @keywords internal
.plotGateUv <- function(
  ind,
  indLab,
  .data,
  marker,
  chnl,
  pop,
  excMin,
  axisLab,
  showGate,
  pathProject,
  minCell,
  bias,
  combnExc,
  chnlGate,
  markerGate,
  gateTypeCytPos,
  gateTypeSinglePos,
  mult,
  gateUnsMethod
) {
  varLoop <- if (!is.null(marker)) marker else chnl
  pList <- lapply(varLoop, function(v) {
    markerCurr <- if (!is.null(marker)) v else NULL
    chnlCurr <- if (!is.null(chnl)) v else NULL
    .plotGateUvMarker(
      marker = markerCurr,
      chnl = chnlCurr,
      ind = ind,
      .data = .data,
      pop = pop,
      excMin = excMin,
      indLab = indLab,
      axisLab = axisLab,
      showGate = showGate,
      pathProject = pathProject,
      minCell = minCell,
      bias = bias,
      combnExc = combnExc,
      chnlGate = chnlGate,
      markerGate = markerGate,
      gateTypeCytPos = gateTypeCytPos,
      gateTypeSinglePos = gateTypeSinglePos,
      mult = mult,
      gateUnsMethod = gateUnsMethod
    )
  }) |>
    stats::setNames(.plotGetLab(varLoop, axisLab))
  pList <- pList[vapply(pList, Negate(is.null), logical(1))]
  if (length(pList) == 0L) {
    return(NULL)
  }
  pList
}

#' @keywords internal
.plotGateUvMarker <- function(
  marker,
  chnl,
  pop,
  ind,
  .data,
  excMin,
  indLab,
  axisLab,
  showGate,
  pathProject,
  minCell,
  bias,
  combnExc,
  chnlGate,
  markerGate,
  gateTypeCytPos,
  gateTypeSinglePos,
  mult,
  gateUnsMethod
) {
  pathBwProject <- file.path(
    pathProject,
    "intermediateData",
    "init",
    chnl,
    "ind",
    ind[length(ind)],
    "bwCpUnsLoc.rds"
  )
  bw <- tryCatch(readRDS(pathBwProject), error = function(e) "nrd0")
  plotTbl <- .plotGateUvMarkerGetPlotTbl(
    marker = marker,
    chnl = chnl,
    pop = pop,
    ind = ind,
    .data = .data,
    excMin = excMin,
    bw = bw,
    indLab = indLab,
    minCell = minCell,
    pathProject = pathProject,
    bias = bias,
    combnExc = combnExc,
    chnlGate = chnlGate,
    markerGate = markerGate,
    gateTypeCytPos = gateTypeCytPos,
    gateTypeSinglePos = gateTypeSinglePos,
    mult = mult,
    gateUnsMethod = gateUnsMethod
  )
  if (is.null(plotTbl)) {
    return(NULL)
  }
  .plotGateUvMarkerPlot(
    plotTbl = plotTbl,
    excMin = excMin,
    ind = ind,
    indLab = indLab,
    pop = pop,
    marker = marker,
    chnl = chnl,
    axisLab = axisLab,
    showGate = showGate,
    pathProject = pathProject
  )
}

#' @keywords internal
.plotGateUvMarkerGetPlotTbl <- function(
  ind,
  .data,
  marker,
  chnl,
  pop,
  excMin,
  bw,
  indLab,
  minCell,
  pathProject,
  bias,
  combnExc,
  chnlGate,
  markerGate,
  gateTypeCytPos,
  gateTypeSinglePos,
  mult,
  gateUnsMethod
) {
  plotTblList <- lapply(seq_along(ind), function(i) {
    plotTbl <- .plotGateUvMarkerGetPlotTblInd(
      ind = ind[[i]],
      .data = .data,
      pop = pop,
      marker = marker,
      chnl = chnl,
      excMin = excMin,
      bw = bw,
      minCell = minCell,
      pathProject = pathProject,
      bias = bias,
      combnExc = combnExc,
      chnlGate = chnlGate,
      markerGate = markerGate,
      gateTypeCytPos = gateTypeCytPos,
      gateTypeSinglePos = gateTypeSinglePos,
      mult = mult,
      gateUnsMethod = gateUnsMethod
    )
    if (is.null(plotTbl)) {
      return(NULL)
    }
    plotTbl[, "ind"] <- as.character(ind[[i]])
    plotTbl[, "indLab"] <- .plotGetLab(as.character(ind[[i]]), indLab, i)
    plotTbl
  })
  plotTblList <- plotTblList[
    vapply(plotTblList, Negate(is.null), logical(1))
  ]
  if (length(plotTblList) == 0L) {
    return(NULL)
  }
  Reduce(rbind, plotTblList)
}

#' @keywords internal
.plotGateUvMarkerGetPlotTblInd <- function(
  ind,
  .data,
  marker,
  chnl,
  pop,
  excMin,
  bw,
  minCell,
  pathProject,
  bias,
  combnExc,
  chnlGate,
  markerGate,
  gateTypeCytPos,
  gateTypeSinglePos,
  mult,
  gateUnsMethod
) {
  exTbl <- getStimExpr(
    pathProject = pathProject,
    .data = .data,
    pop = pop,
    ind = ind,
    chnl = chnl,
    marker = marker,
    bias = bias,
    excMin = excMin,
    combnExc = combnExc,
    chnlGate = chnlGate,
    markerGate = markerGate,
    gateTypeCytPos = gateTypeCytPos,
    gateTypeSinglePos = gateTypeSinglePos,
    mult = mult,
    gateUnsMethod = gateUnsMethod
  )
  if (nrow(exTbl) < minCell) {
    return(NULL)
  }
  .var <- if (!is.null(marker)) marker else chnl
  densObjRaw <- density(exTbl[[.var]], na.rm = TRUE, bw = bw)
  plotTbl <- tibble::tibble(x = densObjRaw$x, y = densObjRaw$y)
  .plotGateUvMarkerAddAdj(
    excMin = excMin,
    plotTbl = plotTbl,
    densObjRaw = densObjRaw,
    exTbl = exTbl
  )
}

#' @keywords internal
.plotGateUvMarkerAddAdj <- function(
  excMin,
  plotTbl,
  densObjRaw,
  exTbl
) {
  if (!excMin) {
    return(NULL)
  }
  probGMin <- attr(exTbl, "probGMin")[[1]][[1]][[1]]
  plotTbl[, "type"] <- "raw"
  densObjAdj <- densObjRaw
  densObjAdj$y <- densObjAdj$y * probGMin
  plotTblAdj <- tibble::tibble(
    x = densObjAdj$x,
    y = densObjAdj$y,
    type = "adj"
  )
  plotTbl |>
    dplyr::bind_rows(plotTblAdj)
}

#' @keywords internal
.plotGateUvMarkerPlot <- function(
  plotTbl,
  excMin,
  ind,
  indLab,
  marker,
  chnl,
  pop,
  axisLab,
  showGate,
  pathProject
) {
  p <- .plotGateUvMarkerPlotInit(plotTbl, excMin, ind, indLab)
  p <- .plotAddAxisTitle(p, marker, chnl, axisLab)
  p <- p + ggplot2::labs(y = "Density")
  .var <- if (!is.null(marker)) marker else chnl
  p <- .plotAddTitle(p, .var, NULL, axisLab)
  p <- .plotAddGate(p, ind, marker, chnl, pop, pathProject, showGate)
  p
}

#' @keywords internal
.plotGateUvMarkerPlotInit <- function(plotTbl, excMin, ind, indLab) {
  alphaLabVec <- c("raw" = 0.5, "adj" = 1)
  p <- if (excMin) {
    if (length(ind) > 1L) {
      ggplot(plotTbl, aes(x = x, y = y, alpha = type, color = indLab)) +
        scale_alpha_manual(values = alphaLabVec)
    } else {
      ggplot(plotTbl, aes(x = x, y = y, alpha = type)) +
        scale_alpha_manual(values = alphaLabVec)
    }
  } else {
    if (length(ind) > 1L) {
      ggplot(plotTbl, aes(x = x, y = y, colour = indLab)) +
        scale_alpha_manual(values = alphaLabVec)
    } else {
      ggplot(plotTbl, aes(x = x, y = y)) +
        scale_alpha_manual(values = alphaLabVec)
    }
  }
  p +
    geom_line() +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "x") +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white")
    )
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
