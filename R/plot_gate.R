#' Plot stimulation gate
#'
#' Plot bivariate hex and univariate density plots for batches of samples, along
#' with their gates.
#'
#' @param ind numeric vector. Specifies indices in `.data` to plot.
#' @param .data GatingSet, flowSet, cytoset, flowFrame, cytoframe, character,
#'   list, data.frame or NULL Cytometry input as accepted by [gateStim()], in
#'   the same sample order used for gating. NULL uses saved expression where
#'   supported.
#' @param pathProject character.
#' Path to the project directory used for `gateStim`.
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
#' @examples
#' # Create example data and run gating
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' pathProject <- file.path(dirname(exampleData$pathGs), "stimgate")
#'
#' # Run gating
#' gateStim(
#'   .data = gs,
#'   pathProject = pathProject,
#'   popGate = "root",
#'   batchList = exampleData$batchList,
#'   marker = exampleData$marker
#' )
#'
#' # Create plots
#' if (requireNamespace("hexbin", quietly = TRUE)) {
#'   plots <- plotStim(
#'     ind = exampleData$batchList[[1]], # indices in `gs` to plot
#'     .data = gs, # GatingSet
#'     pathProject = pathProject,
#'     marker = exampleData$marker,
#'     grid = TRUE
#'   )
#' }
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
  if (!is.null(.data)) {
    .data <- .asStimGatingSet(.data, pop)
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
