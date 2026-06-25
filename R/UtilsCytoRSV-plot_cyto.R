plotCyto <- function(
  data,
  marker,
  lab = NULL,
  coordEqual = TRUE,
  limitsExpand = NULL,
  limitsEqual = FALSE,
  fontSize = 14,
  excMin = FALSE,
  geomUni = "histogram",
  ...
) {
  # @description Plot a hex-plot with suitable default for a single
  # sample for cytometry (CyTOF, flow) data.
  #
  # @param data dataframe. Columns are markers (or channels - it doesn't matter,
  # as long as you know which), and rows are cells.
  # @param marker character vector. Columns of \code{data} to plot. The first
  # element is plotted on the x-axis.
  # @param lab named character vector. If not \code{NULL}, then
  # the axis titles for marker are selected using it.
  # @param limitsExpand list. If not \code{NULL},
  # then it is (effectively) passed onto \code{ggplot2::limitsExpand} to
  # ensure that certain values are included in the plot (such as, for example, 0
  # if that is the minimum value possible
  # but it may not be plotted).
  # If not named, then
  # must consist of one numeric vector that will
  # then force all values in the numeric value
  # to be included in the plot. If named, then
  # must have names \code{x} and/or \code{y},
  # with the elements again being numeric vectors
  # that must be included in plot.
  # @param limitsEqual logical. If \code{TRUE},
  # then the ranges on the x- and y-axes
  # must be equal. Effectively applied after
  # expand_grid is applied. Default is \code{FALSE}.
  # @param fontSize integer. Font size to be passed on to
  # \code{cowplot::theme_cowplot(fontSize = <fontSize>)}.
  # @param coordEqual logical. If \code{TRUE},
  # then the \code{coordEqual} ggplot2 function is applied to
  # the plot, making units take up the same visual space on the x-
  # and y-axes.
  # Note that this will cause plots to note be able to be aligned
  # using functions like \code{cowplot::plot_grid} and
  # \code{patchwork::align_plots}.
  # Default is \code{TRUE}.
  # @param excMin logical.
  # If \code{TRUE}, then cells with expression equal to the minimum
  # value of one or both of the variables plotted are excluded.
  # Useful for CyTOF data.
  # Default is \code{FALSE}.
  # @param geomUni
  # "density" or "histogram".
  # Specifies ggplot2 geom to use
  # for univariate data.
  # Default is \code{"histogram"}.
  # @param ...
  # arguments passed to \code{ggplot2::geom_hex} (for bivariate data)
  # or whichever geom is used for univariate data.
  #
  # @return A ggplot object showing cytometry data visualization.
  # checks
  .plotCytoCheck(data = data, lab = lab, fontSize = fontSize)

  # prep
  prepList <- .plotCytoPrep(
    marker = marker,
    lab = lab,
    data = data,
    excMin = excMin
  )

  .plotCytoPlot(
    nMarker = prepList$nMarker,
    plotTbl = prepList$plotTbl,
    marker = prepList$marker,
    fontSize = fontSize,
    coordEqual = coordEqual,
    limitsExpand = limitsExpand,
    limitsEqual = limitsEqual,
    geomUni = geomUni,
    ...
  )
}

#' @keywords internal
.plotCytoCheck <- function(data, lab, fontSize) {
  if (!is.data.frame(data)) {
    stop("data must be a dataframe")
  }
  if (!is.null(lab)) {
    if (!is.character(lab) || is.null(names(lab))) {
      stop("lab must be a named character vector (if not NULL)")
    }
  }
  if (!is.numeric(fontSize)) {
    stop("fontSize must be numeric")
  }
  invisible(TRUE)
}

#' @keywords internal
.plotCytoPrep <- function(marker, lab, data, excMin) {
  nMarker <- min(2, length(marker))
  marker <- marker[seq_len(nMarker)]

  plotTbl <- .plotCytoPrepPlotTbl(
    marker = marker,
    data = data,
    nMarker = nMarker,
    excMin = excMin
  )

  # axis labels
  if (!is.null(lab)) {
    marker <- lab[marker]
  }

  list(
    "nMarker" = nMarker,
    "marker" = marker,
    "plotTbl" = plotTbl
  )
}

#' @keywords internal
.plotCytoPrepPlotTbl <- function(marker, data, nMarker, excMin) {
  # plotTbl
  plotTbl <- data[, marker, drop = FALSE]
  colnames(plotTbl) <- c("V1", "V2")[seq_len(nMarker)]
  .plotCytoPrepPlotTblExcMin(
    excMin = excMin,
    plotTbl = plotTbl,
    nMarker = nMarker
  )
}

#' @keywords internal
.plotCytoPrepPlotTblExcMin <- function(excMin, plotTbl, nMarker) {
  if (!excMin) {
    return(plotTbl)
  }
  plotTbl <- plotTbl |>
    dplyr::filter(
      V1 > min(V1) # nolint
    )
  if (nMarker == 2) {
    plotTbl <- plotTbl |>
      dplyr::filter(
        V2 > min(V2) # nolint
      )
  }
  plotTbl
}

#' @keywords internal
.plotCytoPlot <- function(
  nMarker,
  plotTbl,
  marker,
  fontSize,
  coordEqual,
  limitsExpand,
  limitsEqual,
  geomUni,
  ...
) {
  # base plot
  switch(
    nMarker,
    .plotCytoPlotUni(
      geomUni = geomUni,
      plotTbl = plotTbl,
      fontSize = fontSize,
      marker = marker,
      limitsExpand = limitsExpand,
      ...
    ),
    .plotCytoPlotBiv(
      plotTbl = plotTbl,
      fontSize = fontSize,
      marker = marker,
      coordEqual = coordEqual,
      limitsExpand = limitsExpand,
      limitsEqual = limitsEqual,
      ...
    )
  )
}


#' @keywords internal
.plotCytoPlotUni <- function(
  geomUni,
  plotTbl,
  fontSize,
  marker,
  geomUniGg,
  limitsExpand,
  ...
) {
  geomUniGg <- .plotCytoPlotUniGeom(geomUni, ...)

  p <- .plotCytoPlotUniBase(
    plotTbl = plotTbl,
    fontSize = fontSize,
    marker = marker,
    geomUniGg = geomUniGg
  )

  .plotCytoPlotUniAxes(
    p = p,
    limitsExpand = limitsExpand
  )
}

#' @keywords internal
.plotCytoPlotUniGeom <- function(geomUni, ...) {
  switch(
    geomUni,
    "histogram" = do.call(
      ggplot2::geom_histogram,
      list(...)
    ),
    "density" = do.call(
      ggplot2::geom_density,
      list(...)
    ),
    stop("geom_unit value of ", geomUni, " not recognised")
  )
}

#' @keywords internal
.plotCytoPlotUniBase <- function(plotTbl, fontSize, marker, geomUniGg) {
  ggplot(
    # nolint
    plotTbl, # nolint
    aes(x = V1) # nolint
  ) +
    cowplot::theme_cowplot(fontSize) +
    cowplot::background_grid(major = "x") +
    labs(x = marker[1]) + # nolint
    geomUniGg
}


#' @keywords internal
.plotCytoPlotUniAxes <- function(p, limitsExpand) {
  .plotCytoPlotUniAxesExpand(
    p = p,
    limitsExpand = limitsExpand
  )
}

#' @keywords internal
.plotCytoPlotUniAxesExpand <- function(p, limitsExpand) {
  # return now if axisLimits fn not required
  if (is.null(limitsExpand)) {
    return(p)
  }
  axisLimits(
    p = p,
    limitsExpand = limitsExpand
  )
}

#' @keywords internal
.plotCytoPlotBiv <- function(
  plotTbl,
  fontSize,
  marker,
  coordEqual,
  limitsExpand,
  limitsEqual,
  ...
) {
  p <- .plotCytoPlotBivBase(
    plotTbl = plotTbl,
    fontSize = fontSize,
    marker = marker,
    ...
  )

  .plotCytoBivAxes(
    p = p,
    coordEqual = coordEqual,
    limitsExpand = limitsExpand,
    limitsEqual = limitsEqual
  )
}

#' @keywords internal
.plotCytoPlotBivBase <- function(plotTbl, fontSize, marker, ...) {
  ggplot(
    # nolint
    plotTbl,
    aes(x = V1, y = V2) # nolint
  ) +
    cowplot::theme_cowplot(fontSize) +
    theme(
      # nolint
      plot.background = element_rect(fill = "white"), # nolint
      panel.background = element_rect(fill = "white") # nolint
    ) +
    geom_hex(...) + # nolint
    scale_fill_viridis_c(
      # nolint
      trans = "log10",
      name = "Count"
    ) +
    cowplot::background_grid(major = "xy") +
    labs(x = marker[1], y = marker[2]) # nolint
}

#' @keywords internal
.plotCytoBivAxes <- function(p, coordEqual, limitsExpand, limitsEqual) {
  if (coordEqual) {
    p <- p + coord_equal()
  }

  # return now if axisLimits fn not required
  if (is.null(limitsExpand) && !limitsEqual) {
    return(p)
  }

  .plotCytoBivAxesExpand(
    p = p,
    limitsExpand = limitsExpand,
    limitsEqual = limitsEqual
  )
}

#' @keywords internal
.plotCytoBivAxesExpand <- function(p, limitsEqual, limitsExpand) {
  axisLimits(
    p = p,
    limitsExpand = limitsExpand,
    limitsEqual = limitsEqual
  )
}
