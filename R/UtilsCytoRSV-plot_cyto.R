#' @title Plot cytometry expression data
#'
#' @description Plot a hex-plot with suitable default for a single
#' sample for cytometry (CyTOF, flow) data.
#'
#' @param data dataframe. Columns are markers (or channels - it doesn't matter,
#' as long as you know which), and rows are cells.
#' @param marker character vector. Columns of \code{data} to plot. The first
#' element is plotted on the x-axis.
#' @param lab named character vector. If not \code{NULL}, then
#' the axis titles for marker are selected using it.
#' @param limits_expand list. If not \code{NULL},
#' then it is (effectively) passed onto \code{ggplot2::limits_expand} to
#' ensure that certain values are included in the plot (such as, for example, 0
#' if that is the minimum value possible
#' but it may not be plotted).
#' If not named, then
#' must consist of one numeric vector that will
#' then force all values in the numeric value
#' to be included in the plot. If named, then
#' must have names \code{x} and/or \code{y},
#' with the elements again being numeric vectors
#' that must be included in plot.
#' @param limits_equal logical. If \code{TRUE},
#' then the ranges on the x- and y-axes
#' must be equal. Effectively applied after
#' expand_grid is applied. Default is \code{FALSE}.
#' @param font_size integer. Font size to be passed on to
#' \code{cowplot::theme_cowplot(font_size = <font_size>)}.
#' @param coord_equal logical. If \code{TRUE},
#' then the \code{coord_equal} ggplot2 function is applied to
#' the plot, making units take up the same visual space on the x-
#' and y-axes.
#' Note that this will cause plots to note be able to be aligned
#' using functions like \code{cowplot::plot_grid} and
#' \code{patchwork::align_plots}.
#' Default is \code{TRUE}.
#' @param exc_min logical.
#' If \code{TRUE}, then cells with expression equal to the minimum
#' value of one or both of the variables plotted are excluded.
#' Useful for CyTOF data.
#' Default is \code{FALSE}.
#' @param geom_uni
#' "density" or "histogram".
#' Specifies ggplot2 geom to use
#' for univariate data.
#' Default is \code{"histogram"}.
#' @param ... arguments passed to \code{ggplot2::geom_hex} (for bivariate data)
#' or whichever geom is used for univariate data.
#'
#' @import ggplot2
#'
#' @details Uses the following defaults:
#' \itemize{
#'    \item Plots a hex plot if \code{marker} param has length 2
#'    \item Uses theme \code{cowplot::theme_cowplot} with background grid
#'    \code{cowplot::background_grid(major = 'xy')}.
#'    \item If a hex plot, then fill colour is defined by
#'    \code{scale_fill_viridis_c} with params \code{trans = 'log10'}.
#'    \item axis units are equal between the axes
#' }
#'
#' @import ggplot2
#'
#' @examples
#' if (!requireNamespace("flowCore", quietly = TRUE)) {
#'   if (!requireNamespace("BiocManager", quietly = TRUE)) {
#'     install.packages("BiocManager")
#'   }
#'   BiocManager::install("flowCore")
#' }
#' data("GvHD", package = "flowCore")
#' ex_tbl <- flowCore::exprs(GvHD[[1]])
#' marker <- c("FL2-H", "FL3-H")
#' plot_cyto(
#'   data = ex_tbl,
#'   marker = marker
#' )
#' lab_vec <- chnl_lab(data = GvHD)
#' plot_cyto(
#'   data = ex_tbl,
#'   marker = marker,
#'   lab = lab_vec
#' )
plot_cyto <- function(data, marker, lab = NULL,
                      coord_equal = TRUE,
                      limits_expand = NULL, limits_equal = FALSE,
                      font_size = 14, exc_min = FALSE,
                      geom_uni = "histogram", ...) {
  # checks
  .plot_cyto_check(data = data, lab = lab, font_size = font_size)

  # prep
  prep_list <- .plot_cyto_prep(
    marker = marker, lab = lab,
    data = data, exc_min = exc_min
  )

  .plot_cyto_plot(
    n_marker = prep_list$n_marker,
    plot_tbl = prep_list$plot_tbl,
    marker = prep_list$marker,
    font_size = font_size,
    coord_equal = coord_equal,
    limits_expand = limits_expand,
    limits_equal = limits_equal,
    geom_uni = geom_uni,
    ...
  )
}

.plot_cyto_check <- function(data, lab, font_size) {
  if (!is.data.frame(data)) stop("data must be a dataframe")
  if (!is.null(lab)) {
    if (!is.character(lab) || is.null(names(lab))) {
      stop("lab must be a named character vector (if not NULL)")
    }
  }
  if (!is.numeric(font_size)) {
    stop("font_size must be numeric")
  }
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    .install_pkg_cran("hexbin") # nolint
  }
  invisible(TRUE)
}

.plot_cyto_prep <- function(marker,
                            lab,
                            data,
                            exc_min) {
  n_marker <- min(2, length(marker))
  marker <- marker[seq_len(n_marker)]

  plot_tbl <- .plot_cyto_prep_plot_tbl(
    marker = marker,
    data = data,
    n_marker = n_marker,
    exc_min = exc_min
  )

  # axis labels
  if (!is.null(lab)) {
    marker <- lab[marker]
  }

  list(
    "n_marker" = n_marker,
    "marker" = marker,
    "plot_tbl" = plot_tbl
  )
}

.plot_cyto_prep_plot_tbl <- function(marker,
                                     data,
                                     n_marker,
                                     exc_min) {
  # plot_tbl
  plot_tbl <- data[, marker, drop = FALSE]
  colnames(plot_tbl) <- c("V1", "V2")[seq_len(n_marker)]
  .plot_cyto_prep_plot_tbl_exc_min(
    exc_min = exc_min,
    plot_tbl = plot_tbl,
    n_marker = n_marker
  )
}

.plot_cyto_prep_plot_tbl_exc_min <- function(exc_min,
                                             plot_tbl,
                                             n_marker) {
  if (!exc_min) {
    return(plot_tbl)
  }
  plot_tbl <- plot_tbl |>
    dplyr::filter(
      V1 > min(V1) # nolint
    )
  if (n_marker == 2) {
    plot_tbl <- plot_tbl |>
      dplyr::filter(
        V2 > min(V2) # nolint
      )
  }
  plot_tbl
}

.plot_cyto_plot <- function(n_marker,
                            plot_tbl,
                            marker,
                            font_size,
                            coord_equal,
                            limits_expand,
                            limits_equal,
                            geom_uni,
                            ...) {
  # base plot
  switch(n_marker,
    .plot_cyto_plot_uni(
      geom_uni = geom_uni,
      plot_tbl = plot_tbl,
      font_size = font_size,
      marker = marker,
      limits_expand = limits_expand,
      ...
    ),
    .plot_cyto_plot_biv(
      plot_tbl = plot_tbl,
      font_size = font_size,
      marker = marker,
      coord_equal = coord_equal,
      limits_expand = limits_expand,
      limits_equal = limits_equal,
      ...
    )
  )
}


.plot_cyto_plot_uni <- function(geom_uni,
                                plot_tbl,
                                font_size,
                                marker,
                                geom_uni_gg,
                                limits_expand,
                                ...) {
  geom_uni_gg <- .plot_cyto_plot_uni_geom(geom_uni, ...)

  p <- .plot_cyto_plot_uni_base(
    plot_tbl = plot_tbl,
    font_size = font_size,
    marker = marker,
    geom_uni_gg = geom_uni_gg
  )

  .plot_cyto_plot_uni_axes(
    p = p, limits_expand = limits_expand
  )
}

.plot_cyto_plot_uni_geom <- function(geom_uni, ...) {
  switch(geom_uni,
    "histogram" = do.call(
      ggplot2::geom_histogram,
      list(...)
    ),
    "density" = do.call(
      ggplot2::geom_density,
      list(...)
    ),
    stop("geom_unit value of ", geom_uni, " not recognised")
  )
}

.plot_cyto_plot_uni_base <- function(plot_tbl,
                                     font_size,
                                     marker,
                                     geom_uni_gg) {
  ggplot( # nolint
    plot_tbl, # nolint
    aes(x = V1) # nolint
  ) +
    cowplot::theme_cowplot(font_size) +
    cowplot::background_grid(major = "x") +
    labs(x = marker[1]) + # nolint
    geom_uni_gg
}


.plot_cyto_plot_uni_axes <- function(p, limits_expand) {
  .plot_cyto_plot_uni_axes_expand(
    p = p, limits_expand = limits_expand
  )
}

.plot_cyto_plot_uni_axes_expand <- function(p, limits_expand) {
  # return now if axis_limits fn not required
  if (is.null(limits_expand)) {
    return(p)
  }
  axis_limits(
    p = p, limits_expand = limits_expand
  )
}

.plot_cyto_plot_biv <- function(plot_tbl,
                                font_size,
                                marker,
                                coord_equal,
                                limits_expand,
                                limits_equal,
                                ...) {
  p <- .plot_cyto_plot_biv_base(
    plot_tbl = plot_tbl,
    font_size = font_size,
    marker = marker,
    ...
  )

  .plot_cyto_biv_axes(
    p = p,
    coord_equal = coord_equal,
    limits_expand = limits_expand,
    limits_equal = limits_equal
  )
}

.plot_cyto_plot_biv_base <- function(plot_tbl,
                                     font_size,
                                     marker,
                                     ...) {
  ggplot( # nolint
    plot_tbl, aes(x = V1, y = V2) # nolint
  ) +
    cowplot::theme_cowplot(font_size) +
    theme( # nolint
      plot.background = element_rect(fill = "white"), # nolint
      panel.background = element_rect(fill = "white") # nolint
    ) +
    geom_hex(...) + # nolint
    scale_fill_viridis_c( # nolint
      trans = "log10",
      name = "Count"
    ) +
    cowplot::background_grid(major = "xy") +
    labs(x = marker[1], y = marker[2]) # nolint
}

.plot_cyto_biv_axes <- function(p,
                                coord_equal,
                                limits_expand,
                                limits_equal) {
  if (coord_equal) p <- p + coord_equal()

  # return now if axis_limits fn not required
  if (is.null(limits_expand) && !limits_equal) {
    return(p)
  }

  .plot_cyto_biv_axes_expand(
    p = p,
    limits_expand = limits_expand,
    limits_equal = limits_equal
  )
}

.plot_cyto_biv_axes_expand <- function(p,
                                       limits_equal,
                                       limits_expand) {
  axis_limits(
    p = p, limits_expand = limits_expand, limits_equal = limits_equal
  )
}
