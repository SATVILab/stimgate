#' @title Manage axis limits
#'
#' @description Manage axis limits.
#'  Fix axis limits to be equal between x- and y-axes,
#' and/or expand axis coordinates.
#' The primary use of `axis_limits`
#' is forcing the x- and y-axes
#' to have the same limits "automatically"
#' (i.e. by inspecting the `ggplot` object,
#' thus not requiring the user to manually
#' calculate limits to pass to `ggplot2::expand_limits`).
#'
#' @param p object of class 'ggplot'. Limits are adjusted for this plot.
#' @param limits_expand list. If not \code{NULL},
#' then it is (effectively) passed onto \code{ggplot2::expand_limits} to
#' ensure that certain values are included in the plot (such as, for example, 0
#' if that is the minimum value possible but it may not be plotted). If not named, then # nolint
#' must consist of one numeric vector that will then force all values in the numeric value # nolint
#' to be included in the plot. If named, then must have names \code{x} and/or \code{y}, # nolint
#' with the elements again being numeric vectors that must be included in plot.
#' @param limits_equal logical. If \code{TRUE}, then the ranges on the x- and y-axes # nolint
#' must be equal. Effectively applied after expand_grid is applied. Default is \code{FALSE}. # nolint
#'
#' @export
#'

#'
#' @examples
#' data("cars", package = "datasets")
#' library(ggplot2)
#' p <- ggplot(cars, aes(speed, dist)) +
#'   geom_point()
#'
#' axis_limits(
#'   p,
#'   limits_equal = TRUE
#' )
#'
#' # both axes
#' axis_limits(
#'   p,
#'   limits_expand = list(200)
#' )
#' # x only
#' axis_limits(
#'   p,
#'   limits_expand = list(x = 75)
#' )
#' # y only
#' axis_limits(
#'   p,
#'   limits_expand = list(y = 200)
#' )
#' # lower and upper expansion
#' axis_limits(
#'   p,
#'   limits_expand = list(
#'     y = c(-50, 200),
#'     x = c(-10, 75)
#'   )
#' )
#'
#' # note that when fixing range and expanding, range is fixed
#' # after expansions are applied, so effectively the larger expansions apply to both.
#' # compare the following output to the previous output:
#' axis_limits(
#'   p,
#'   limits_expand = list(
#'     y = c(-50, 200),
#'     x = c(-10, 75)
#'   ),
#'   limits_equal = TRUE
#' )
axis_limits <- function(p,
                        limits_expand = NULL,
                        limits_equal = FALSE) {

  # initial check
  # ------------------------

  if (!is.logical(limits_equal)) {
    stop("limits_equal must be logical")
  }

  # do nothing
  if (is.null(limits_expand) && !limits_equal) {
    return(p)
  }

  # checks
  # ----------------------------

  if (!identical(
    class(p),
    c("gg", "ggplot")
  )) {
    stop("p must be of class c('gg', 'ggplot')")
  }

  if (!is.null(limits_expand)) {
    if (!is.list(limits_expand)) {
      stop("limits_expand must be a list (if not NULL)")
    }
    if (length(limits_expand) == 2 && is.null(names(limits_expand))) {
      stop("limits_expand must be named if of length 2")
    }
    if (length(limits_expand) > 2) {
      stop("limits_expand must have length 1 or 2 (if not NULL)")
    }
    if (!is.null(names(limits_expand))) {
      if (length(setdiff(names(limits_expand), c("x", "y"))) > 0) {
        stop("limits_expand must have names of 'x' and/or 'y' (if named)")
      }
    }
    class_input <- purrr::map_lgl(limits_expand, is.numeric) |> all()
    if (!class_input) {
      stop("input to limits_expand must be numeric (if limits_expand not NULL)")
    }
  }

  # prep
  # -------------------



  # ===================
  # adjustments
  # ===================

  # calc ranges in advance if needed
  # --------------------
  if (limits_equal) {
    plot_tbl <- p$data
    x_var <- as.character(rlang::get_expr(p$mapping$x))
    y_var <- as.character(rlang::get_expr(p$mapping$y))
    range_x <- range(plot_tbl[[x_var]])
    range_y <- range(plot_tbl[[y_var]])
    range <- c(
      min(range_x[1], range_y[1]),
      c(max(range_x[2], range_y[2]))
    )
  }

  # tidy limits_expand if provided
  # ------------------

  # ensure that limits_expand is named if
  # it's specified
  if (!is.null(limits_expand)) {
    if (is.null(names(limits_expand))) {
      limits_expand <- list(
        x = limits_expand[[1]],
        y = limits_expand[[1]]
      )
    }
    # ensure that limits_expand consists of
    # two sorted (not strictly) variables
    for (i in seq_along(limits_expand)) {
      limits_expand[[i]] <- c(
        min(limits_expand[[i]]),
        max(limits_expand[[i]])
      )
    }
  }

  # put limits_expand together with limits_equal,
  # if provided
  if (is.null(limits_expand)) {
    # we know now that limits_equal is true
    limits_expand <- list(
      x = range,
      y = range
    )
  } else {
    # limits_equal may or may not be true
    if (limits_equal) {
      limits_expand_all <- limits_expand |>
        unlist()
      lims <- c(
        min(range, limits_expand_all),
        max(range, limits_expand_all)
      )
      for (i in seq_along(limits_expand)) {
        limits_expand[[i]] <- lims
      }
      if (length(limits_expand) == 1) {
        nm <- setdiff(c("x", "y"), names(limits_expand))
        limits_expand <- limits_expand |>
          append(list(range) |> stats::setNames(nm))
      }
    }
  }

  limits_expand_arg <- purrr::map_chr(seq_along(limits_expand), function(i) {
    vals <- paste0(limits_expand[[i]], collapse = ", ")
    paste0(names(limits_expand)[i], " = c(", vals, ")")
  }) |>
    paste0(collapse = ", ")

  parse_text <- paste0("p <- p + ggplot2::expand_limits(", limits_expand_arg, ")")
  env <- environment()
  eval(parse(text = parse_text), envir = env)

  p
}
