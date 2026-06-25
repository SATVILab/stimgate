#' Manage axis limits
#'
#' Manage axis limits.
#'  Fix axis limits to be equal between x- and y-axes,
#' and/or expand axis coordinates.
#' The primary use of `axisLimits`
#' is forcing the x- and y-axes
#' to have the same limits "automatically"
#' (i.e. by inspecting the `ggplot` object,
#' thus not requiring the user to manually
#' calculate limits to pass to `ggplot2::expand_limits`).
#'
#' @param p object of class 'ggplot' Limits are adjusted for this plot.
#' @param limitsExpand list or NULL If not NULL,
#' then it is (effectively) passed onto ggplot2::expand_limits to
#' ensure that certain values are included in the plot (such as, for
#' example, 0 if that is the minimum value possible but it may not be
#' plotted). If not named, then must consist of one numeric vector
#' that will then force all values in the numeric value to be included
#' in the plot. If named, then must have names x and/or y, with the
#' elements again being numeric vectors that must be included in plot.
#' Default is NULL.
#' @param limitsEqual logical If TRUE, then the ranges on the x- and
#' y-axes must be equal. Effectively applied after expand_grid is
#' applied. Default is FALSE.
#'
#' @return A ggplot object with adjusted axis limits.
#'
#' @export
#'
#' @examples
#' data("cars", package = "datasets")
#' library(ggplot2)
#' p <- ggplot(cars, aes(speed, dist)) +
#'   geom_point()
#'
#' axisLimits(
#'   p,
#'   limitsEqual = TRUE
#' )
#'
#' # both axes
#' axisLimits(
#'   p,
#'   limitsExpand = list(200)
#' )
#' # x only
#' axisLimits(
#'   p,
#'   limitsExpand = list(x = 75)
#' )
#' # y only
#' axisLimits(
#'   p,
#'   limitsExpand = list(y = 200)
#' )
#' # lower and upper expansion
#' axisLimits(
#'   p,
#'   limitsExpand = list(
#'     y = c(-50, 200),
#'     x = c(-10, 75)
#'   )
#' )
#'
#' # note that when fixing range and expanding, range is fixed
#' # after expansions are applied, so effectively the larger
#' # expansions apply to both.
#' # compare the following output to the previous output:
#' axisLimits(
#'   p,
#'   limitsExpand = list(
#'     y = c(-50, 200),
#'     x = c(-10, 75)
#'   ),
#'   limitsEqual = TRUE
#' )
axisLimits <- function(p, limitsExpand = NULL, limitsEqual = FALSE) {
  # initial check
  # ------------------------

  if (!is.logical(limitsEqual)) {
    stop("limitsEqual must be logical")
  }

  # do nothing
  if (is.null(limitsExpand) && !limitsEqual) {
    return(p)
  }

  # checks
  # ----------------------------

  if (!length(intersect(class(p), c("gg", "ggplot"))) == 2) {
    stop("p must be of class 'ggplot' and 'gg'")
  }

  if (!is.null(limitsExpand)) {
    if (!is.list(limitsExpand)) {
      stop("limitsExpand must be a list (if not NULL)")
    }
    if (length(limitsExpand) == 2 && is.null(names(limitsExpand))) {
      stop("limitsExpand must be named if of length 2")
    }
    if (length(limitsExpand) > 2) {
      stop("limitsExpand must have length 1 or 2 (if not NULL)")
    }
    if (!is.null(names(limitsExpand))) {
      if (length(setdiff(names(limitsExpand), c("x", "y"))) > 0) {
        stop("limitsExpand must have names of 'x' and/or 'y' (if named)")
      }
    }
    classInput <- purrr::map_lgl(limitsExpand, is.numeric) |> all()
    if (!classInput) {
      stop("input to limitsExpand must be numeric (if limitsExpand not NULL)")
    }
  }

  # prep
  # -------------------

  # ===================
  # adjustments
  # ===================

  # calc ranges in advance if needed
  # --------------------
  if (limitsEqual) {
    plotTbl <- p$data
    xVar <- as.character(rlang::get_expr(p$mapping$x))
    yVar <- as.character(rlang::get_expr(p$mapping$y))
    rangeX <- range(plotTbl[[xVar]])
    rangeY <- range(plotTbl[[yVar]])
    range <- c(
      min(rangeX[1], rangeY[1]),
      c(max(rangeX[2], rangeY[2]))
    )
  }

  # tidy limitsExpand if provided
  # ------------------

  # ensure that limitsExpand is named if
  # it's specified
  if (!is.null(limitsExpand)) {
    if (is.null(names(limitsExpand))) {
      limitsExpand <- list(
        x = limitsExpand[[1]],
        y = limitsExpand[[1]]
      )
    }
    # ensure that limitsExpand consists of
    # two sorted (not strictly) variables
    for (i in seq_along(limitsExpand)) {
      limitsExpand[[i]] <- c(
        min(limitsExpand[[i]]),
        max(limitsExpand[[i]])
      )
    }
  }

  # put limitsExpand together with limitsEqual,
  # if provided
  if (is.null(limitsExpand)) {
    # we know now that limitsEqual is true
    limitsExpand <- list(
      x = range,
      y = range
    )
  } else {
    # limitsEqual may or may not be true
    if (limitsEqual) {
      limitsExpandAll <- limitsExpand |>
        unlist()
      lims <- c(
        min(range, limitsExpandAll),
        max(range, limitsExpandAll)
      )
      for (i in seq_along(limitsExpand)) {
        limitsExpand[[i]] <- lims
      }
      if (length(limitsExpand) == 1) {
        nm <- setdiff(c("x", "y"), names(limitsExpand))
        limitsExpand <- limitsExpand |>
          append(list(range) |> stats::setNames(nm))
      }
    }
  }

  limitsExpandArg <- purrr::map_chr(seq_along(limitsExpand), function(i) {
    vals <- paste0(limitsExpand[[i]], collapse = ", ")
    paste0(names(limitsExpand)[i], " = c(", vals, ")")
  }) |>
    paste0(collapse = ", ")

  parseText <- paste0(
    "p <- p + ggplot2::expand_limits(",
    limitsExpandArg,
    ")"
  )
  env <- environment()
  eval(parse(text = parseText), envir = env)

  p
}
