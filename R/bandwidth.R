# Shared ordinary scalar bandwidth selection helpers.

#' Calculate an ordinary scalar bandwidth
#' @keywords internal
.bwCalcOne <- function(x, bwMtd, bwAdj = 1, bwNcellMin = NULL, bwNcellMax = NULL) {
  if (!is.character(bwMtd) || length(bwMtd) != 1L ||
      is.na(bwMtd) || !bwMtd %in% .verifyBwMtds) {
    stop("Unsupported `bwMtd`; use one of: ", paste(.verifyBwMtds, collapse = ", "))
  }
  xBw <- .bwCalcOneSampleOrdinary(x, bwNcellMin, bwNcellMax)
  bwOut <- .bwCalcOneBase(xBw, bwMtd)
  if (!is.finite(bwOut) || bwOut <= 0) {
    return(NA_real_)
  }
  as.numeric(bwOut)[1] * bwAdj
}

#' @keywords internal
.bwCalcOneSampleOrdinary <- function(
  x,
  bwNcellMin = NULL,
  bwNcellMax = NULL
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(x)
  }

  bwNcellMin <- .bwAsSafeSampleN(
    bwNcellMin,
    default = NULL,
    lower = 2L
  )

  bwNcellMax <- .bwAsSafeSampleN(
    bwNcellMax,
    default = NULL,
    lower = 2L
  )

  if (!is.null(bwNcellMin) && length(x) < bwNcellMin) {
    sdX <- .bwRobustSd(x)
    x <- sample(x, replace = TRUE, size = bwNcellMin) +
      stats::rnorm(bwNcellMin, mean = 0, sd = sdX / 10)
  }

  if (!is.null(bwNcellMax) && length(x) > bwNcellMax) {
    x <- sample(x, size = bwNcellMax, replace = FALSE)
  }

  x
}

#' @keywords internal
.bwCalcOneBase <- function(x, bwMtd) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(NA_real_)
  }

  bwCalc <- switch(bwMtd,
    "nrd0" = try(stats::bw.nrd0(x), silent = TRUE),
    "sj" = try(stats::bw.SJ(x), silent = TRUE),
    {
      derivOrder <- suppressWarnings(
        as.numeric(gsub("^hpi", "", bwMtd))
      )

      if (!is.finite(derivOrder)) {
        return(NA_real_)
      }

      try(
        suppressWarnings(
          ks::hpi(x, deriv.order = derivOrder)
        ),
        silent = TRUE
      )
    }
  )

  if (
    inherits(bwCalc, "try-error") ||
      !is.finite(bwCalc) ||
      bwCalc <= 0
  ) {
    bwCalc <- try(stats::bw.nrd0(x), silent = TRUE)
  }

  if (
    inherits(bwCalc, "try-error") ||
      !is.finite(bwCalc) ||
      bwCalc <= 0
  ) {
    return(NA_real_)
  }

  as.numeric(bwCalc)[1]
}

#' @keywords internal
.bwRobustSd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(.Machine$double.eps)
  }

  iqrX <- diff(stats::quantile(x, c(0.75, 0.25), na.rm = TRUE))
  sdX <- abs(iqrX) / 1.5

  if (!is.finite(sdX) || sdX <= 0) {
    sdX <- stats::sd(x, na.rm = TRUE)
  }

  if (!is.finite(sdX) || sdX <= 0) {
    sdX <- .Machine$double.eps
  }

  sdX
}

#' @keywords internal
.bwAsSafeSampleN <- function(
  x,
  default = NULL,
  lower = 0L,
  upper = .Machine$integer.max
) {
  if (is.null(x) || length(x) == 0L) {
    return(default)
  }

  x <- suppressWarnings(as.numeric(x)[1])

  if (!is.finite(x)) {
    return(default)
  }

  x <- floor(x)
  x <- max(as.numeric(lower), min(x, as.numeric(upper)))

  as.integer(x)
}
