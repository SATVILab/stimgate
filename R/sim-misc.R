.simMiscGetTrans <- function(transformation) {
  switch(
    transformation,
    "gamma" = calc_gamma,
    "gaussian" = calc_gaussian,
    "skew" = calc_skew,
    stop("transformation not recognised")
  )
}