.simMiscGetTrans <- function(transformation) {
  switch(
    transformation,
    "gamma" = calc_gamma,
    "gaussian" = calc_gaussian,
    "skew" = calc_skew,
    stop("transformation not recognised")
  )
}

calc_gamma <- function(x) {
  gamma(1 + abs(x / 4))
}
calc_gaussian <- function(x) {
  x
}
calc_skew <- function(x, epsilon = 2, delta = 1) {
  cluster_mean <- mean(x)
  cluster_sd <- sd(x)
  weight <- 1 / (1 + exp(1.5 * (cluster_mean - 3.5)))
  epsilon <- epsilon * weight
  out <- sinh(epsilon + delta * asinh(x))
  gamma_divisor <- rgamma(length(x), shape = 5, rate = 5)
  out_fat <- out / sqrt(gamma_divisor)
  (out - mean(out)) / sd(out) * sd(x) + cluster_mean
}

