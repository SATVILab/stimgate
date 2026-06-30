.simMiscTagTrans <- function(fun, name) {
  attr(fun, "sim_transformation") <- name
  fun
}

.simMiscGetTrans <- function(transformation) {
  switch(
    transformation,
    "gamma" = .simMiscTagTrans(calc_gamma, "gamma"),
    "gamma_fixed_mean_and_spread" = .simMiscTagTrans(
      calc_gamma_fixed_mean_and_spread,
      "gamma_fixed_mean_and_spread"
    ),
    "gaussian" = .simMiscTagTrans(calc_gaussian, "gaussian"),
    "skew" = .simMiscTagTrans(calc_skew, "skew"),
    stop("transformation not recognised")
  )
}

calc_gamma <- function(x) {
  gamma(1 + abs(x / 4))
}

calc_gamma_fixed_mean_and_spread <- function(x) {
  cluster_mean <- mean(x)
  cluster_sd <- sd(x)
  out <- calc_gamma(x)
  (out - mean(out)) / sd(out) * cluster_sd + cluster_mean
}

calc_gaussian <- function(x) {
  x
}

calc_skew <- function(x, epsilon = 0.5, delta = 1) {
  cluster_mean <- mean(x)
  weight <- 1 / (1 + exp(1.5 * (cluster_mean - 3.5)))
  epsilon <- epsilon * weight
  out <- sinh(epsilon + delta * asinh(x))
  gamma_divisor <- rgamma(length(x), shape = 5, rate = 5)
  out / sqrt(gamma_divisor)
}


set.seed(3)
rnorm_vec <- rnorm(1e3)
sample_skew_corrected(rnorm_vec)
sample_skew_corrected(calc_skew(rnorm_vec))
