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

.simMiscGetMeanPosTbl <- function() {
  tibble::tribble(
    ~transformation , ~mean_pos_setting , ~mean_pos ,
    "gaussian"      , "low"             , 4.5       ,
    "gaussian"      , "high"            , 8         ,
    "skew"          , "low"             , 6         ,
    "skew"          , "high"            , 8.5       ,
    "gamma"         , "low"             , 4         ,
    "gamma"         , "high"            , 7
  ) |>
    dplyr::mutate(
      mean_pos_setting = factor(
        .data$mean_pos_setting,
        levels = c("low", "high")
      )
    )
}

.simMiscGetStimColVec <- function() {
  uns_nm_vec <- c("uns", "unstim", "unstimulated")
  stim_nm_vec <- c("stim", "stimulated")
  c(
    rep("#0072B2", length(uns_nm_vec)), # blue
    rep("#D55E00", length(stim_nm_vec)) # vermillion
  ) |>
    stats::setNames(c(uns_nm_vec, stim_nm_vec))
}

.simGetCores <- function(n_cores = NULL) {
  if (!is.null(n_cores)) {
    return(n_cores)
  }
  n_tasks <- Sys.getenv("SLURM_NTASKS")
  if (!nzchar(n_tasks)) {
    return(future::availableCores() - 1L)
  }
  n_tasks <- as.integer(n_tasks)
  if (is.na(n_tasks)) {
    return(future::availableCores() - 1L)
  }
  job_name <- Sys.getenv("SLURM_JOB_NAME") |>
    as.character()
  if (!grepl("VScode$", job_name)) {
    n_tasks
  } else {
    n_tasks - 1L
  }
}

.simMiscGetTransPretty <- function() {
  c("gaussian" = "Gaussian", "gamma" = "Gamma", "skew" = "Skew")
}

.simMiscGetMeanPosPretty <- function() {
  c("low" = "Low", "high" = "High")
}
