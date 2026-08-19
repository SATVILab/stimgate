root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_cyt <- file.path(root_dir, "scripts", "r", "functionsForBenchmarking-Cyt.R")
script_trans <- file.path(root_dir, "scripts", "r", "sim-trans.R")

settings_tbl <- tibble::tibble(
  n_cell = 5000,
  background_relative_to_response = 0.2,
  prob_exact = TRUE,
  mixture_type = "gaussianOnly",
  cluster_perturbation_sd = 0,
  cov_ev_min = 1.5,
  cov_ev_max = 1.5
)

test_that("the experiment-based helper keeps QMD gn/gp semantics and stays distinct from the condition-based helper", {
  if (!file.exists(script_trans)) {
    stop("Expected analysis helper not found: ", script_trans)
  }

  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_cyt, local = env)
  source(script_trans, local = env)

  old_out <- env$sim_trans_univariate_one(
    transformation = "gaussian",
    n_cell = settings_tbl$n_cell[[1]],
    mean_pos = 1.5,
    prob_response = 0.01,
    background_relative_to_response = settings_tbl$background_relative_to_response[[1]],
    prob_exact = settings_tbl$prob_exact[[1]],
    mixture_type = settings_tbl$mixture_type[[1]],
    cluster_perturbation_sd = settings_tbl$cluster_perturbation_sd[[1]],
    cov_ev_min = settings_tbl$cov_ev_min[[1]],
    cov_ev_max = settings_tbl$cov_ev_max[[1]]
  )

  expect_true(all(c("negative", "response") %in% unique(old_out$response_class)))

  experiment_out <- env$sim_trans_univariate_experiment_one(
    transformation = "gaussian",
    mean_pos = 1.5,
    prob_response = 0.01,
    settings = settings_tbl
  )

  expect_setequal(unique(experiment_out$condition), c("unstimulated", "stimulated"))
  expect_true(all(experiment_out$response_class %in% c("gn", "gp")))
  expect_true(any(experiment_out$response_class == "gn"))
  expect_true(any(experiment_out$response_class == "gp"))
  expect_false(identical(
    env$sim_trans_univariate_one,
    env$sim_trans_univariate_experiment_one
  ))
})

test_that("make_density_tbl handles gamma fixed-range and non-gamma densities without browser debugging", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_cyt, local = env)
  source(script_trans, local = env)

  lines <- readLines(script_trans, warn = FALSE)
  expect_false(any(grepl("browser\\s*\\(", lines)))

  gamma_tbl <- tibble::tibble(
    mean_pos_setting = "low",
    mean_pos = 1,
    prob_response = 0.01,
    transformation = "gamma",
    condition = "stimulated",
    F1 = c(
      stats::runif(200, min = -0.25, max = 0.5),
      stats::runif(200, min = 1.2, max = 2.7)
    )
  )

  gamma_out <- env$make_density_tbl(gamma_tbl, gamma_range = c(-0.25, 2.7))
  expect_gt(nrow(gamma_out), 0)
  expect_true(all(is.finite(gamma_out$F1)))
  expect_true(all(is.finite(gamma_out$density)))

  gaussian_tbl <- dplyr::mutate(gamma_tbl, transformation = "gaussian")
  gaussian_out <- env$make_density_tbl(gaussian_tbl, gamma_range = c(-0.25, 2.7))
  expect_gt(nrow(gaussian_out), 0)
  expect_true(all(is.finite(gaussian_out$F1)))
  expect_true(all(is.finite(gaussian_out$density)))
})
