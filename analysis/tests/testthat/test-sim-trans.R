root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_trans <- file.path(root_dir, "scripts", "r", "sim-trans.R")

settings_tbl <- tibble::tibble(
  n_cell = 5000L,
  background_relative_to_response = 0.2,
  prob_exact = TRUE,
  mixture_type = "gaussianOnly",
  cluster_perturbation_sd = 0,
  cov_ev_min = 1.5,
  cov_ev_max = 1.5
)

test_that("transformation simulation helpers source cleanly without legacy functionsForBenchmarking-Cyt.R", {
  for (f in c(script_misc, script_trans)) {
    if (!file.exists(f)) stop("Expected analysis helper not found: ", f)
  }

  env <- new.env(parent = getNamespace("stimgate"))
  expect_no_error(source(script_misc, local = env))
  expect_no_error(source(script_trans, local = env))
  expect_false(exists("simCytExperiment", envir = env, inherits = FALSE))
  expect_false(exists("simCytCondition", envir = env, inherits = FALSE))
})

test_that("analysis/1-sim-trans.qmd does not source functionsForBenchmarking-Cyt.R", {
  qmd_path <- file.path(root_dir, "analysis", "1-sim-trans.qmd")
  expect_true(file.exists(qmd_path))

  lines <- readLines(qmd_path, warn = FALSE)
  expect_false(
    any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
    info = "analysis/1-sim-trans.qmd should not source functionsForBenchmarking-Cyt.R"
  )
})

test_that("sim_trans_univariate_experiment_one fixed-seed parity matches direct simcyto for gaussian, gamma and skew", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_trans, local = env)

  as_tbl_from_sim <- function(sim_obj, transformation, mean_pos) {
    purrr::map_dfr(seq_along(sim_obj$flowFrameList), function(ind) {
      tibble::tibble(
        transformation = transformation,
        mean_pos = mean_pos,
        condition = c("unstimulated", "stimulated")[[ind]],
        response_class = sim_obj$labelsList[[ind]],
        F1 = as.numeric(flowCore::exprs(sim_obj$flowFrameList[[ind]])[, "F1"])
      )
    })
  }

  run_case <- function(seed, transformation, mean_pos, prob_response) {
    prob_response_uns <- prob_response * settings_tbl$background_relative_to_response[[1]]
    captured_sim <- NULL
    orig_simcyto_experiment <- simcyto::simCytExperiment

    set.seed(seed)
    helper_tbl <- testthat::with_mocked_bindings(
      simCytExperiment = function(...) {
        out <- orig_simcyto_experiment(...)
        captured_sim <<- out
        out
      },
      .package = "simcyto",
      env$sim_trans_univariate_experiment_one(
        transformation = transformation,
        mean_pos = mean_pos,
        prob_response = prob_response,
        settings = settings_tbl
      )
    )

    expect_type(captured_sim, "list")
    expect_true(all(c("flowFrameList", "labelsList") %in% names(captured_sim)))

    set.seed(seed)
    direct_sim <- simcyto::simCytExperiment(
      nSample = 1L,
      nMarker = 1L,
      nCondition = 2L,
      nCluster = 2L,
      nCellByCondition = c(settings_tbl$n_cell[[1]], settings_tbl$n_cell[[1]]),
      transformationFunc = env$.simMiscGetTrans(transformation),
      mixtureType = settings_tbl$mixture_type[[1]],
      meanExprMat = matrix(c(0, mean_pos), byrow = TRUE, ncol = 1),
      clusterLabelVec = c("gn", "gp"),
      probVecUns = c(1 - prob_response_uns, prob_response_uns),
      probExact = settings_tbl$prob_exact[[1]],
      probResponseVecByStimCondition = list(c(-prob_response, prob_response)),
      samplePerturbationSd = 0,
      conditionPerturbationSd = 0,
      clusterPerturbationSd = settings_tbl$cluster_perturbation_sd[[1]],
      covEvMin = settings_tbl$cov_ev_min[[1]],
      covEvMax = settings_tbl$cov_ev_max[[1]]
    )

    direct_tbl <- as_tbl_from_sim(direct_sim, transformation, mean_pos)

    expect_equal(helper_tbl, direct_tbl, tolerance = 1e-12)
    expect_equal(captured_sim$labelsList, direct_sim$labelsList)

    count_helper <- helper_tbl |>
      dplyr::count(.data$condition, .data$response_class, name = "n_cell") |>
      dplyr::arrange(.data$condition, .data$response_class)
    count_direct <- direct_tbl |>
      dplyr::count(.data$condition, .data$response_class, name = "n_cell") |>
      dplyr::arrange(.data$condition, .data$response_class)
    expect_equal(count_helper, count_direct)
    expect_equal(sum(count_helper$n_cell), 2L * settings_tbl$n_cell[[1]])

    moments_helper <- helper_tbl |>
      dplyr::group_by(.data$condition, .data$response_class) |>
      dplyr::summarise(
        mean_F1 = mean(.data$F1),
        sd_F1 = stats::sd(.data$F1),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$condition, .data$response_class)
    moments_direct <- direct_tbl |>
      dplyr::group_by(.data$condition, .data$response_class) |>
      dplyr::summarise(
        mean_F1 = mean(.data$F1),
        sd_F1 = stats::sd(.data$F1),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$condition, .data$response_class)
    expect_equal(moments_helper, moments_direct, tolerance = 1e-12)
    expect_true(
      moments_helper$mean_F1[moments_helper$condition == "stimulated" &
        moments_helper$response_class == "gp"] >
        moments_helper$mean_F1[moments_helper$condition == "unstimulated" &
          moments_helper$response_class == "gn"]
    )

    dens_input_helper <- helper_tbl |>
      dplyr::mutate(
        mean_pos_setting = "low",
        prob_response = prob_response
      )
    dens_input_direct <- direct_tbl |>
      dplyr::mutate(
        mean_pos_setting = "low",
        prob_response = prob_response
      )

    gamma_range <- range(dens_input_helper$F1)
    if (identical(transformation, "gamma")) {
      gamma_range <- c(
        gamma_range[1] - 0.05 * abs(diff(gamma_range)),
        gamma_range[2]
      )
    }

    dens_helper <- env$make_density_tbl(dens_input_helper, gamma_range = gamma_range) |>
      dplyr::arrange(.data$condition, .data$F1)
    dens_direct <- env$make_density_tbl(dens_input_direct, gamma_range = gamma_range) |>
      dplyr::arrange(.data$condition, .data$F1)
    expect_equal(dens_helper, dens_direct, tolerance = 1e-12)
    expect_true(all(is.finite(dens_helper$density)))
  }

  run_case(seed = 101L, transformation = "gaussian", mean_pos = 4.5, prob_response = 0.002)
  run_case(seed = 202L, transformation = "gamma", mean_pos = 4, prob_response = 0.002)
  run_case(seed = 303L, transformation = "skew", mean_pos = 6, prob_response = 0.002)
})

test_that("make_density_tbl handles gamma fixed-range and non-gamma densities without browser debugging", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
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
