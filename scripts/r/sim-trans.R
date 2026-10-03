sim_trans_univariate_experiment_one <- function(
    transformation,
    mean_pos,
    prob_response,
    settings = NULL) {
  if (is.null(settings)) {
    stop(
      "settings must be supplied explicitly to sim_trans_univariate_experiment_one()."
    )
  }

  n_cell <- settings$n_cell[[1]]
  background_relative_to_response <- settings$background_relative_to_response[[1]]
  prob_response_uns <- prob_response * background_relative_to_response

  transformation_fun <- .simMiscGetTrans(transformation)
  attr(transformation_fun, "sim_transformation") <- transformation

  out <- simcyto::simCytExperiment(
    nSample = 1L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nCellByCondition = c(n_cell, n_cell),
    transformationFunc = transformation_fun,
    mixtureType = settings$mixture_type[[1]],
    meanExprMat = matrix(
      c(0, mean_pos),
      byrow = TRUE,
      ncol = 1
    ),
    clusterLabelVec = c("gn", "gp"),
    probVecUns = c(
      1 - prob_response_uns,
      prob_response_uns
    ),
    probExact = settings$prob_exact[[1]],
    probResponseVecByStimCondition = list(
      c(-prob_response, prob_response)
    ),
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = settings$cluster_perturbation_sd[[1]],
    covEvMin = settings$cov_ev_min[[1]],
    covEvMax = settings$cov_ev_max[[1]]
  )

  flow_frame_list <- out$flowFrameList
  labels_list <- out$labelsList

  purrr::map_dfr(seq_along(flow_frame_list), function(ind) {
    tibble::tibble(
      transformation = transformation,
      mean_pos = mean_pos,
      condition = c("unstimulated", "stimulated")[[ind]],
      response_class = labels_list[[ind]],
      F1 = as.numeric(flowCore::exprs(flow_frame_list[[ind]])[, "F1"])
    )
  })
}

make_density_tbl <- function(.data, gamma_range, n = 2048) {
  .data |>
    dplyr::group_by(
      .data$mean_pos_setting,
      .data$mean_pos,
      .data$prob_response,
      .data$transformation,
      .data$condition
    ) |>
    dplyr::group_modify(function(.x, .y) {
      x <- .x$F1
      x <- x[is.finite(x)]

      if (length(x) < 2L || length(unique(x)) < 2L) {
        return(tibble::tibble(F1 = numeric(), density = numeric()))
      }

      bw <- switch(as.character(.y$transformation)[[1]],
        "gamma" = 0.025,
        0.15
      )

      if (.y$transformation[[1]] == "gamma") {
        dens <- stats::density(
          x,
          n = n,
          bw = bw,
          from = gamma_range[1],
          to = gamma_range[2]
        )
      } else {
        dens <- stats::density(x, n = n, bw = bw)
      }

      tibble::tibble(
        F1 = dens$x,
        density = dens$y
      )
    }) |>
    dplyr::ungroup()
}

# Set the RNG for one simulation row so that any single row can be
# regenerated with identical random numbers.
sim_trans_set_row_seed <- function(simulation_seed, row_index) {
  set.seed(
    as.integer(simulation_seed) + as.integer(row_index) - 1L,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
}

# Provenance stored next to the cached table: everything the table depends on.
sim_trans_cache_settings <- function(
    main_settings,
    mean_pos_settings_tbl,
    transformation_vec,
    simulation_seed) {
  list(
    main_settings = as.data.frame(main_settings),
    mean_pos_settings_tbl = as.data.frame(mean_pos_settings_tbl),
    transformation_vec = as.character(transformation_vec),
    simulation_seed = as.integer(simulation_seed)
  )
}

sim_trans_write_cache <- function(uni_tbl, settings, path) {
  .write_rds_atomic(list(settings = settings, uni_tbl = uni_tbl), path)
}

# Read the cached table, stopping if it was made with different settings.
sim_trans_read_cache <- function(path, settings) {
  if (!file.exists(path)) {
    stop("Cached uni_tbl.rds not found. Run simulations first.")
  }
  cached <- readRDS(path)
  if (
    !is.list(cached) ||
      !all(c("settings", "uni_tbl") %in% names(cached))
  ) {
    stop(
      "Cached uni_tbl.rds at ", path, " has no stored settings (legacy ",
      "cache). Rerun with run_simulations = TRUE."
    )
  }
  if (!isTRUE(all.equal(cached$settings, settings))) {
    stop(
      "Cached uni_tbl.rds at ", path, " was made with different settings ",
      "(main_settings, mean_pos settings, transformations or ",
      "simulation_seed) from the current ones. Rerun with ",
      "run_simulations = TRUE."
    )
  }
  cached$uni_tbl
}
