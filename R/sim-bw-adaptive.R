.bwDiagGetTrans <- function(transformation) {
  if (exists(".simMiscGetTrans", mode = "function")) {
    return(.simMiscGetTrans(transformation))
  }

  switch(
    transformation,
    "gamma" = calc_gamma,
    "gamma_fixed_mean_and_spread" = calc_gamma_fixed_mean_and_spread,
    "gaussian" = calc_gaussian,
    "skew" = calc_skew,
    stop("Transformation not recognised: ", transformation)
  )
}

.bwDiagExcMin <- function(x, exc_min = TRUE) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  n_init <- length(x)

  if (!isTRUE(exc_min) || n_init == 0L) {
    return(list(x = x, prob_g_min = 1))
  }

  x_out <- x[x > min(x, na.rm = TRUE)]
  list(
    x = x_out,
    prob_g_min = length(x_out) / n_init
  )
}

.bwDiagCapForCpUnsLoc <- function(x_stim, x_uns, cap_stim_range = TRUE) {
  x_stim <- suppressWarnings(as.numeric(x_stim))
  x_uns <- suppressWarnings(as.numeric(x_uns))
  x_stim <- x_stim[is.finite(x_stim)]
  x_uns <- x_uns[is.finite(x_uns)]

  if (!isTRUE(cap_stim_range) || length(x_stim) < 2L) {
    return(list(
      x_stim = x_stim,
      x_uns = x_uns,
      max_dens_x = NA_real_
    ))
  }

  range_stim <- range(x_stim, na.rm = TRUE)
  max_dens_x <- max(x_stim, na.rm = TRUE) - 0.05 * diff(range_stim)

  list(
    x_stim = pmin(x_stim, max_dens_x),
    x_uns = pmin(x_uns, max_dens_x),
    max_dens_x = max_dens_x
  )
}

.bwDiagSimulateOne <- function(
  n_sample,
  n_marker,
  n_condition,
  n_cluster,
  n_cell_stim,
  ncell_uns_relative_to_stim,
  prob_response,
  background_relative_to_response,
  mean_pos,
  transformation,
  prob_exact,
  bias_uns,
  exc_min,
  cap_stim_range,
  cov_ev_min,
  cov_ev_max
) {
  n_cell_uns <- round(n_cell_stim * ncell_uns_relative_to_stim)
  n_cell_by_condition <- c(n_cell_uns, n_cell_stim)
  transformation_func <- .bwDiagGetTrans(transformation)

  mean_expr_mat <- matrix(
    c(0, mean_pos),
    byrow = TRUE,
    ncol = 1
  )

  prob_response_uns <- prob_response * background_relative_to_response
  prob_vec_uns <- c(1 - prob_response_uns, prob_response_uns)
  prob_response_vec_by_stim_condition <- list(c(-prob_response, prob_response))

  out_list_experiment <- simCytExperiment(
    nSample = n_sample,
    nMarker = n_marker,
    nCondition = n_condition,
    nCluster = n_cluster,
    nCellByCondition = n_cell_by_condition,
    transformationFunc = transformation_func,
    mixtureType = "gaussianOnly",
    meanExprMat = mean_expr_mat,
    clusterLabelVec = c("gn", "gp"),
    probVecUns = prob_vec_uns,
    probExact = prob_exact,
    probResponseVecByStimCondition = prob_response_vec_by_stim_condition,
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    covEvMin = cov_ev_min,
    covEvMax = cov_ev_max
  )

  flow_frame_list <- out_list_experiment[["flowFrameList"]]
  ind_uns <- 1L
  ind_stim <- 2L

  x_uns_raw <- as.numeric(flowCore::exprs(flow_frame_list[[ind_uns]])[, "F1"])
  x_stim_raw <- as.numeric(flowCore::exprs(flow_frame_list[[ind_stim]])[, "F1"])

  # cp_uns_loc applies the unstim bias before density estimation.
  x_uns_raw <- x_uns_raw + (bias_uns %||% 0)

  stim_min_obj <- .bwDiagExcMin(x_stim_raw, exc_min = exc_min)
  uns_min_obj <- .bwDiagExcMin(x_uns_raw, exc_min = exc_min)

  cap_obj <- .bwDiagCapForCpUnsLoc(
    x_stim = stim_min_obj$x,
    x_uns = uns_min_obj$x,
    cap_stim_range = cap_stim_range
  )

  list(
    x_stim = cap_obj$x_stim,
    x_uns = cap_obj$x_uns,
    prob_g_min_stim = stim_min_obj$prob_g_min,
    prob_g_min_uns = uns_min_obj$prob_g_min,
    max_dens_x = cap_obj$max_dens_x,
    labels_list = out_list_experiment[["labelsList"]]
  )
}

.bwDiagClipBw <- function(
  bw,
  bw_min = -Inf,
  bw_max = Inf,
  bw_fallback = NA_real_
) {
  bw <- suppressWarnings(as.numeric(bw)[1])
  if (!is.finite(bw) || bw <= 0) {
    bw <- suppressWarnings(as.numeric(bw_fallback)[1])
  }
  if (!is.finite(bw) || bw <= 0) {
    return(NA_real_)
  }
  if (is.finite(bw_min)) {
    bw <- max(bw_min, bw)
  }
  if (is.finite(bw_max)) {
    bw <- min(bw_max, bw)
  }
  bw
}

.bwDiagBwOne <- function(
  x,
  bw_mtd,
  bw_adaptive = FALSE,
  bw_adj = 1,
  bw_ncell_min = NULL,
  bw_ncell_max = NULL,
  bw_min = -Inf,
  bw_max = Inf,
  bw_fallback = NA_real_,
  norm_peak_frac = 0.1,
  norm_peak_min_rel = 0.75,
  norm_extra_frac = 0.2,
  norm_extra_max = Inf,
  norm_extra_jitter_frac = 0.25,
  norm_density_n = 512L,
  norm_excess_bw_mtd = "hpi3",
  norm_excess_ncell = 10000L,
  norm_adaptive_ncell = 2500L,
  norm_mtd = "moments"
) {
  use_adaptive <- isTRUE(bw_adaptive) && grepl("Norm$", bw_mtd)

  bw_args <- list(
    x = x,
    bwMtd = bw_mtd,
    bwAdj = bw_adj,
    bwNcellMin = bw_ncell_min,
    bwNcellMax = bw_ncell_max,
    normPeakFrac = norm_peak_frac,
    normPeakMinRel = norm_peak_min_rel,
    normExtraFrac = norm_extra_frac,
    normExtraMax = norm_extra_max,
    normExtraJitterFrac = norm_extra_jitter_frac,
    normDensityN = norm_density_n,
    normExcessBwMtd = norm_excess_bw_mtd,
    normExcessNcell = norm_excess_ncell,
    normAdaptiveNcell = norm_adaptive_ncell,
    normMtd = norm_mtd,
    adaptive = use_adaptive
  )

  # Keeps the notebook usable if it is accidentally run against an older helper:
  # unsupported arguments are dropped before calling .bwCalcOne().
  bw_args <- bw_args[intersect(names(bw_args), names(formals(.bwCalcOne)))]

  out <- tryCatch(
    do.call(.bwCalcOne, bw_args),
    error = function(e) {
      structure(NA_real_, error = e$message)
    }
  )

  if (is.list(out) && all(c("bin", "bw") %in% names(out))) {
    out$bw <- pmax(suppressWarnings(as.numeric(out$bw)), .Machine$double.eps)
    return(out)
  }

  .bwDiagClipBw(
    bw = out,
    bw_min = bw_min,
    bw_max = bw_max,
    bw_fallback = bw_fallback
  )
}

.bwDiagGrid <- function(x, n = 512L, pad_frac = 0.15) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  n <- max(16L, as.integer(n))

  range_vec <- range(x, na.rm = TRUE)
  range_width <- diff(range_vec)
  if (!is.finite(range_width) || range_width <= 0) {
    range_width <- stats::sd(x, na.rm = TRUE)
  }
  if (!is.finite(range_width) || range_width <= 0) {
    range_width <- 1
  }

  pad_frac <- suppressWarnings(as.numeric(pad_frac)[1])
  if (!is.finite(pad_frac) || pad_frac < 0) {
    pad_frac <- 0.15
  }

  pad <- pad_frac * range_width
  seq(range_vec[[1]] - pad, range_vec[[2]] + pad, length.out = n)
}

.bwDiagBwOnGrid <- function(bw_obj, grid, fallback = NA_real_) {
  grid <- suppressWarnings(as.numeric(grid))

  if (is.list(bw_obj) && all(c("bin", "bw") %in% names(bw_obj))) {
    return(
      stats::approx(
        x = suppressWarnings(as.numeric(bw_obj$bin)),
        y = suppressWarnings(as.numeric(bw_obj$bw)),
        xout = grid,
        rule = 2
      )$y
    )
  }

  bw_scalar <- .bwDiagClipBw(
    bw = bw_obj,
    bw_fallback = fallback
  )
  rep(bw_scalar, length(grid))
}

.bwDiagRepairBwGrid <- function(bw_grid) {
  bw_grid <- suppressWarnings(as.numeric(bw_grid))
  good <- is.finite(bw_grid) & bw_grid > 0
  if (!any(good)) {
    return(rep(NA_real_, length(bw_grid)))
  }
  bw_grid[!good] <- stats::median(bw_grid[good], na.rm = TRUE)
  pmax(bw_grid, .Machine$double.eps)
}

.bwDiagTrapz <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L || length(x) != length(y)) {
    return(NA_real_)
  }
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

.bwDiagDensityOnGrid <- function(
  x,
  grid,
  bw_grid,
  normalise = TRUE,
  prob_g_min = NULL
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  grid <- suppressWarnings(as.numeric(grid))
  bw_grid <- suppressWarnings(as.numeric(bw_grid))

  ok <- is.finite(grid) & is.finite(bw_grid) & bw_grid > 0
  grid <- grid[ok]
  bw_grid <- bw_grid[ok]

  if (
    length(x) < 2L ||
      length(unique(x)) < 2L ||
      length(grid) < 2L ||
      length(grid) != length(bw_grid)
  ) {
    return(NULL)
  }

  y <- vapply(
    seq_along(grid),
    function(i) {
      mean(stats::dnorm(grid[[i]], mean = x, sd = bw_grid[[i]]))
    },
    numeric(1)
  )
  y <- pmax(y, 0)

  area_before <- .bwDiagTrapz(grid, y)
  if (isTRUE(normalise) && is.finite(area_before) && area_before > 0) {
    y <- y / area_before
  }

  prob_g_min <- suppressWarnings(as.numeric(prob_g_min)[1])
  if (is.finite(prob_g_min)) {
    prob_g_min <- max(0, min(1, prob_g_min))
    y <- y * prob_g_min
  }

  tibble::tibble(
    x = grid,
    dens = y,
    bw = bw_grid,
    area_before_normalise = area_before,
    area_after = .bwDiagTrapz(grid, y)
  )
}

.bwDiagComponentRow <- function(bw_obj, arm, bw_mtd, bw_mode, n_cell) {
  if (is.list(bw_obj) && all(c("bin", "bw") %in% names(bw_obj))) {
    return(tibble::tibble(
      arm = arm,
      bw_mtd = bw_mtd,
      bw_mode = bw_mode,
      n_cell = n_cell,
      bw_scalar = NA_real_,
      bw_core = suppressWarnings(as.numeric(bw_obj$bwCore)[1]),
      bw_extra = suppressWarnings(as.numeric(bw_obj$bwExtra)[1]),
      n_adaptive = suppressWarnings(as.integer(bw_obj$nAdaptive)[1]),
      threshold_x = suppressWarnings(as.numeric(bw_obj$coreObj$thresholdX)[1]),
      peak_x = suppressWarnings(as.numeric(bw_obj$coreObj$peakX)[1]),
      adaptive_returned = TRUE
    ))
  }

  tibble::tibble(
    arm = arm,
    bw_mtd = bw_mtd,
    bw_mode = bw_mode,
    n_cell = n_cell,
    bw_scalar = suppressWarnings(as.numeric(bw_obj)[1]),
    bw_core = NA_real_,
    bw_extra = NA_real_,
    n_adaptive = NA_integer_,
    threshold_x = NA_real_,
    peak_x = NA_real_,
    adaptive_returned = FALSE
  )
}

.bwDiagRunOne <- function(
  x_stim,
  x_uns,
  prob_g_min_stim,
  prob_g_min_uns,
  bw_mtd,
  bw_adaptive,
  bw_adj,
  bw_ncell_min,
  bw_ncell_max,
  bw_min,
  bw_max,
  bw_fallback,
  norm_peak_frac,
  norm_peak_min_rel,
  norm_extra_frac,
  norm_extra_max,
  norm_extra_jitter_frac,
  norm_density_n,
  norm_excess_bw_mtd,
  norm_excess_ncell,
  norm_adaptive_ncell,
  norm_mtd,
  bw_adaptive_density_n,
  bw_adaptive_pad_frac
) {
  bw_mode <- if (isTRUE(bw_adaptive) && grepl("Norm$", bw_mtd)) {
    "adaptive"
  } else {
    "constant"
  }

  grid <- .bwDiagGrid(
    x = c(x_stim, x_uns),
    n = bw_adaptive_density_n,
    pad_frac = bw_adaptive_pad_frac
  )

  bw_stim_obj <- .bwDiagBwOne(
    x = x_stim,
    bw_mtd = bw_mtd,
    bw_adaptive = identical(bw_mode, "adaptive"),
    bw_adj = bw_adj,
    bw_ncell_min = bw_ncell_min,
    bw_ncell_max = bw_ncell_max,
    bw_min = bw_min,
    bw_max = bw_max,
    bw_fallback = bw_fallback,
    norm_peak_frac = norm_peak_frac,
    norm_peak_min_rel = norm_peak_min_rel,
    norm_extra_frac = norm_extra_frac,
    norm_extra_max = norm_extra_max,
    norm_extra_jitter_frac = norm_extra_jitter_frac,
    norm_density_n = norm_density_n,
    norm_excess_bw_mtd = norm_excess_bw_mtd,
    norm_excess_ncell = norm_excess_ncell,
    norm_adaptive_ncell = norm_adaptive_ncell,
    norm_mtd = norm_mtd
  )

  bw_uns_obj <- .bwDiagBwOne(
    x = x_uns,
    bw_mtd = bw_mtd,
    bw_adaptive = identical(bw_mode, "adaptive"),
    bw_adj = bw_adj,
    bw_ncell_min = bw_ncell_min,
    bw_ncell_max = bw_ncell_max,
    bw_min = bw_min,
    bw_max = bw_max,
    bw_fallback = bw_fallback,
    norm_peak_frac = norm_peak_frac,
    norm_peak_min_rel = norm_peak_min_rel,
    norm_extra_frac = norm_extra_frac,
    norm_extra_max = norm_extra_max,
    norm_extra_jitter_frac = norm_extra_jitter_frac,
    norm_density_n = norm_density_n,
    norm_excess_bw_mtd = norm_excess_bw_mtd,
    norm_excess_ncell = norm_excess_ncell,
    norm_adaptive_ncell = norm_adaptive_ncell,
    norm_mtd = norm_mtd
  )

  component_tbl <- dplyr::bind_rows(
    .bwDiagComponentRow(bw_stim_obj, "stim", bw_mtd, bw_mode, length(x_stim)),
    .bwDiagComponentRow(bw_uns_obj, "unstim", bw_mtd, bw_mode, length(x_uns))
  )

  bw_stim_grid <- .bwDiagRepairBwGrid(.bwDiagBwOnGrid(
    bw_stim_obj,
    grid,
    bw_fallback
  ))
  bw_uns_grid <- .bwDiagRepairBwGrid(.bwDiagBwOnGrid(
    bw_uns_obj,
    grid,
    bw_fallback
  ))

  if (identical(bw_mode, "constant")) {
    bw_shared_scalar <- min(
      component_tbl$bw_scalar[is.finite(component_tbl$bw_scalar)],
      na.rm = TRUE
    )
    if (!is.finite(bw_shared_scalar)) {
      bw_shared_scalar <- bw_fallback
    }
    bw_shared_grid <- rep(bw_shared_scalar, length(grid))
    dens_stim_weight <- rep(NA_real_, length(grid))
    dens_uns_weight <- rep(NA_real_, length(grid))
  } else {
    dens_stim_pre <- .bwDiagDensityOnGrid(
      x = x_stim,
      grid = grid,
      bw_grid = bw_stim_grid,
      normalise = TRUE,
      prob_g_min = prob_g_min_stim
    )
    dens_uns_pre <- .bwDiagDensityOnGrid(
      x = x_uns,
      grid = grid,
      bw_grid = bw_uns_grid,
      normalise = TRUE,
      prob_g_min = prob_g_min_uns
    )

    denom <- dens_stim_pre$dens + dens_uns_pre$dens
    bw_shared_grid <- ifelse(
      is.finite(denom) & denom > 0,
      (dens_stim_pre$dens * bw_stim_grid + dens_uns_pre$dens * bw_uns_grid) /
        denom,
      rowMeans(cbind(bw_stim_grid, bw_uns_grid), na.rm = TRUE)
    )
    bw_shared_grid <- .bwDiagRepairBwGrid(bw_shared_grid)
    dens_stim_weight <- dens_stim_pre$dens
    dens_uns_weight <- dens_uns_pre$dens
  }

  dens_stim <- .bwDiagDensityOnGrid(
    x = x_stim,
    grid = grid,
    bw_grid = bw_shared_grid,
    normalise = TRUE,
    prob_g_min = prob_g_min_stim
  ) |>
    dplyr::mutate(arm = "stim", density_stage = "final_shared")

  dens_uns <- .bwDiagDensityOnGrid(
    x = x_uns,
    grid = grid,
    bw_grid = bw_shared_grid,
    normalise = TRUE,
    prob_g_min = prob_g_min_uns
  ) |>
    dplyr::mutate(arm = "unstim", density_stage = "final_shared")

  bw_grid_tbl <- tibble::tibble(
    x = grid,
    bw_stim_own = bw_stim_grid,
    bw_uns_own = bw_uns_grid,
    bw_shared = bw_shared_grid,
    dens_stim_weight = dens_stim_weight,
    dens_uns_weight = dens_uns_weight
  ) |>
    tidyr::pivot_longer(
      cols = c("bw_stim_own", "bw_uns_own", "bw_shared"),
      names_to = "bw_curve",
      values_to = "bw"
    ) |>
    dplyr::mutate(
      bw_mtd = bw_mtd,
      bw_mode = bw_mode,
      config = paste0(bw_mtd, " / ", bw_mode)
    )

  density_tbl <- dplyr::bind_rows(dens_stim, dens_uns) |>
    dplyr::mutate(
      bw_mtd = bw_mtd,
      bw_mode = bw_mode,
      config = paste0(bw_mtd, " / ", bw_mode)
    )

  component_tbl <- component_tbl |>
    dplyr::mutate(config = paste0(bw_mtd, " / ", bw_mode))

  list(
    component = component_tbl,
    bw_grid = bw_grid_tbl,
    density = density_tbl
  )
}

.scalar_or_null <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NULL)
  }
  x <- x[[1]]
  if (length(x) == 0L || is.na(x)) {
    return(NULL)
  }
  x
}

.clean_file_token <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

.first_existing <- function(.data, cols, default = NA_real_) {
  col <- cols[cols %in% names(.data)][1]
  if (length(col) == 0L || is.na(col)) {
    return(rep(default, nrow(.data)))
  }
  .data[[col]]
}

.add_missing_cols <- function(.data, cols) {
  for (nm in names(cols)) {
    if (!nm %in% names(.data)) {
      .data[[nm]] <- cols[[nm]]
    }
  }
  .data
}

.truth_from_labels <- function(
  labelsList,
  nSample,
  nCondition,
  nCellStim,
  nCellUns
) {
  purrr::map_df(seq_len(nSample), function(sample_curr) {
    ind_uns <- (sample_curr - 1L) * nCondition + 1L
    ind_stim <- seq.int(ind_uns + 1L, sample_curr * nCondition)

    label_vec_uns <- labelsList[[ind_uns]]
    prop_uns_truth <- sum(grepl("^gp$", label_vec_uns)) / length(label_vec_uns)

    purrr::map_df(ind_stim, function(ind) {
      label_vec_stim <- labelsList[[ind]]
      prop_stim_truth <- sum(grepl("^gp$", label_vec_stim)) /
        length(label_vec_stim)

      tibble::tibble(
        sample = as.character(sample_curr),
        ind = as.character(ind),
        chnl = "F1",
        nCellStimTruth = nCellStim,
        nCellUnsTruth = nCellUns,
        propStimTruth = prop_stim_truth,
        propUnsTruth = prop_uns_truth,
        propRespTruth = prop_stim_truth - prop_uns_truth
      )
    })
  })
}

.est_from_gate_stats <- function(pathProject, nSample, nCondition) {
  path_gate_stats <- file.path(pathProject, "gateStats.rds")
  if (!file.exists(path_gate_stats)) {
    return(tibble::tibble())
  }

  gate_stats <- readRDS(path_gate_stats)
  if (!is.data.frame(gate_stats) || nrow(gate_stats) == 0L) {
    return(tibble::tibble())
  }

  gate_stats <- gate_stats |>
    .add_missing_cols(list(
      ind = NA_character_,
      chnl = "F1"
    ))

  gate_stats |>
    dplyr::mutate(
      ind = as.character(.data$ind),
      sample = as.character(
        ((suppressWarnings(as.numeric(.data$ind)) - 1) %/% nCondition) + 1
      ),
      chnl = as.character(.data$chnl %||% "F1"),
      propRespEst = suppressWarnings(as.numeric(.first_existing(
        dplyr::cur_data_all(),
        c("propBs", "propBsEst", "propRespEst")
      ))),
      threshold = suppressWarnings(as.numeric(.first_existing(
        dplyr::cur_data_all(),
        c("gate", "gateCyt", "threshold")
      ))),
      nCellStim = suppressWarnings(as.numeric(.first_existing(
        dplyr::cur_data_all(),
        c("nCellStim", "n_cell_stim"),
        default = NA_real_
      ))),
      nCellUns = suppressWarnings(as.numeric(.first_existing(
        dplyr::cur_data_all(),
        c("nCellUns", "n_cell_uns"),
        default = NA_real_
      )))
    ) |>
    dplyr::filter(!is.na(.data$ind)) |>
    dplyr::filter(
      (suppressWarnings(as.numeric(.data$ind)) - 1) %% nCondition != 0
    ) |>
    dplyr::select(
      sample,
      ind,
      chnl,
      propRespEst,
      threshold,
      nCellStim,
      nCellUns,
      dplyr::everything()
    )
}

.run_manual_adaptive_freq_bs <- function(
  transformation,
  prob_response,
  n_cell,
  mean_pos_setting,
  mean_pos,
  bias_uns,
  sample_perturbation_sd,
  condition_perturbation_sd,
  cluster_perturbation_sd,
  background_relative_to_response,
  n_cell_uns_relative_to_stim,
  method_family,
  method_id,
  bwAdaptive,
  bw,
  bwAdaptiveCore,
  bwAdaptiveExtra,
  bwAdaptiveCrossover,
  bwAdaptiveTransitionWidth,
  sim_id,
  nSample = n_sample_run,
  nIter = n_iter_run,
  nMarker = n_marker_run,
  nCondition = n_condition_run,
  nCluster = n_cluster_run
) {
  bw_setting <- .scalar_or_null(bw)
  bwAdaptiveCore <- .scalar_or_null(bwAdaptiveCore)
  bwAdaptiveExtra <- .scalar_or_null(bwAdaptiveExtra)
  bwAdaptiveCrossover <- .scalar_or_null(bwAdaptiveCrossover)
  bwAdaptiveTransitionWidth <- .scalar_or_null(bwAdaptiveTransitionWidth) %||% 0

  purrr::map_df(seq_len(nIter), function(iter_num) {
    n_cell_uns <- round(n_cell * n_cell_uns_relative_to_stim)
    n_cell_by_condition <- c(n_cell_uns, n_cell)

    transformation_func <- .simMiscGetTrans(transformation)

    mean_expr_mat <- matrix(
      c(0, mean_pos),
      byrow = TRUE,
      ncol = 1
    )

    cluster_label_vec <- c("gn", "gp")
    prob_response_uns <- prob_response * background_relative_to_response
    prob_vec_uns <- c(1 - prob_response_uns, prob_response_uns)
    prob_response_vec_by_stim_condition <- list(c(
      -prob_response,
      prob_response
    ))

    out_list_experiment <- simCytExperiment(
      nSample = nSample,
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nCellByCondition = n_cell_by_condition,
      transformationFunc = transformation_func,
      mixtureType = "gaussianOnly",
      meanExprMat = mean_expr_mat,
      clusterLabelVec = cluster_label_vec,
      probVecUns = prob_vec_uns,
      probExact = prob_exact,
      probResponseVecByStimCondition = prob_response_vec_by_stim_condition,
      conditionPerturbationSd = condition_perturbation_sd,
      clusterPerturbationSd = cluster_perturbation_sd,
      samplePerturbationSd = sample_perturbation_sd,
      covEvMin = cov_ev_min,
      covEvMax = cov_ev_max
    )

    fs <- as(out_list_experiment[["flowFrameList"]], "flowSet")
    gs <- flowWorkspace::GatingSet(fs)
    labels_list <- out_list_experiment[["labelsList"]]

    path_project <- file.path(
      tempdir(),
      "stimgate",
      "sim-bw",
      "manual-adaptive-perf",
      paste0(
        "sim-",
        sim_id,
        "-iter-",
        iter_num,
        "-pid-",
        Sys.getpid(),
        "-",
        format(Sys.time(), "%Y%m%d%H%M%S"),
        "-",
        sample.int(1e7, 1)
      )
    )

    on.exit(
      {
        if (dir.exists(path_project)) {
          unlink(path_project, recursive = TRUE)
        }
      },
      add = TRUE
    )

    if (dir.exists(path_project)) {
      unlink(path_project, recursive = TRUE)
    }
    dir.create(path_project, recursive = TRUE, showWarnings = FALSE)

    batch_list <- lapply(seq_len(nSample), function(i) {
      seq((i - 1L) * nCondition + 1L, i * nCondition)
    })

    Sys.setenv("STIMGATE_INTERMEDIATE" = "TRUE")
    invisible(gateStim(
      .data = gs,
      pathProject = path_project,
      popGate = "root",
      batchList = batch_list,
      marker = paste0("MarkerF", seq_len(nMarker)),
      calcCytPosGates = calc_cyt_pos_gates,
      calcSinglePosGates = calc_single_pos_gates,
      biasUns = bias_uns,
      bw = bw_setting,
      bwFallback = if (!is.null(bw_setting)) bw_setting else "auto",
      bwMin = bw_min,
      bwMax = bw_max,
      bwMtd = bw_mtd,
      bwAdj = bw_adj,
      bwNcellMin = bw_ncell_min,
      bwNcellMax = bw_ncell_max,
      bwCluster = bw_cluster,
      bwAdaptive = isTRUE(bwAdaptive),
      bwAdaptiveDensityN = bw_adaptive_density_n,
      bwAdaptivePadFrac = bw_adaptive_pad_frac,
      bwAdaptiveCore = bwAdaptiveCore,
      bwAdaptiveExtra = bwAdaptiveExtra,
      bwAdaptiveCrossover = bwAdaptiveCrossover,
      bwAdaptiveTransitionWidth = bwAdaptiveTransitionWidth,
      normPeakFrac = norm_peak_frac,
      normPeakMinRel = norm_peak_min_rel,
      normExtraFrac = norm_extra_frac,
      normExtraMax = norm_extra_max,
      normExtraJitterFrac = norm_extra_jitter_frac,
      normLambda = norm_lambda,
      normDensityN = norm_density_n,
      normExcessBwMtd = norm_excess_bw_mtd,
      normExcessNcell = norm_excess_ncell,
      normAdaptiveNcell = norm_adaptive_ncell,
      normMtd = norm_mtd,
      minCell = min_cell,
      maxPosProbX = max_pos_prob_x,
      gateCombn = gate_combn,
      tolClust = tol_clust
    ))

    truth_tbl <- .truth_from_labels(
      labelsList = labels_list,
      nSample = nSample,
      nCondition = nCondition,
      nCellStim = n_cell,
      nCellUns = n_cell_uns
    )

    est_tbl <- .est_from_gate_stats(
      pathProject = path_project,
      nSample = nSample,
      nCondition = nCondition
    )

    if (nrow(est_tbl) == 0L) {
      est_tbl <- truth_tbl |>
        dplyr::transmute(
          sample = .data$sample,
          ind = .data$ind,
          chnl = .data$chnl,
          propRespEst = NA_real_,
          threshold = NA_real_,
          nCellStim = n_cell,
          nCellUns = n_cell_uns
        )
    }

    truth_tbl |>
      dplyr::left_join(
        est_tbl,
        by = c("sample", "ind", "chnl"),
        suffix = c("", "EstTbl")
      ) |>
      dplyr::mutate(
        iter = iter_num,
        transformation = transformation,
        prob_response = prob_response,
        n_cell = n_cell,
        mean_pos_setting = mean_pos_setting,
        mean_pos = mean_pos,
        bias_uns = bias_uns,
        sample_perturbation_sd = sample_perturbation_sd,
        condition_perturbation_sd = condition_perturbation_sd,
        cluster_perturbation_sd = cluster_perturbation_sd,
        background_relative_to_response = background_relative_to_response,
        n_cell_uns_relative_to_stim = n_cell_uns_relative_to_stim,
        method_family = method_family,
        method_id = method_id,
        bwAdaptive = isTRUE(bwAdaptive),
        bw = bw_setting %||% NA_real_,
        bwAdaptiveCore = bwAdaptiveCore %||% NA_real_,
        bwAdaptiveExtra = bwAdaptiveExtra %||% NA_real_,
        bwAdaptiveCrossover = bwAdaptiveCrossover %||% NA_real_,
        bwAdaptiveTransitionWidth = bwAdaptiveTransitionWidth %||% NA_real_,
        sim_id = sim_id
      ) |>
      dplyr::select(
        transformation,
        prob_response,
        n_cell,
        mean_pos_setting,
        mean_pos,
        bias_uns,
        sample_perturbation_sd,
        condition_perturbation_sd,
        cluster_perturbation_sd,
        background_relative_to_response,
        n_cell_uns_relative_to_stim,
        method_family,
        method_id,
        bwAdaptive,
        bw,
        bwAdaptiveCore,
        bwAdaptiveExtra,
        bwAdaptiveCrossover,
        bwAdaptiveTransitionWidth,
        sim_id,
        iter,
        sample,
        ind,
        chnl,
        propRespTruth,
        propRespEst,
        threshold,
        dplyr::everything()
      )
  })
}
