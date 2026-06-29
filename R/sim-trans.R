sim_trans_univariate_one <- function(
  transformation,
  n_cell = main_settings$n_cell[[1]],
  mean_pos = main_settings$mean_pos[[1]],
  prob_response = main_settings$prob_response[[1]],
  background_relative_to_response = main_settings$background_relative_to_response[[
    1
  ]],
  prob_exact = main_settings$prob_exact[[1]],
  mixture_type = main_settings$mixture_type[[1]],
  cluster_perturbation_sd = main_settings$cluster_perturbation_sd[[1]],
  cov_ev_min = main_settings$cov_ev_min[[1]],
  cov_ev_max = main_settings$cov_ev_max[[1]]
) {
  trans_func <- .simMiscGetTrans(transformation)

  prob_background <- prob_response * background_relative_to_response
  prob_uns <- c(1 - prob_background, prob_background)
  prob_stim <- prob_uns + c(-prob_response, prob_response)

  mean_expr_mat <- matrix(c(0, mean_pos), ncol = 1)
  cluster_label_vec <- c("negative", "response")

  cond_list <- list(
    unstimulated = prob_uns,
    stimulated = prob_stim
  )

  purrr::imap_dfr(cond_list, function(prob_vec, condition_nm) {
    out <- simCytCondition(
      nMarker = 1L,
      nCell = n_cell,
      transformationFunc = trans_func,
      mixtureType = mixture_type,
      meanExprMat = mean_expr_mat,
      clusterLabelVec = cluster_label_vec,
      probVec = prob_vec,
      probExact = prob_exact,
      clusterPerturbationSd = cluster_perturbation_sd,
      covEvMin = cov_ev_min,
      covEvMax = cov_ev_max
    )

    tibble::tibble(
      transformation = transformation,
      condition = condition_nm,
      response_class = out$conditionLabels,
      F1 = as.numeric(out$conditionMatrix[, "F1"])
    )
  })
}

sim_trans_bivariate_one <- function(
  transformation,
  n_cell = main_settings$n_cell[[1]],
  mean_pos = main_settings$mean_pos[[1]],
  prob_response = main_settings$prob_response[[1]],
  background_relative_to_response = main_settings$background_relative_to_response[[
    1
  ]],
  prob_exact = main_settings$prob_exact[[1]],
  mixture_type = main_settings$mixture_type[[1]],
  cluster_perturbation_sd = main_settings$cluster_perturbation_sd[[1]],
  cov_ev_min = main_settings$cov_ev_min[[1]],
  cov_ev_max = main_settings$cov_ev_max[[1]]
) {
  trans_func <- .simMiscGetTrans(transformation)

  prob_background <- prob_response * background_relative_to_response
  prob_background_each <- prob_background / 3
  prob_response_each <- prob_response / 3

  prob_vec <- c(
    1 - prob_background - prob_response,
    prob_background_each + prob_response_each,
    prob_background_each + prob_response_each,
    prob_background_each + prob_response_each
  )

  mean_expr_mat <- matrix(
    c(
      0,
      0,
      mean_pos,
      0,
      0,
      mean_pos,
      mean_pos,
      mean_pos
    ),
    byrow = TRUE,
    ncol = 2
  )

  cluster_label_vec <- c(
    "negative",
    "F1 response",
    "F2 response",
    "F1 and F2 response"
  )

  out <- simCytCondition(
    nMarker = 2L,
    nCell = n_cell,
    transformationFunc = trans_func,
    mixtureType = mixture_type,
    meanExprMat = mean_expr_mat,
    clusterLabelVec = cluster_label_vec,
    probVec = prob_vec,
    probExact = prob_exact,
    clusterPerturbationSd = cluster_perturbation_sd,
    covEvMin = cov_ev_min,
    covEvMax = cov_ev_max
  )

  tibble::tibble(
    transformation = transformation,
    condition = "stimulated",
    response_class = out$conditionLabels
  ) |>
    dplyr::bind_cols(as.data.frame(out$conditionMatrix))
}

sim_trans_downsample_for_display <- function(
  .data,
  background_n = 12000,
  response_n = Inf
) {
  background_tbl <- .data |>
    dplyr::filter(.data$response_class == "negative")

  response_tbl <- .data |>
    dplyr::filter(.data$response_class != "negative")

  if (nrow(background_tbl) > background_n) {
    background_tbl <- dplyr::slice_sample(background_tbl, n = background_n)
  }

  if (is.finite(response_n) && nrow(response_tbl) > response_n) {
    response_tbl <- dplyr::slice_sample(response_tbl, n = response_n)
  }

  dplyr::bind_rows(background_tbl, response_tbl)
}

plot_univariate_transformation <- function(plot_tbl, transformation) {
  trans_curr <- transformation
  x_ref <- .simMiscGetTrans(trans_curr)(
    c(0, main_settings$mean_pos[[1]])
  )

  plot_tbl |>
    dplyr::filter(.data$transformation == trans_curr) |>
    ggplot(
      aes(
        x = .data$F1,
        colour = .data$condition,
        fill = .data$condition
      )
    ) +
    geom_density(alpha = 0.2, linewidth = 0.6) +
    geom_rug(
      data = function(x) dplyr::filter(x, .data$response_class != "negative"),
      aes(x = .data$F1),
      sides = "b",
      alpha = 0.35,
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = x_ref, linetype = "dotted") +
    scale_y_sqrt() +
    labs(
      title = trans_curr,
      x = "F1 expression",
      y = "Density, square-root scale",
      colour = NULL,
      fill = NULL
    ) +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "xy", minor = "none") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

plot_bivariate_transformation <- function(plot_tbl, transformation) {
  trans_curr <- transformation
  ref <- .simMiscGetTrans(trans_curr)(
    c(0, main_settings$mean_pos[[1]])
  )

  plot_tbl |>
    dplyr::filter(.data$transformation == trans_curr) |>
    sim_trans_downsample_for_display() |>
    ggplot(
      aes(
        x = .data$F1,
        y = .data$F2,
        colour = .data$response_class
      )
    ) +
    geom_point(size = 0.25, alpha = 0.35) +
    geom_hline(yintercept = ref, linetype = "dotted") +
    geom_vline(xintercept = ref, linetype = "dotted") +
    coord_equal() +
    labs(
      title = trans_curr,
      x = "F1 expression",
      y = "F2 expression",
      colour = NULL
    ) +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "xy", minor = "none") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}
