.get_cp_cluster_control_update <- function(control) {
  if (!"min_threshold_frac" %in% names(control)) {
    control[["min_threshold_frac"]] <- 0.8
  }
  if (!"min_threshold_quant" %in% names(control)) {
    control[["min_threshold_quant"]] <- 0.1
  }
  control
}

.get_cp_cluster_gate_stats_tbl_update <- function(gate_stats_tbl,
                                                  debug = FALSE) {
  .debug(debug, "Updating gate statistics table") # nolint
  gate_stats_tbl |>
    dplyr::mutate(
      prop_stim_pos = pmax(count_stim, 1) / n_cell_stim, # nolint
      prop_uns_pos = pmax(count_uns, 1) / n_cell_uns, # nolint
      prop_stim_sd = sqrt(
        prop_stim_pos * (1 - prop_stim_pos) / n_cell_stim # nolint
      ),
      prop_uns_sd = sqrt(
        prop_uns_pos * (1 - prop_uns_pos) / n_cell_uns # nolint
      ),
      prop_bs_sd = sqrt(
        prop_stim_sd^2 + prop_uns_sd^2 # nolint
      )
    )
}

.get_cp_cluster_cp_get_min <- function(gate_tbl, gate_tbl_ctrl) {
  suppressWarnings(
    min(gate_tbl$gate_cyt, gate_tbl$gate, gate_tbl_ctrl$gate, na.rm = TRUE)
  )
}

.get_cp_cluster_cp_get_max <- function(gate_tbl, gate_tbl_ctrl) {
  suppressWarnings(
    max(gate_tbl$gate, gate_tbl$gate_single, gate_tbl_ctrl$gate, na.rm = TRUE)
  )
}

.get_cp_cluster_prop_bs_by_cp_tbl_obj <- function(gs,
                                                  gate_tbl,
                                                  ind_batch_list,
                                                  ind_in_batch_lab_vec,
                                                  ind_in_batch_uns,
                                                  high,
                                                  cut,
                                                  pop_gate,
                                                  data_name,
                                                  calc_cyt_pos_gates,
                                                  cp_min,
                                                  max_cp,
                                                  gate_stats_tbl,
                                                  filter_other_cyt_pos,
                                                  debug) {
  .debug(debug, "Getting prop_bs_by_cp_tbl object") # nolint
  # statistics
  # ----------------

  data_list_obj <- .get_prop_bs_by_cp_tbl_data_list(
    gs = gs,
    gate_tbl = gate_tbl,
    ind_batch_list = ind_batch_list,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    ind_in_batch_uns = ind_in_batch_uns,
    high = high,
    cut = cut,
    pop_gate = pop_gate,
    data_name = data_name,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    max_cp = max_cp,
    filter_other_cyt_pos = filter_other_cyt_pos,
    cp_min = cp_min,
    debug = debug
  )

  prop_bs_by_cp_tbl <- .get_prop_bs_by_cp_tbl_actual(
    data_list = data_list_obj[["data_list"]],
    cp_min = cp_min,
    max_cp = max_cp,
    gate_stats_tbl = gate_stats_tbl,
    chnl = cut,
    debug = debug,
    ind_batch_list = ind_batch_list,
    ind_in_batch_uns = ind_in_batch_uns
  )

  list(
    "prop_bs_by_cp_tbl" = prop_bs_by_cp_tbl,
    "expr_max" = data_list_obj[["expr_max"]],
    "expr_min" = data_list_obj[["expr_min"]]
  )
}

.get_prop_bs_by_cp_tbl_data_list <- function(gs,
                                             gate_tbl,
                                             ind_batch_list,
                                             ind_in_batch_lab_vec,
                                             ind_in_batch_uns,
                                             high,
                                             cut,
                                             pop_gate,
                                             data_name,
                                             calc_cyt_pos_gates,
                                             max_cp,
                                             filter_other_cyt_pos,
                                             cp_min,
                                             debug) {
  .debug(debug, "Getting data list") # nolint
  data_list <- purrr::map(ind_batch_list, function(ind_batch) {
    ex_list <- .get_ex_list( # nolint
      data = gs, # nolint
      ind_batch = ind_batch,
      ind_in_batch_gate = seq_along(ind_in_batch_lab_vec),
      ind_in_batch_uns = ind_in_batch_uns,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec,
      pop = pop_gate,
      cut = names(high), high = NULL,
      data_name = data_name
    )

    # range of expressions
    # (used later when calculating how "flat" a derivative is)
    expr_range_tbl <- purrr::map_df(
      seq_along(ex_list),
      function(i) {
        purrr::map_df(names(high), function(chnl_ind) {
          if (!chnl_ind %in% colnames(ex_list[[i]])) {
            return(NULL)
          }
          quant_vec <- quantile(
            ex_list[[i]][[chnl_ind]], c(0.0025, 0.999)
          )
          tibble::tibble(
            lb = quant_vec[[1]],
            ub = quant_vec[[2]]
          )
        })
      }
    )

    expr_min <- quantile(expr_range_tbl[["lb"]], 0.25)
    expr_max <- max(expr_range_tbl[["ub"]])

    # filter to yield cells negative for all cytokine combinations
    # except possible this cytokine single-positive
    if (filter_other_cyt_pos) {
      ex_list_filter <- purrr::map(seq_along(ex_list), function(i) {
        if (i == ind_in_batch_uns) {
          return(ex_list[[i]])
        }
        gate_tbl_ind <- gate_tbl |>
          dplyr::filter(ind == ex_list[[i]]$ind[1]) # nolint

        pos_ind_vec_but_single_pos_curr <-
          .get_pos_ind_but_single_pos_for_one_cyt( # nolint
            ex = ex_list[[i]],
            gate_tbl = gate_tbl_ind,
            chnl_single_exc = cut,
            chnl = NULL,
            gate_type_cyt_pos = ifelse(calc_cyt_pos_gates,
              "cyt", "base"
            ),
            gate_type_single_pos = "base"
          )
        ex_list[[i]][!pos_ind_vec_but_single_pos_curr, , drop = FALSE]
      }) |>
        stats::setNames(names(ex_list))
    } else {
      ex_list_filter <- ex_list
    }

    # this is just to save memory
    # for when calculating the per-threshold performance
    out_tbl <- ex_list_filter |>
      purrr::map(function(x) {
        x <- x |>
          dplyr::mutate(n_cell = nrow(x))
        x_out <- x |>
          dplyr::filter(x[["cut"]] >= min(.env$cp_min, max(x[["cut"]]))) # nolint
        if (nrow(x_out) == 0) {
          x_out <- x[1, ] |>
            dplyr::select(batch:stim) # nolint
          x_add <- x[1, ] |>
            dplyr::select(-c(batch:stim)) # nolint
          x_add[1, ] <- NA
          x_out <- x_out |>
            dplyr::bind_cols(x_add)
        }
        x_out
      })
    list(
      "out_tbl" = out_tbl,
      "expr_min" = expr_min,
      "expr_max" = expr_max
    )
  })
  expr_min_vec <- sapply(data_list, function(x) x$expr_min)
  expr_max_vec <- sapply(data_list, function(x) x$expr_max)
  expr_min <- min(expr_min_vec)
  expr_max <- max(
    max(expr_max_vec),
    max_cp + 0.2 * (max(expr_max_vec) - expr_min)
  )
  data_list <- lapply(data_list, function(x) x$out_tbl) |>
    purrr::flatten()
  list(
    "data_list" = data_list,
    "expr_min" = expr_min,
    "expr_max" = expr_max
  )
}




.get_prop_bs_by_cp_tbl_actual <- function(data_list,
                                          cp_min,
                                          max_cp,
                                          gate_stats_tbl,
                                          debug,
                                          chnl,
                                          ind_batch_list,
                                          ind_in_batch_uns) {
  .debug(debug, "Getting prop_bs_by_cp_tbl") # nolint
  cp_par_list <- .get_prop_bs_by_cp_tbl_actual_prep_cp(cp_min, max_cp)
  purrr::map(seq_along(data_list), function(i) {
    .get_prop_bs_by_cp_tbl_actual_progress(debug, i, data_list)
    ex_list <- .get_prop_bs_by_cp_tbl_actual_ex_get(
      data_list, i, ind_batch_list, ind_in_batch_uns
    )
    if (is.null(ex_list)) {
      return(NULL)
    }
    .get_prop_bs_by_cp_tbl(
      ex_stim = ex_list$stim,
      ex_uns = ex_list$uns,
      cp_seq = cp_par_list[["seq"]],
      gate_stats_tbl = gate_stats_tbl,
      chnl = chnl,
      debug = debug
    )
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()
}

.get_prop_bs_by_cp_tbl_prep <- function(gate_stats_tbl,
                                        ex_stim,
                                        ex_uns,
                                        cp_seq,
                                        chnl) {
  gate_stats_tbl_curr <- gate_stats_tbl |>
    dplyr::filter(.data$ind == ex_stim$ind[1]) # nolint
  count_stim_vec <- rep(NA, length(cp_seq))
  count_uns_vec <- rep(NA, length(cp_seq))
  for (i in seq_along(cp_seq)) {
    cp <- cp_seq[i]
    count_stim_vec[i] <- sum(ex_stim[[chnl]] > cp) # nolint
    count_uns_vec[i] <- sum(ex_uns[[chnl]] > cp) # nolint
  }
  prop_bs_sd <- gate_stats_tbl_curr$prop_bs_sd
  prop_bs_orig <- gate_stats_tbl_curr$prop_bs
  list(
    "gate_stats" = gate_stats_tbl_curr,
    "count_stim" = count_stim_vec,
    "count_uns" = count_uns_vec,
    "bs_sd" = prop_bs_sd,
    "bs_orig" = prop_bs_orig
  )
}

.get_prop_bs_by_cp_tbl <- function(ex_stim,
                                   ex_uns,
                                   cp_seq,
                                   gate_stats_tbl,
                                   chnl,
                                   debug) {
  par_list <- .get_prop_bs_by_cp_tbl_prep(
    gate_stats_tbl, ex_stim, ex_uns, cp_seq, chnl
  )

  .get_prop_bs_by_cp_tbl_init(ex_stim, ex_uns, par_list, cp_seq, chnl) |>
    .get_prop_bs_by_cp_tbl_calc(ex_stim$n_cell[1], ex_uns$n_cell[1])
}

.get_prop_bs_by_cp_tbl_actual_prep_cp <- function(cp_min, cp_max) {
  cp_range <- c(cp_min, cp_max)
  cp_seq_vec <- seq(cp_range[1], cp_range[2], length.out = 1e2)
  list("range" = cp_range, "seq" = cp_seq_vec)
}

.get_prop_bs_by_cp_tbl_init <- function(ex_stim,
                                        ex_uns,
                                        par_list,
                                        cp_seq,
                                        chnl) {
  tibble::tibble(
    ind = ex_stim$ind[1], stim = ex_stim$stim[1],
    prop_bs_orig = par_list[["bs_orig"]], prop_bs_sd = par_list[["bs_sd"]],
    cp = cp_seq, max_expr = max(ex_stim[[chnl]], ex_uns[[chnl]]), # nolint
    count_stim_cp = par_list[["count_stim"]],
    count_uns_cp = par_list[["count_uns"]]
  )
}

.get_prop_bs_by_cp_tbl_calc <- function(.data, n_cell_stim, n_cell_uns) {
  .data |>
    dplyr::mutate(
      prop_stim_cp = count_stim_cp / n_cell_stim, # nolint
      prop_uns_cp = count_uns_cp / n_cell_uns # nolint
    ) |>
    dplyr::mutate(
      prop_bs_cp = prop_stim_cp - prop_uns_cp, # nolint
      prop_bs_cp_diff = prop_bs_cp - prop_bs_orig, # nolint
      prop_bs_cp_diff_sd = prop_bs_cp_diff / prop_bs_sd, # nolint
      prop_stim_pos_cp = pmax(count_stim_cp, 1) / n_cell_stim,
      prop_uns_pos_cp = pmax(count_uns_cp, 1) / n_cell_uns,
      prop_stim_sd_cp = sqrt(
        prop_stim_pos_cp * (1 - prop_stim_pos_cp) / n_cell_stim # nolint
      ),
      prop_uns_sd_cp = sqrt(
        prop_uns_pos_cp * (1 - prop_uns_pos_cp) / n_cell_uns # nolint
      ),
      prop_bs_sd_cp = sqrt(
        prop_stim_sd_cp^2 + prop_uns_sd_cp^2 # nolint
      ),
      prop_bs_cp_diff_sd_max = prop_bs_cp_diff /
        pmax(prop_bs_sd, prop_bs_sd_cp) # nolint
    )
}

.get_prop_bs_by_cp_tbl_actual_progress <- function(debug, i, data_list) {
  if (i %% 20 == 0) {
    .debug(debug, paste0("Processing ", i, " of ", length(data_list))) # nolint
  }
}

.get_prop_bs_by_cp_tbl_actual_progress <- function(debug, i, data_list) {
  if (i %% 20 == 0) {
    .debug(debug, paste0("Processing ", i, " of ", length(data_list))) # nolint
  }
}

.get_prop_bs_by_cp_tbl_actual_ex_get <- function(data_list,
                                                 i,
                                                 ind_batch_list,
                                                 ind_in_batch_uns) {
  ex <- data_list[[i]]
  if (.get_prop_bs_by_cp_return_early_stim(ex)) {
    return(NULL)
  }

  ex_uns <- .get_cp_cluster_tbl_ex_uns(
    ind_stim = ex$ind[1], ind_batch_list = ind_batch_list,
    ind_in_batch_uns = ind_in_batch_uns, data_list = data_list
  )

  if (.get_prop_bs_by_cp_return_early_uns(ex_uns)) {
    return(NULL)
  }

  list("stim" = ex, "uns" = ex_uns)
}

.get_prop_bs_by_cp_return_early_stim <- function(ex_stim) {
  if (is.na(ex_stim$cut[1])) {
    return(TRUE)
  }
  ex_stim$stim[1] == "uns"
}

.get_cp_cluster_tbl_ex_uns <- function(ind_stim,
                                       ind_batch_list,
                                       ind_in_batch_uns,
                                       data_list) {
  ex_uns_ind <- .get_cp_cluster_tbl_ex_uns_ind(
    ind_stim, ind_batch_list, ind_in_batch_uns
  )
  .get_cp_cluster_tbl_ex_uns_get(data_list, ex_uns_ind)
}
.get_cp_cluster_tbl_ex_uns_ind <- function(ind_stim,
                                           ind_batch_list,
                                           ind_in_batch_uns) {
  batch_vec <- .get_batch_from_ind(ind_stim, ind_batch_list) # nolint
  batch_vec[ind_in_batch_uns]
}

.get_cp_cluster_tbl_ex_uns_get <- function(data_list,
                                           ex_uns_ind) {
  data_list[[
    purrr::map_lgl(data_list, ~ .env$ex_uns_ind == .x$ind)
  ]]
}

.get_prop_bs_by_cp_return_early_uns <- function(ex_uns) {
  is.na(ex_uns$cut[1])
}

.get_cp_cluster_dens_tbl_get <- function(ind_batch_list,
                                         gs,
                                         ind_in_batch_lab_vec,
                                         ind_in_batch_uns,
                                         high,
                                         data_name,
                                         filter_other_cyt_pos,
                                         calc_cyt_pos_gates,
                                         cut,
                                         expr_min,
                                         expr_max,
                                         pop_gate,
                                         gate_tbl,
                                         control,
                                         bw,
                                         debug) {
  .debug(debug, "Getting density table") # nolint
  purrr::map_df(ind_batch_list, function(ind_batch) {
    .debug(debug, paste0("Processing batch ", ind_batch)) # nolint
    ex_list <- .get_ex_list( # nolint
      data = gs, ind_batch = ind_batch,
      ind_in_batch_gate = seq_along(
        ind_in_batch_lab_vec
      ),
      ind_in_batch_uns = ind_in_batch_uns,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec,
      pop = pop_gate,
      cut = names(high), high = NULL,
      data_name = data_name
    )

    # filter to yield cells negative for all cytokine combinations
    # except possible this cytokine single-positive
    if (filter_other_cyt_pos) {
      .debug(debug, "Filtering other cytokine positive cells") # nolint
      ex_list_filter <- purrr::map(seq_along(ex_list), function(i) {
        if (i == ind_in_batch_uns) {
          return(ex_list[[i]])
        }

        gate_tbl_ind <- gate_tbl |>
          dplyr::filter(ind == ex_list[[i]]$ind[1]) # nolint

        pos_ind_vec_but_single_pos_curr <-
          .get_pos_ind_but_single_pos_for_one_cyt( # nolint
            ex = ex_list[[i]],
            gate_tbl = gate_tbl_ind,
            chnl_single_exc = cut,
            chnl = NULL,
            gate_type_cyt_pos = ifelse(calc_cyt_pos_gates,
              "cyt", "base"
            ),
            gate_type_single_pos = "base"
          )
        ex_list[[i]][!pos_ind_vec_but_single_pos_curr, , drop = FALSE]
      }) |>
        stats::setNames(names(ex_list))
    } else {
      ex_list_filter <- ex_list
    }

    ex_list_filter <- ex_list_filter[-ind_in_batch_uns]

    # calculate x_ind so as to make filtering
    # easier later
    min_threshold_gate_quant <- quantile(
      gate_tbl$gate, control$min_threshold_quant,
      na.rm = TRUE
    )
    min_threshold <- control$min_threshold_frac * min_threshold_gate_quant

    purrr::map_df(ex_list_filter, function(x) {
      if (length(x[[cut]][x[[cut]] > min(x[[cut]])]) < 3) { # nolint
        return(tibble::tibble(
          batch_sh = x$batch_sh[1],
          stim = x$stim[1],
          ind = x$ind[1],
          y = rep(NA, 512),
          x = paste0("x", seq.int(from = 1, to = 512))
        ) |>
          tidyr::pivot_wider(
            names_from = x,
            values_from = y # nolint
          ))
      }
      dens <- density(
        x[[cut]][x[[cut]] > min(x[[cut]])], # nolint
        from = expr_min, to = expr_max, bw = bw
      )
      tibble::tibble(
        batch_sh = x$batch_sh[1],
        stim = x$stim[1],
        ind = x$ind[1],
        y = dens$y, # nolint
        x = dens$x,
        x_ind = paste0("x", seq_len(length(dens$y)))
      ) |>
        dplyr::filter(
          x <= min_threshold
        ) |>
        dplyr::select(-x) |>
        dplyr::mutate(
          y = y / sum(y)
        ) |>
        tidyr::pivot_wider(
          names_from = x_ind, # nolint
          values_from = y
        )
    })
  }) |>
    dplyr::filter(!is.na(x1)) # nolint
}

.get_cp_cluster_n_clus <- function(dens_tbl) {
  max_cluster <- min(6, nrow(dens_tbl) / 3) |>
    floor() |>
    max(1)
  if (max_cluster == 1L) {
    return(1L)
  }
  dens_mat <- dens_tbl[, grepl("^x\\d+", colnames(dens_tbl))] |>
    as.matrix()
  clus_gap_obj <- cluster::clusGap(
    dens_mat,
    FUNcluster = kmeans,
    K.max = max_cluster,
    B = 50
  )
  cluster::maxSE(
    clus_gap_obj$Tab[, "gap"], clus_gap_obj$Tab[, "SE.sim"]
  )
}

.get_cp_cluster_plot_check_1 <- function(dens_tbl) {
  dens_plot <- dens_tbl |>
    # dplyr::group_by(grp) |>
    # dplyr::slice(1) |>
    # dplyr::ungroup() |>
    tidyr::pivot_longer(
      names_to = "x",
      values_to = "y",
      x1:x512 # nolint
    ) |>
    dplyr::group_by(grp, ind) |> # nolint
    dplyr::mutate(x = x_vec) |> # nolint
    dplyr::ungroup()


  ggplot( # nolint
    dens_plot,
    aes(x, y, col = grp, group = ind) # nolint
  ) +
    geom_line() # nolint
}

.get_cp_cluster_plot_check_1lse <- function(prop_bs_by_cp_tbl) {
  data_plot <- prop_bs_by_cp_tbl |>
    dplyr::group_by(grp, cp) |> # nolint
    dplyr::summarise(
      prop_l1se = sum(prop_bs_cp_diff_sd < 1) / dplyr::n() # nolint
    ) |>
    dplyr::ungroup()

  ggplot(data_plot, aes(x = cp, y = prop_l1se, col = grp)) + # nolint
    geom_line() + # nolint
    geom_smooth( # nolint
      method = "gam",
      formula = y ~ s(x, bs = "cs"),
      se = FALSE
    )
}

.get_cp_cluster_clus <- function(dens_tbl, n_clus) {
  dens_mat <- dens_tbl[, grepl("^x\\d+", colnames(dens_tbl))] |>
    as.matrix()
  stats::kmeans(dens_mat, centers = n_clus)$cluster
}

.get_cp_cluster_data_mod_pre <- function(prop_bs_by_cp_tbl) {
  data_mod_pre <- prop_bs_by_cp_tbl |>
    tidyr::pivot_longer(-c(ind:prop_bs_cp_diff_sd_max), # nolint
      names_to = "grp", values_to = "grp_level"
    )

  # so here we just get each individual's
  # group for each number of groups.
  # I wonder if the need to slice is coming
  # from
  data_mod_filter_vec <- data_mod_pre |>
    dplyr::group_by(ind, grp, grp_level) |> # nolint
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::group_by(grp, grp_level) |>
    dplyr::arrange(grp, grp_level, ind) |>
    dplyr::summarise(ind_vec = paste0(ind, collapse = "_")) |>
    dplyr::ungroup() |>
    dplyr::group_by(ind_vec) |> # nolint
    dplyr::slice(1) |>
    dplyr::mutate(grp_grp_level = paste0(grp, "_", grp_level)) |>
    dplyr::pull("grp_grp_level")

  data_mod <- data_mod_pre |>
    dplyr::filter(
      paste0(grp, "_", grp_level) %in% data_mod_filter_vec # nolint
    ) |> # nolint
    dplyr::group_by(grp, grp_level, cp) |> # nolint
    dplyr::summarise(
      prop_l1se =
        sum(prop_bs_cp_diff_sd < 1, na.rm = TRUE) / # nolint
          sum(!is.na(prop_bs_cp_diff_sd))
    ) |>
    dplyr::ungroup()

  data_mod |>
    dplyr::group_by(grp, grp_level) |> # nolint
    dplyr::filter(sum(is.na(prop_l1se)) != dplyr::n()) |> # nolint
    dplyr::ungroup()
}

.get_cp_cluster_cp_grp_lab_vec_get <- function(prop_bs_by_cp_tbl,
                                               expr_min,
                                               expr_max,
                                               debug = FALSE) {
  data_mod <- .get_cp_cluster_data_mod_pre(
    prop_bs_by_cp_tbl = prop_bs_by_cp_tbl
  )
  purrr::map(unique(data_mod$grp), function(n_grp_curr) {
    .debug(debug, paste0("Processing cluster ", n_grp_curr)) # nolint
    data_mod_curr <- data_mod |>
      dplyr::filter(.data$grp == n_grp_curr) # nolint
    purrr::map( # nolint
      unique(data_mod_curr$grp_level), function(grp_level) {
        data_mod_curr_grp <- data_mod_curr |>
          dplyr::filter(grp_level == .env$grp_level) # nolint
        data_mod_curr_grp_not_na <- data_mod_curr_grp |>
          dplyr::filter(!is.na(prop_l1se)) # nolint
        if (nrow(data_mod_curr_grp_not_na) <= 5) {
          return(max(data_mod_curr_grp$cp, na.rm = TRUE))
        }
        fit <- try(
          suppressWarnings(
            mgcv::gam(prop_l1se ~ s(cp, bs = "cs"),
              data = data_mod_curr_grp_not_na,
              family = mgcv::betar(link = "logit"),
              maxit = 20
            )
          ),
          silent = TRUE
        )
        if (inherits(fit, "try-error")) {
          return(max(data_mod_curr_grp$cp, na.rm = TRUE))
        }
        data_mod_curr_crp <- data_mod_curr_grp |> # nolint
          dplyr::mutate(
            pred = predict(
              fit, data_mod_curr_grp,
              type = "response"
            )
          )
        cp_range <- range(data_mod_curr_grp$cp)
        min_cp_permitted <- data_mod_curr_grp$cp[
          data_mod_curr_grp$prop_l1se > 0.5
        ] |>
          min(na.rm = TRUE)
        if (abs(min_cp_permitted) == Inf) {
          return(max(data_mod_curr_grp$cp, na.rm = TRUE))
        }
        data_pred <- tibble::tibble(
          cp = seq(min_cp_permitted, cp_range[2],
            length.out = 1e5
          )
        )
        data_pred <- data_pred |> dplyr::mutate(
          pred = predict(fit, newdata = data_pred, type = "response")
        )
        # okay, so we calculate the rate of
        # increase of the derivative
        data_pred_der <- data_pred[-1, ] |>
          dplyr::mutate(
            der = (pred - data_pred$pred[-nrow(data_pred)]) / # nolint
              (cp - data_pred$cp[-nrow(data_pred)]) # nolint
          )
        # we also then calculate
        # a 1/100th of the range
        # we divide this into 1/1000t
        max_der <- 0.1 / diff(c(expr_max, expr_min))
        # for each one percent change in cp_range,
        # there is allowed a 0.001 change in prob
        # issue may occur if cp_range is
        # very limited - perhaps better to use expression range?
        # we then choose any parts where the derivative is
        # positive and is sufficiently flat.
        cp_der <- data_pred_der |>
          dplyr::filter(der > 0) |> # nolint
          dplyr::filter(der < max_der) |>
          dplyr::pull("cp") |>
          min()
        cp_cdf <- data_pred_der |>
          dplyr::filter(pred > 0.85) |> # nolint
          dplyr::pull("cp") |>
          min()
        min(cp_der, cp_cdf)
        # need to decide what to do if it doesn't work
        # cluster hierarchically, and then  adjust?
      }
    ) |>
      stats::setNames(paste0(n_grp_curr, "_", unique(data_mod_curr$grp_level)))
  }) |>
    purrr::flatten() |>
    unlist()
}

.get_cp_cluster_gate_summ_stat_tbl_get <- function(gate_tbl,
                                                   cut,
                                                   grp_ind_lab_vec,
                                                   debug) {
  .debug( # nolint
    debug, "Getting quantiles of original gates per clustered observations" # nolint
  )
  if ("chnl" %in% names(gate_tbl)) {
    gate_tbl <- gate_tbl |>
      dplyr::filter(.data$chnl == .env$cut) # nolint
  }
  gate_tbl |>
    dplyr::mutate(grp = grp_ind_lab_vec[as.character(ind)]) |> # nolint
    dplyr::group_by(grp) |> # nolint
    dplyr::summarise(
      gate_05 = quantile(gate, 0.05, na.rm = TRUE), # nolint
      gate_10 = quantile(gate, 0.1, na.rm = TRUE), # nolint
      gate_25 = quantile(gate, 0.25, na.rm = TRUE), # nolint
      gate_50 = median(gate, na.rm = TRUE), # nolint
      gate_75 = quantile(gate, 0.75, na.rm = TRUE), # nolint
      gate_90 = quantile(gate, 0.9, na.rm = TRUE), # nolint
      gate_95 = quantile(gate, 0.95, na.rm = TRUE), # nolint
      .groups = "drop"
    )
}

.get_cp_cluster_cp_join_get <- function(prop_bs_by_cp_tbl,
                                        cp_grp_lab_vec,
                                        debug) {
  .debug(debug, "Getting cp_join") # nolint
  prop_bs_by_cp_tbl |>
    dplyr::group_by(ind) |> # nolint
    dplyr::mutate(
      cp_join = cp_grp_lab_vec[paste0("grp_", grp)] # nolint
    ) |>
    dplyr::ungroup()
}

.get_cp_cluster_gate_tbl_chnl_get <- function(gate_tbl, chnl) {
  if ("chnl" %in% colnames(gate_tbl)) {
    gate_tbl_chnl <- gate_tbl |>
      dplyr::filter(.data$chnl == .env$chnl) # nolint
  } else {
    gate_tbl_chnl <- gate_tbl
  }
  gate_tbl_chnl
}

.get_cp_cluster_cp_tbl_add_info <- function(cp_tbl,
                                            gate_stats_tbl,
                                            gate_summ_stat_tbl,
                                            gate_tbl_ctrl,
                                            gate_tbl_chnl,
                                            debug) {
  .debug(debug, "Adding information to cp table") # nolint
  cp_tbl |>
    dplyr::left_join(
      gate_tbl_chnl |>
        dplyr::select(gate, ind) |> # nolint
        dplyr::rename(cp_orig = gate),
      by = "ind"
    ) |>
    dplyr::left_join(
      gate_stats_tbl |>
        dplyr::select(ind, freq_bs, freq_stim) |> # nolint
        dplyr::rename(
          freq_orig = freq_stim,
          freq_bs_orig = freq_bs
        ),
      by = "ind"
    ) |>
    dplyr::left_join(gate_summ_stat_tbl,
      by = "grp"
    ) |>
    dplyr::left_join(
      gate_tbl_ctrl |>
        dplyr::select(ind, gate) |>
        dplyr::filter(!is.na(gate)) |>
        # dplyr::group_by(ind) |>
        # dplyr::filter(gate == min(gate, na.rm = TRUE)) |>
        dplyr::rename(cp_tg_ctrl = gate),
      by = "ind"
    )
}

.get_cp_cluster_cp_tbl_add_orig_quant_min <- function(cp_tbl, debug) {
  .debug(debug, "Adding original and minimum quantile threshold") # nolint
  cp_tbl |>
    dplyr::mutate(
      cp_orig_quant_min = pmax(pmin(cp_orig, max_expr), gate_05) # nolint
    )
}
.get_cp_cluster_cp_join_lse_get <- function(cp_tbl, debug) {
  .debug(debug, "Getting cp_join_lse") # nolint
  cp_tbl <- cp_tbl |>
    dplyr::group_by(ind) |> # nolint
    dplyr::mutate(
      lse_orig = prop_bs_cp_diff_sd < 0.01, # nolint
      cp_join_lse = min(cp[lse_orig & cp >= cp_join]), # nolint
      cp_join_lse_orig = pmin(cp_join_lse, cp_orig_quant_min) # nolint
    ) |>
    dplyr::ungroup()
}

.get_cp_cluster_cp_join_tg_get <- function(cp_tbl, debug) {
  .debug(debug, "Getting cp_join_tg") # nolint
  cp_tbl |>
    dplyr::group_by(ind) |> # nolint
    dplyr::mutate(
      cp_join_tg = min(cp[ # nolint
        cp >= cp_join & cp >= cp_tg_ctrl & # nolint
          prop_bs_cp_diff_sd <= 2 # nolint
      ]),
      cp_join_tg = ifelse(
        is.na(cp_join_tg), cp_join_lse, cp_join_tg # nolint
      ),
      cp_join_tg_orig = pmin(cp_join_tg, cp_orig_quant_min) # nolint
    ) |>
    dplyr::ungroup()
}

.get_cp_cluster_cp_lse_orig_mean <- function(cp_tbl, debug) {
  .debug(debug, "Getting cp_join_lse_orig_mean") # nolint
  cp_tbl |>
    dplyr::mutate(
      cp_join_lse_orig_mean = pmin(
        cp_orig_quant_min, # nolint
        cp_join_lse_orig + # nolint
          (cp_orig_quant_min - cp_join_lse_orig) / 2
      )
    )
}

.get_cp_cluster_cp_join_tg_orig_mean <- function(cp_tbl, debug) {
  .debug(debug, "Getting cp_join_tg_orig_mean") # nolint
  cp_tbl |>
    dplyr::mutate(
      cp_join_tg_orig_mean = pmin(
        cp_orig_quant_min, # nolint
        cp_join_tg_orig + (cp_orig_quant_min - cp_join_tg_orig) / 2 # nolint
      )
    )
}

.get_cp_cluster_cp_join_lse_orig_mean_tg <- function(cp_tbl, debug) {
  .debug(debug, "Getting cp_join_lse_orig_mean_tg") # nolint
  cp_tbl |>
    dplyr::mutate(
      cp_join_lse_orig_mean_tg = pmin(
        cp_join_lse_orig_mean, cp_join_tg_orig # nolint
      )
    )
}

.get_cp_cluster_cp_filter_above_cp_join_lse_orig_mean_tg <-
  function(cp_tbl, debug) {
    .debug( # nolint
      debug, "Filtering above cp_join_lse_orig_mean_tg"
    ) # nolint
    cp_tbl |>
      dplyr::filter(cp >= cp_join_lse_orig_mean_tg) |> # nolint
      dplyr::group_by(ind) |> # nolint
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::select(
        grp, ind, cp_orig_quant_min, # nolint
        cp_join, cp_join_lse, cp_join_lse_orig, # nolint
        cp_join_lse_orig_mean, # nolint
        cp_join_tg_orig, cp_join_tg_orig_mean, # nolint
        cp_join_lse_orig_mean_tg,
        prop_bs_orig, prop_bs_cp_diff, # nolint
        prop_bs_cp_diff_sd, prop_bs_cp # nolint
      ) |>
      dplyr::ungroup() |>
      dplyr::group_by(grp) |>
      dplyr::arrange(dplyr::desc(cp_orig_quant_min)) |>
      dplyr::ungroup()
  }

.get_cp_cluster_cp_impute_missing_batch <- function(cp_tbl,
                                                    cut,
                                                    gate_tbl,
                                                    dens_tbl,
                                                    debug) {
  .debug( # nolint
    debug, "considering imputing missing thresholds by batch"
  ) # nolint
  ind_with_missing_gates <- setdiff(
    gate_tbl$ind, cp_tbl$ind[!is.na(cp_tbl$cp_join_tg_orig)]
  )
  ind_with_missing_gates <- ind_with_missing_gates[
    !is.na(ind_with_missing_gates)
  ]
  if (length(ind_with_missing_gates) == 0L) {
    .debug(debug, "no missing thresholds by batch") # nolint
    return(cp_tbl)
  }
  .debug(debug, "imputing missing thresholds by batch") # nolint
  for (ind_curr in ind_with_missing_gates) {
    batch <- gate_tbl |>
      dplyr::filter(ind == ind_curr) |> # nolint
      dplyr::pull("batch")

    gate_tbl_ind <- gate_tbl |>
      dplyr::filter(ind == ind_curr) # nolint
    batch <- gate_tbl_ind |>
      dplyr::pull("batch")

    gate_tbl_batch_ind_vec <- gate_tbl |>
      dplyr::filter(batch == .env$batch) |> # nolint
      dplyr::pull("ind")
    gate_tbl_batch_ind_vec <- gate_tbl_batch_ind_vec[
      !is.na(gate_tbl_batch_ind_vec)
    ]

    cp_tbl_batch <- cp_tbl |>
      dplyr::filter(ind %in% gate_tbl_batch_ind_vec) # nolint

    dens_tbl_ind <- dens_tbl |>
      dplyr::filter(ind == ind_curr) # nolint

    # this means no samples had estimates for this batch
    cp_vec_imp <- cp_tbl_batch$cp_join_tg_orig[
      !is.na(cp_tbl_batch$cp_join_tg_orig)
    ]

    if (length(cp_vec_imp) == 0) next

    gate_impute <- .combine_cp( # nolint
      stats::setNames(cp_vec_imp, paste0("a", seq_along(cp_vec_imp))),
      gate_tbl$gate_combn[1]
    )[[1]][[1]]

    cp_tbl_add <- cp_tbl[1, ] |>
      tidyr::pivot_longer(cols = -grp) |> # nolint
      dplyr::mutate(value = NA_real_) |>
      tidyr::pivot_wider(id_cols = grp) |>
      dplyr::mutate(
        grp = ifelse(
          nrow(dens_tbl_ind) > 0, dens_tbl_ind$grp[1], NA
        ),
        ind = ind_curr,
        cp_join_tg_orig = gate_impute
      )

    cp_tbl <- cp_tbl |> dplyr::filter(ind != ind_curr) # nolint

    cp_tbl <- cp_tbl |> dplyr::bind_rows(cp_tbl_add)
  }
  cp_tbl <- cp_tbl |> dplyr::arrange(ind) # nolint
  cp_tbl
}

.get_cp_cluster_impute_missing_ind <- function(cp_tbl,
                                               cut,
                                               gate_tbl,
                                               dens_tbl,
                                               debug) {
  .debug(debug, "considering imputing missing thresholds individually") # nolint
  ind_with_missing_gates <- setdiff(
    gate_tbl$ind, cp_tbl$ind[!is.na(cp_tbl$cp_join_tg_orig)]
  )
  ind_with_missing_gates <- ind_with_missing_gates[
    !is.na(ind_with_missing_gates)
  ]
  if (length(ind_with_missing_gates) == 0L) {
    .debug(debug, "no missing thresholds individually") # nolint
    return(cp_tbl)
  }
  .debug(debug, "imputing missing thresholds individually") # nolint
  for (ind_curr in ind_with_missing_gates) {
    dens_tbl_ind <- dens_tbl |>
      dplyr::select(ind, grp) |> # nolint
      dplyr::filter(ind == ind_curr)
    # this means that there was no density estimates for this sample
    if (nrow(dens_tbl_ind) == 0) next
    gate_impute <- cp_tbl |>
      dplyr::filter(grp == dens_tbl_ind$grp[1]) |> # nolint
      dplyr::pull("cp_join_tg_orig") |>
      quantile(0.75)

    cp_tbl_add <- cp_tbl[1, ] |>
      tidyr::pivot_longer(cols = -grp) |> # nolint
      dplyr::mutate(value = NA_real_) |>
      tidyr::pivot_wider(id_cols = grp) |>
      dplyr::mutate(
        grp = dens_tbl_ind$grp[1],
        ind = dens_tbl_ind$ind[1],
        cp_join_tg_orig = gate_impute
      )

    cp_tbl <- cp_tbl |> dplyr::filter(ind != ind_curr) # nolint

    cp_tbl <- cp_tbl |> dplyr::bind_rows(cp_tbl_add)
  }
  cp_tbl |> dplyr::arrange(ind) # nolint
}

.get_cp_cluster_impute_missing_final <- function(cp_tbl,
                                                 cut,
                                                 gate_tbl,
                                                 dens_tbl,
                                                 debug) {
  .debug(debug, "considering imputing missing thresholds finally") # nolint
  ind_with_missing_gates <- setdiff(
    gate_tbl$ind, cp_tbl$ind[!is.na(cp_tbl$cp_join_tg_orig)]
  )
  ind_with_missing_gates <- ind_with_missing_gates[
    !is.na(ind_with_missing_gates)
  ]
  if (length(ind_with_missing_gates) == 0L) {
    .debug(debug, "no missing thresholds finally") # nolint
    return(cp_tbl)
  }
  .debug(debug, "imputing missing thresholds finally") # nolint
  for (ind_curr in ind_with_missing_gates) {
    dens_tbl_ind <- dens_tbl |>
      dplyr::select(ind, grp) |> # nolint
      dplyr::filter(ind == ind_curr)
    # this means that there was no density estimates for this sample
    if (nrow(dens_tbl_ind) == 0) next
    gate_impute <- cp_tbl |>
      dplyr::filter(grp == dens_tbl_ind$grp[1]) |> # nolint
      dplyr::pull("cp_join_tg_orig") |>
      quantile(0.75)

    cp_tbl_add <- cp_tbl[1, ] |>
      tidyr::pivot_longer(cols = -grp) |> # nolint
      dplyr::mutate(value = NA_real_) |>
      tidyr::pivot_wider(id_cols = grp) |>
      dplyr::mutate(
        grp = dens_tbl_ind$grp[1],
        ind = dens_tbl_ind$ind[1],
        cp_join_tg_orig = gate_impute
      )

    cp_tbl <- cp_tbl |> dplyr::filter(ind != ind_curr) # nolint

    cp_tbl <- cp_tbl |> dplyr::bind_rows(cp_tbl_add)
  }
  cp_tbl |> dplyr::arrange(ind) # nolint
}

.get_cp_cluster_impute_missing_final_batch <- function(cp_tbl,
                                                       cut,
                                                       gate_tbl,
                                                       dens_tbl,
                                                       debug) {
  .debug( # nolint
    debug, "considering imputing missing thresholds finally by batch"
  )
  ind_with_missing_gates <- setdiff(
    gate_tbl$ind, cp_tbl$ind[!is.na(cp_tbl$cp_join_tg_orig)]
  )
  ind_with_missing_gates <- ind_with_missing_gates[
    !is.na(ind_with_missing_gates)
  ]
  if (length(ind_with_missing_gates) == 0L) { # nolint
    .debug(debug, "no missing thresholds finally by batch") # nolint
    return(cp_tbl)
  }
  .debug(debug, "imputing missing thresholds finally by batch") # nolint
  for (ind_curr in ind_with_missing_gates) {
    gate_tbl_ind <- gate_tbl |> dplyr::filter(ind == ind_curr) # nolint
    batch <- gate_tbl_ind |>
      dplyr::pull("batch")

    gate_tbl_batch_ind_vec <- gate_tbl |>
      dplyr::filter(batch == .env$batch) |> # nolint
      dplyr::pull("ind")
    gate_tbl_batch_ind_vec <- gate_tbl_batch_ind_vec[
      !is.na(gate_tbl_batch_ind_vec)
    ]

    cp_tbl_batch <- cp_tbl |>
      dplyr::filter(ind %in% gate_tbl_batch_ind_vec) # nolint

    dens_tbl_ind <- dens_tbl |>
      dplyr::filter(ind == ind_curr) # nolint

    # this means no samples had estimates for this batch
    cp_vec_imp <- cp_tbl_batch$cp_join_tg_orig[
      !is.na(cp_tbl_batch$cp_join_tg_orig)
    ]

    if (length(cp_vec_imp) == 0) next

    gate_impute <- .combine_cp( # nolint
      stats::setNames(cp_vec_imp, paste0("a", seq_along(cp_vec_imp))),
      gate_tbl$gate_combn[1]
    )[[1]][[1]]

    cp_tbl_add <- cp_tbl[1, ] |>
      tidyr::pivot_longer(cols = -grp) |> # nolint
      dplyr::mutate(value = NA_real_) |>
      tidyr::pivot_wider(id_cols = grp) |>
      dplyr::mutate(
        grp = ifelse(nrow(dens_tbl_ind) > 0, dens_tbl_ind$grp[1], NA),
        ind = ind_curr,
        cp_join_tg_orig = gate_impute
      )

    cp_tbl <- cp_tbl |> dplyr::filter(ind != ind_curr) # nolint

    cp_tbl <- cp_tbl |> dplyr::bind_rows(cp_tbl_add)
  }
  cp_tbl <- cp_tbl |> dplyr::arrange(ind) # nolint
}

.get_cp_cluster_plot_thresholds <- function(cp_tbl,
                                            params,
                                            marker_list,
                                            chnl,
                                            gs,
                                            ind_in_batch_lab_vec,
                                            ind_batch_list,
                                            ind_in_batch_uns,
                                            debug) {
  # diff by cp
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = cp_orig_quant_min - cp_join_lse_orig_mean) |> # nolint
    dplyr::filter(diff == max(diff)) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)

  # diff by prop
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = prop_bs_cp_diff - prop_bs_cp_diff) |> # nolint
    dplyr::filter(diff == max(abs(diff))) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)

  cp_vec <- stats::setNames(sel_tbl$cp_join_lse_orig, sel_tbl$grp)

  # tg
  # -------------------------

  # diff by cp
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = cp_orig_quant_min - cp_join_tg_orig_mean) |> # nolint
    dplyr::filter(diff == max(diff)) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)

  # diff by prop
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = prop_bs_cp_diff - prop_bs_cp_diff) |> # nolint
    dplyr::filter(diff == max(abs(diff))) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)

  # diff by cp
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = cp_orig_quant_min - cp_join_lse_orig_mean_tg) |> # nolint
    dplyr::filter(diff == max(diff)) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)

  # diff by prop
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = prop_bs_cp_diff - prop_bs_cp_diff) |> # nolint
    dplyr::filter(diff == max(abs(diff))) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)


  cp_vec <- stats::setNames(sel_tbl$cp_join_lse_orig_mean_tg, sel_tbl$grp)

  plot_list <- purrr::map(sel_tbl$ind, function(i) {
    ind_batch <- which(
      purrr::map_lgl(
        params$ind_batch_list,
        function(ind_batch) i %in% ind_batch
      )
    )
    ind_batch <- params$ind_batch_list[[ind_batch]]
    ex_tbl_all <- .get_ex_list( # nolint
      data = gs, ind_batch = ind_batch,
      ind_in_batch_gate = 1:5, ind_in_batch_uns = 5,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec,
      pop = "root",
      cut = purrr::map_chr(marker_list, function(x) x$cut), high = NULL,
      data_name = purrr::map
    ) |>
      dplyr::bind_rows()
    ex_uns <- ex_tbl_all |> dplyr::filter(stim == "uns") # nolint
    ex_stim <- ex_tbl_all |> dplyr::filter(.data$ind == .env$i) # nolint
    cut_alt <- ifelse(chnl == "Ho165Di", "Nd146Di", "Ho165Di")
    ex_uns <- ex_uns[ex_uns[[chnl]] > min(ex_uns[[chnl]]) |
      ex_uns[[cut_alt]] > min(ex_uns[[cut_alt]]), ]
    ex_stim <- ex_stim[ex_stim[[chnl]] > min(ex_stim[[chnl]]) |
      ex_stim[[cut_alt]] > min(ex_stim[[cut_alt]]), ]
    cn_vec_uns <- colnames(ex_uns)
    chnl_ind <- which(cn_vec_uns == chnl)
    cn_vec_uns[chnl_ind] <- "x"
    cut_alt_ind <- which(cn_vec_uns == cut_alt)
    cn_vec_uns[cut_alt_ind] <- "y"
    colnames(ex_uns) <- cn_vec_uns

    cn_vec_stim <- colnames(ex_stim)
    chnl_ind <- which(cn_vec_stim == chnl)
    cn_vec_stim[chnl_ind] <- "x"
    cut_alt_ind <- which(cn_vec_stim == cut_alt)
    cn_vec_stim[cut_alt_ind] <- "y"
    colnames(ex_stim) <- cn_vec_stim

    ex_plot <- ex_uns |>
      dplyr::bind_rows(ex_stim)
    ggplot(ex_plot, aes(x = x, y = y)) + # nolint
      geom_hex() + # nolint
      scale_fill_viridis_c(trans = "log10") + # nolint
      facet_wrap(~stim, ncol = 1) + # nolint
      geom_vline( # nolint
        xintercept = cp_vec[sel_tbl$grp[sel_tbl$ind == i]],
        col = "red", size = 2
      ) +
      coord_cartesian( # nolint
        xlim = c(0, 7),
        ylim = c(0, 7)
      ) +
      labs(x = chnl, y = cut_alt) # nolint
  })

  cowplot::plot_grid(
    plotlist = plot_list, ncol = 2,
    labels = names(cp_vec)
  )
}
