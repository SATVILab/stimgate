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
                                                  .debug = FALSE) {
  .debug_msg(.debug, "Updating gate statistics table") # nolint
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

.get_cp_cluster_prop_bs_by_cp_tbl_obj <- function(.data,
                                                  gate_tbl,
                                                  ind_batch_list,
                                                  chnl_cut,
                                                  pop_gate,
                                                  calc_cyt_pos_gates,
                                                  cp_min,
                                                  max_cp,
                                                  gate_stats_tbl,
                                                  filter_other_cyt_pos,
                                                  .debug) {
  .debug_msg(.debug, "Getting prop_bs_by_cp_tbl object") # nolint
  # statistics
  # ----------------

  data_list_obj <- .get_prop_bs_by_cp_tbl_data_list(
    .data = .data,
    gate_tbl = gate_tbl,
    ind_batch_list = ind_batch_list,
    chnl_cut,
    pop_gate = pop_gate,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    max_cp = max_cp,
    filter_other_cyt_pos = filter_other_cyt_pos,
    cp_min = cp_min,
    .debug = .debug
  )

  prop_bs_by_cp_tbl <- .get_prop_bs_by_cp_tbl_actual(
    data_list = data_list_obj[["data_list"]],
    cp_min = cp_min,
    max_cp = max_cp,
    gate_stats_tbl = gate_stats_tbl,
    .debug = .debug,
    ind_batch_list = ind_batch_list
  )

  list(
    "prop_bs_by_cp_tbl" = prop_bs_by_cp_tbl,
    "expr_max" = data_list_obj[["expr_max"]],
    "expr_min" = data_list_obj[["expr_min"]]
  )
}



.get_prop_bs_by_cp_tbl_data_list <- function(.data,
                                             gate_tbl,
                                             ind_batch_list,
                                             chnl_cut,
                                             pop_gate,
                                             calc_cyt_pos_gates,
                                             max_cp,
                                             filter_other_cyt_pos,
                                             cp_min,
                                             .debug) {
  .debug_msg(.debug, "Getting .data list") # nolint
  data_list <- .get_prop_bs_by_cp_tbl_data_list_init(
    ind_batch_list = ind_batch_list, .data = .data, pop_gate = pop_gate,
    chnl_cut = chnl_cut, filter_other_cyt_pos = filter_other_cyt_pos,
    calc_cyt_pos_gates = calc_cyt_pos_gates, cp_min = cp_min, .debug = .debug
  )
  .get_prop_bs_by_cp_tbl_data_list_final(data_list, max_cp)
}

.get_prop_bs_by_cp_tbl_data_list_init <- function(ind_batch_list,
                                                  .data,
                                                  pop_gate,
                                                  chnl_cut,
                                                  filter_other_cyt_pos,
                                                  calc_cyt_pos_gates,
                                                  cp_min,
                                                  .debug) {
  purrr::map(seq_along(ind_batch_list), function(i) {
    ind_batch <- ind_batch_list[[i]]
    ex_list <- .get_ex_list( # nolint
      .data = .data, ind_batch = ind_batch, pop = pop_gate,
      chnl_cut = chnl_cut, batch = names(ind_batch_list)[i]
    )

    # get min and max expression across batch
    # (used downstream, and only cells with
    # expr above min are used for clustering)
    min_max_vec <- .get_prop_bs_by_cp_tbl_data_list_minmax(ex_list)

    # filter to yield cells negative for all cytokine combinations
    # except possible this cytokine single-positive
    ex_list_filter <- .get_prop_bs_by_cp_tbl_data_list_filter_cyt_pos(
      filter_other_cyt_pos = filter_other_cyt_pos, gate_tbl = gate_tbl,
      ex_list = ex_list, calc_cyt_pos_gates = calc_cyt_pos_gates
    )

    # this is just to save memory
    # for when calculating the per-threshold performance
    out_tbl <- .get_prop_bs_by_cp_tbl_data_list_filter_above_min(
      ex_list_filter, cp_min
    )

    list(
      "out_tbl" = out_tbl,
      "expr_min" = min_max_vec[[1]],
      "expr_max" = min_max_vec[[2]]
    )
  })
}

.get_prop_bs_by_cp_tbl_data_list_minmax <- function(ex_list) {
  # range of expressions
  # (used later when calculating how "flat" a derivative is)
  expr_range_tbl <- purrr::map_df(
    seq_along(ex_list),
    function(i) {
      ex <- ex_list[[i]]
      if (nrow(ex) <= 5) {
        return(NULL)
      }
      quant_vec <- quantile(.get_cut(ex), c(0.0025, 0.999))
      tibble::tibble(
        lb = quant_vec[[1]],
        ub = 3 * quant_vec[[2]]
      )
    }
  )

  expr_min <- quantile(expr_range_tbl[["lb"]], 0.0025)
  expr_max <- max(expr_range_tbl[["ub"]])

  c("min" = expr_min, "max" = expr_max)
}

.get_prop_bs_by_cp_tbl_data_list_filter_cyt_pos <- function(filter_other_cyt_pos,
                                                            ex_list,
                                                            gate_tbl,
                                                            calc_cyt_pos_gates) {
  if (!filter_other_cyt_pos) {
    return(ex_list)
  }
  purrr::map(seq_along(ex_list), function(i) {
    if (i == length(ex_list)) {
      return(ex_list[[i]])
    }
    gate_tbl_ind <- gate_tbl |>
      dplyr::filter(ind == attr(ex_list[[i]], "ind")) # nolint

    pos_ind_vec_but_single_pos_curr <-
      .get_pos_ind_but_single_pos_for_one_cyt( # nolint
        ex = ex_list[[i]],
        gate_tbl = gate_tbl_ind,
        chnl_single_exc = attr(ex_list[[i]], "chnl_cut"),
        chnl = NULL,
        gate_type_cyt_pos = ifelse(calc_cyt_pos_gates,
          "cyt", "base"
        ),
        gate_type_single_pos = "base"
      )
    ex_list[[i]][!pos_ind_vec_but_single_pos_curr, , drop = FALSE]
  }) |>
    stats::setNames(names(ex_list))
}

.get_prop_bs_by_cp_tbl_data_list_filter_above_min <- function(ex_list_filter,
                                                              cp_min) {
  ex_list_filter |>
    purrr::map(function(x) {
      attr(x, "n_cell") <- nrow(x)
      x_out <- x |>
        dplyr::filter(.get_cut(x) >= min(.env$cp_min, max(.get_cut(x)))) # nolint
      if (nrow(x_out) == 0) {
        all_cols <- colnames(x)
        batch_idx <- which(all_cols == "batch")
        stim_idx <- which(all_cols == "stim")
        sel_idx <- seq(batch_idx, stim_idx)
        x_out <- x[1, sel_idx, drop = FALSE]
        x_add <- x[1, setdiff(seq_along(x), sel_idx)]
        x_add[1, ] <- NA
        x_out <- x_out |>
          dplyr::bind_cols(x_add)
      }
      x_out
    })
}

.get_prop_bs_by_cp_tbl_data_list_final <- function(data_list, max_cp) {
  expr_min_vec <- vapply(data_list, function(x) x$expr_min, numeric(1))
  expr_max_vec <- vapply(data_list, function(x) x$expr_max, numeric(1))
  expr_min <- min(expr_min_vec, na.rm = TRUE)
  expr_max <- max(
    max(expr_max_vec, na.rm = TRUE),
    max_cp + 0.2 * (max(expr_max_vec, na.rm = TRUE) - expr_min)
  )
  data_list <- lapply(data_list, function(x) x$out_tbl) |>
    purrr::flatten()
  list(
    "data_list" = data_list, "expr_min" = expr_min, "expr_max" = expr_max
  )
}


.get_prop_bs_by_cp_tbl_actual <- function(data_list,
                                          cp_min,
                                          max_cp,
                                          gate_stats_tbl,
                                          .debug,
                                          ind_batch_list,
                                          ind_in_batch_uns) {
  .debug_msg(.debug, "Getting prop_bs_by_cp_tbl") # nolint
  cp_par_list <- .get_prop_bs_by_cp_tbl_actual_prep(cp_min, max_cp)
  purrr::map(seq_along(data_list), function(i) {
    .get_prop_bs_by_cp_tbl_actual_ind(
      .debug = .debug,
      i = i,
      data_list = data_list,
      cp_par_list = cp_par_list,
      gate_stats_tbl = gate_stats_tbl,
      ind_batch_list = ind_batch_list
    )
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()
}

.get_prop_bs_by_cp_tbl_actual_ind <- function(.debug,
                                              i,
                                              cp_par_list,
                                              gate_stats_tbl,
                                              ind_batch_list,
                                              data_list) {
  .get_prop_bs_by_cp_tbl_actual_progress(.debug, i, data_list)
  ex_list <- .get_prop_bs_by_cp_tbl_actual_ex_get(
    data_list, i, ind_batch_list
  )
  if (is.null(ex_list)) {
    return(NULL)
  }
  .get_prop_bs_by_cp_tbl_ind(
    ex_stim = ex_list$stim,
    ex_uns = ex_list$uns,
    cp_seq = cp_par_list[["seq"]],
    gate_stats_tbl = gate_stats_tbl,
    .debug = .debug
  )
}

.get_prop_bs_by_cp_tbl_actual_prep <- function(cp_min, cp_max) {
  cp_range <- c(cp_min, cp_max)
  cp_seq_vec <- seq(cp_range[1], cp_range[2], length.out = 1e2)
  list("range" = cp_range, "seq" = cp_seq_vec)
}

.get_prop_bs_by_cp_tbl_ind <- function(ex_stim,
                                       ex_uns,
                                       cp_seq,
                                       gate_stats_tbl,
                                       .debug) {
  par_list <- .get_prop_bs_by_cp_tbl_ind_prep(
    gate_stats_tbl, ex_stim, ex_uns, cp_seq
  )

  .get_prop_bs_by_cp_tbl_ind_init(ex_stim, ex_uns, par_list, cp_seq) |>
    .get_prop_bs_by_cp_tbl_ind_calc(
      attr(ex_stim, "n_cell"), attr(ex_uns, "n_cell")
      )
}


.get_prop_bs_by_cp_tbl_ind_prep <- function(gate_stats_tbl,
                                            ex_stim,
                                            ex_uns,
                                            cp_seq) {
  gate_stats_tbl_curr <- gate_stats_tbl |>
    dplyr::filter(.data$ind == attr(ex_stim, "ind")) # nolint
  count_stim_vec <- rep(NA, length(cp_seq))
  count_uns_vec <- rep(NA, length(cp_seq))
  for (i in seq_along(cp_seq)) {
    cp <- cp_seq[i]
    count_stim_vec[i] <- sum(.get_cut(ex_stim) > cp) # nolint
    count_uns_vec[i] <- sum(.get_cut(ex_uns) > cp) # nolint
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

.get_prop_bs_by_cp_tbl_ind_init <- function(ex_stim,
                                            ex_uns,
                                            par_list,
                                            cp_seq) {
  tibble::tibble(
    ind = attr(ex_stim, "ind"),
    prop_bs_orig = par_list[["bs_orig"]], prop_bs_sd = par_list[["bs_sd"]],
    cp = cp_seq, max_expr = max(.get_cut(ex_stim), .get_cut(ex_uns)), # nolint
    count_stim_cp = par_list[["count_stim"]],
    count_uns_cp = par_list[["count_uns"]]
  )
}

.get_prop_bs_by_cp_tbl_ind_calc <- function(.data, n_cell_stim, n_cell_uns) {
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

.get_prop_bs_by_cp_tbl_actual_progress <- function(.debug, i, data_list) {
  if (i %% 20 == 0) {
    .debug_msg(.debug, paste0("Processing ", i, " of ", length(data_list))) # nolint
  }
}

.get_prop_bs_by_cp_tbl_actual_progress <- function(.debug, i, data_list) {
  if (i %% 20 == 0) {
    .debug_msg(.debug, paste0("Processing ", i, " of ", length(data_list))) # nolint
  }
}

.get_prop_bs_by_cp_tbl_actual_ex_get <- function(data_list,
                                                 i,
                                                 ind_batch_list) {
  ex <- data_list[[i]]
  if (.get_prop_bs_by_cp_return_early_stim(ex)) {
    return(NULL)
  }
  ex_uns <- data_list[[attr(ex, "ind_uns")]]

  if (.get_prop_bs_by_cp_return_early_uns(ex_uns)) {
    return(NULL)
  }

  list("stim" = ex, "uns" = ex_uns)
}

.get_prop_bs_by_cp_return_early_stim <- function(ex) {
  is.na(.get_cut(ex)[1]) || attr(ex, "is_uns")
}

.get_prop_bs_by_cp_return_early_uns <- function(ex) {
  is.na(.get_cut(ex)[1])
}

.get_cp_cluster_dens_tbl_get_batch_prep <- function(ind_batch) {
  .debug_msg(.debug, paste0("Processing batch ", ind_batch)) # nolint
}


.get_cp_cluster_dens_tbl_get <- function(ind_batch_list,
                                         .data,
                                         filter_other_cyt_pos,
                                         calc_cyt_pos_gates,
                                         chnl_cut,
                                         expr_min,
                                         expr_max,
                                         pop_gate,
                                         gate_tbl,
                                         control,
                                         bw,
                                         .debug) {
  .debug_msg(.debug, "Getting density table") # nolint
  min_threshold <- .get_cp_cluster_dens_tbl_get_min_threshold(
    gate_tbl, control
  )
  dens_tbl <- purrr::map_df(seq_along(ind_batch_list), function(i) {
    ind_batch <- ind_batch_list[[i]]
    ex_list <- .get_cp_cluster_dens_tbl_get_batch_prep_ex_list(
      .data = .data, ind_batch = ind_batch,
      pop_gate = pop_gate, chnl_cut,
      filter_other_cyt_pos = filter_other_cyt_pos,
      gate_tbl = gate_tbl, calc_cyt_pos_gates = calc_cyt_pos_gates,
      control = control, .debug = .debug, batch = names(ind_batch_list)[i]
    )

    purrr::map_df(ex_list, function(x) {
      .get_cp_cluster_dens_tbl_get_actual_ind(
        expr_vec = .get_cut(x), batch = attr(x, "batch"),
        ind = attr(x, "ind"), min_threshold = min_threshold, chnl_cut,
        expr_min = expr_min, expr_max = expr_max, bw = bw, .debug = .debug
      )
    })
  }) |>
    dplyr::filter(!is.na(x1)) # nolint
  not_all_na_vec_ind <- purrr::map_lgl(
    seq_len(ncol(dens_tbl)),
    function(i) !all(is.na(dens_tbl[[i]]))
  )
  dens_tbl[, not_all_na_vec_ind]
}



.get_cp_cluster_dens_tbl_get_actual_ind <- function(expr_vec,
                                                    batch,
                                                    ind,
                                                    min_threshold,
                                                    chnl_cut,
                                                    expr_min,
                                                    expr_max,
                                                    bw,
                                                    .debug) {
  if (.get_cp_cluster_dens_tbl_get_actual_ind_early_return_check(expr_vec)) { # nolint
    return(.get_cp_cluster_dens_tbl_get_actual_ind_early_return(
      batch, ind
    ))
  }
  dens <- density(
    expr_vec[expr_vec > min(expr_vec)],
    from = expr_min, to = expr_max, bw = bw
  )
  tibble::tibble(
    batch = batch[[1]],ind = ind[[1]],
    y = dens[["y"]], x = dens[["x"]],
    x_ind = paste0("x", seq_len(length(dens[["y"]])))
  ) |>
    # extract only left-hand size elements
    dplyr::filter(x <= min_threshold) |> # nolint
    dplyr::select(-x) |>
    dplyr::mutate(y = y / sum(y)) |> # normalise size
    tidyr::pivot_wider(names_from = x_ind, values_from = y) # nolint
}

.get_cp_cluster_dens_tbl_get_actual_ind_early_return_check <-
  function(expr_vec) {
    length(expr_vec[expr_vec > min(expr_vec)]) < 3
  }

.get_cp_cluster_dens_tbl_get_actual_ind_early_return <- function(batch,
                                                                 ind) {
  tibble::tibble(
    batch = batch[1], ind = ind[1],
    y = rep(NA, 512), x = paste0("x", seq.int(from = 1, to = 512))
  ) |>
    tidyr::pivot_wider(names_from = x, values_from = y) # nolint
}


.get_cp_cluster_dens_tbl_get_min_threshold <- function(gate_tbl, control) {
  # get a low gate value, such that
  # when we look at clustering samples together
  # we only use expression levels below this value.
  # the thinking is that this makes the clustering depend
  # on the part that's got most of the cells and exhibits
  # the batch effects, and isn't dependent on how many cells
  # responded.
  # could cause problems if >40% of cells respond,
  # but we can let people switch this off.
  min_threshold_gate_quant <- quantile(
    gate_tbl$gate, control$min_threshold_quant,
    na.rm = TRUE
  )
  control$min_threshold_frac * min_threshold_gate_quant
}

.get_cp_cluster_dens_tbl_get_batch_prep_ex_list <- function(.data,
                                                            ind_batch,
                                                            batch,
                                                            pop_gate,
                                                            chnl_cut,
                                                            filter_other_cyt_pos, # nolint
                                                            gate_tbl,
                                                            calc_cyt_pos_gates,
                                                            control,
                                                            .debug) {
  ex_list <- .get_ex_list( # nolint
    .data = .data, ind_batch = ind_batch,
    pop = pop_gate, chnl_cut = chnl_cut, batch = batch
  )

  if (!filter_other_cyt_pos) {
    return(ex_list[-length(ex_list)])
  }

  # filter to yield cells negative for all cytokine combinations
  # except possible this cytokine single-positive
  .get_cp_cluster_dens_tbl_get_batch_prep_ex_list_filter(
    .debug = .debug, ex_list = ex_list,
    chnl_cut, gate_tbl = gate_tbl, calc_cyt_pos_gates = calc_cyt_pos_gates
  )
}

.get_cp_cluster_dens_tbl_get_batch_prep_ex_list_filter <- function(.debug,
                                                                   ex_list,
                                                                   chnl_cut,
                                                                   gate_tbl,
                                                                   calc_cyt_pos_gates) { # nolint
  .debug_msg(.debug, "Filtering other cytokine positive cells") # nolint
  ex_list_filter <- purrr::map(seq_along(ex_list), function(i) {
    if (i == length(ex_list)) {
      return(ex_list[[i]])
    }
    .get_cp_cluster_dens_tbl_get_batch_prep_ex_list_filter_ind(
      ex_list[[i]], gate_tbl, chnl_cut, calc_cyt_pos_gates
    )
  }) |>
    stats::setNames(names(ex_list))

  ex_list_filter[-length(ex_list_filter)]
}

.get_cp_cluster_dens_tbl_get_batch_prep_ex_list_filter_ind <- function(ex_tbl,
                                                                       gate_tbl,
                                                                       chnl_cut,
                                                                       calc_cyt_pos_gates) { # nolint

  pos_ind_vec_but_single_pos_curr <-
    .get_pos_ind_but_single_pos_for_one_cyt( # nolint
      ex = ex_tbl,
      gate_tbl = gate_tbl[gate_tbl[["ind"]] == attr(ex_tbl, "ind"), ],
      chnl_single_exc = chnl_cut,
      chnl = NULL,
      gate_type_cyt_pos = if (calc_cyt_pos_gates) "cyt" else "base",
      gate_type_single_pos = "base"
    )
  ex_tbl[!pos_ind_vec_but_single_pos_curr, , drop = FALSE]
}

.get_cp_cluster_n_clus <- function(dens_tbl) {
  max_cluster <- min(6, nrow(dens_tbl) / 3) |>
    floor() |>
    max(1)
  if (max_cluster == 1L) {
    return(1L)
  }
  for (i in seq_len(ncol(dens_tbl))) {
    if (any(is.na(dens_tbl[[i]]))) {
      message(i)
    }
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
                                               .debug = FALSE) {
  data_mod <- .get_cp_cluster_data_mod_pre(
    prop_bs_by_cp_tbl = prop_bs_by_cp_tbl
  )
  purrr::map(unique(data_mod$grp), function(n_grp_curr) {
    .debug_msg(.debug, paste0("Processing cluster ", n_grp_curr)) # nolint
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
              .data = data_mod_curr_grp_not_na,
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
          cp = seq(min_cp_permitted, cp_range[2],length.out = 1e5)
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
        # we divide this into 1/1000
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
                                                   chnl_cut,
                                                   grp_ind_lab_vec,
                                                   .debug) {
  .debug_msg( # nolint
    .debug, "Getting quantiles of original gates per clustered observations" # nolint
  )
  if ("chnl" %in% names(gate_tbl)) {
    gate_tbl <- gate_tbl |>
      dplyr::filter(.data$chnl == .env$chnl_cut) # nolint
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
                                        .debug) {
  .debug_msg(.debug, "Getting cp_join") # nolint
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
                                            .debug) {
  .debug_msg(.debug, "Adding information to cp table") # nolint
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

.get_cp_cluster_cp_tbl_add_orig_quant_min <- function(cp_tbl, .debug) {
  .debug_msg(.debug, "Adding original and minimum quantile threshold") # nolint
  cp_tbl |>
    dplyr::mutate(
      cp_orig_quant_min = pmax(pmin(cp_orig, max_expr), gate_05) # nolint
    )
}
.get_cp_cluster_cp_join_lse_get <- function(cp_tbl, .debug) {
  .debug_msg(.debug, "Getting cp_join_lse") # nolint
  cp_tbl <- cp_tbl |>
    dplyr::group_by(ind) |> # nolint
    dplyr::mutate(
      lse_orig = prop_bs_cp_diff_sd < 0.01, # nolint
      cp_join_lse = min(cp[lse_orig & cp >= cp_join]), # nolint
      cp_join_lse_orig = pmin(cp_join_lse, cp_orig_quant_min) # nolint
    ) |>
    dplyr::ungroup()
}

.get_cp_cluster_cp_join_tg_get <- function(cp_tbl, .debug) {
  .debug_msg(.debug, "Getting cp_join_tg") # nolint
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

.get_cp_cluster_cp_lse_orig_mean <- function(cp_tbl, .debug) {
  .debug_msg(.debug, "Getting cp_join_lse_orig_mean") # nolint
  cp_tbl |>
    dplyr::mutate(
      cp_join_lse_orig_mean = pmin(
        cp_orig_quant_min, # nolint
        cp_join_lse_orig + # nolint
          (cp_orig_quant_min - cp_join_lse_orig) / 2
      )
    )
}

.get_cp_cluster_cp_join_tg_orig_mean <- function(cp_tbl, .debug) {
  .debug_msg(.debug, "Getting cp_join_tg_orig_mean") # nolint
  cp_tbl |>
    dplyr::mutate(
      cp_join_tg_orig_mean = pmin(
        cp_orig_quant_min, # nolint
        cp_join_tg_orig + (cp_orig_quant_min - cp_join_tg_orig) / 2 # nolint
      )
    )
}

.get_cp_cluster_cp_join_lse_orig_mean_tg <- function(cp_tbl, .debug) {
  .debug_msg(.debug, "Getting cp_join_lse_orig_mean_tg") # nolint
  cp_tbl |>
    dplyr::mutate(
      cp_join_lse_orig_mean_tg = pmin(
        cp_join_lse_orig_mean, cp_join_tg_orig # nolint
      )
    )
}

.get_cp_cluster_cp_filter_above_cp_join_lse_orig_mean_tg <-
  function(cp_tbl, .debug) {
    .debug_msg( # nolint
      .debug, "Filtering above cp_join_lse_orig_mean_tg"
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
                                                    chnl_cut,
                                                    gate_tbl,
                                                    dens_tbl,
                                                    .debug) {
  .debug_msg( # nolint
    .debug, "considering imputing missing thresholds by batch"
  ) # nolint
  ind_with_missing_gates <- setdiff(
    gate_tbl$ind, cp_tbl$ind[!is.na(cp_tbl$cp_join_tg_orig)]
  )
  ind_with_missing_gates <- ind_with_missing_gates[
    !is.na(ind_with_missing_gates)
  ]
  if (length(ind_with_missing_gates) == 0L) {
    .debug_msg(.debug, "no missing thresholds by batch") # nolint
    return(cp_tbl)
  }
  .debug_msg(.debug, "imputing missing thresholds by batch") # nolint
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
                                               chnl_cut,
                                               gate_tbl,
                                               dens_tbl,
                                               .debug) {
  .debug_msg(.debug, "considering imputing missing thresholds individually") # nolint
  ind_with_missing_gates <- setdiff(
    gate_tbl$ind, cp_tbl$ind[!is.na(cp_tbl$cp_join_tg_orig)]
  )
  ind_with_missing_gates <- ind_with_missing_gates[
    !is.na(ind_with_missing_gates)
  ]
  if (length(ind_with_missing_gates) == 0L) {
    .debug_msg(.debug, "no missing thresholds individually") # nolint
    return(cp_tbl)
  }
  .debug_msg(.debug, "imputing missing thresholds individually") # nolint
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
                                                 chnl_cut,
                                                 gate_tbl,
                                                 dens_tbl,
                                                 .debug) {
  .debug_msg(.debug, "considering imputing missing thresholds finally") # nolint
  ind_with_missing_gates <- setdiff(
    gate_tbl$ind, cp_tbl$ind[!is.na(cp_tbl$cp_join_tg_orig)]
  )
  ind_with_missing_gates <- ind_with_missing_gates[
    !is.na(ind_with_missing_gates)
  ]
  if (length(ind_with_missing_gates) == 0L) {
    .debug_msg(.debug, "no missing thresholds finally") # nolint
    return(cp_tbl)
  }
  .debug_msg(.debug, "imputing missing thresholds finally") # nolint
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
                                                       chnl_cut,
                                                       gate_tbl,
                                                       dens_tbl,
                                                       .debug) {
  .debug_msg( # nolint
    .debug, "considering imputing missing thresholds finally by batch"
  )
  ind_with_missing_gates <- setdiff(
    gate_tbl$ind, cp_tbl$ind[!is.na(cp_tbl$cp_join_tg_orig)]
  )
  ind_with_missing_gates <- ind_with_missing_gates[
    !is.na(ind_with_missing_gates)
  ]
  if (length(ind_with_missing_gates) == 0L) { # nolint
    .debug_msg(.debug, "no missing thresholds finally by batch") # nolint
    return(cp_tbl)
  }
  .debug_msg(.debug, "imputing missing thresholds finally by batch") # nolint
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

