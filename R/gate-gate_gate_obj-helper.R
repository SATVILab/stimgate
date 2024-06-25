.get_pop_man_vec <- function(pop_man_sub,
                             pop_man_match_exact,
                             pop_gate) {
  if (!is.null(pop_man_sub)) {
    if (pop_man_match_exact) {
      pop_man_vec <- paste0(pop_gate, "/", pop_man_sub)
    } else {
      pop_man_vec <- pop_man_sub
    }
  } else {
    pop_man_vec <- NULL
  }
  pop_man_vec
}
.get_gate_obj_delete_old_gates <- function(bias_uns,
                                           data_name) {
  for (bias in bias_uns) {
    dir_save <- file.path(
      tempdir(),
      data_name, paste0("cp_locb", bias, "_plots")
    )
    if (dir.exists(dir_save)) unlink(dir_save, recursive = TRUE)
  }
}

.get_gate_obj_ind_uns_vec_get <- function(ind_batch_list, ind_in_batch_uns) {
  purrr::map(ind_batch_list, function(x) x[ind_in_batch_uns])
}

# Get gates for each sample within each batch
.get_gate_obj_pre_adj_gates_gate <- function(ind_batch_list,
                                             data,
                                             ind_in_batch_gate,
                                             ind_in_batch_uns,
                                             ind_in_batch_lab_vec,
                                             pop_gate,
                                             cut,
                                             high,
                                             gate_combn,
                                             pop_man,
                                             pop_man_match_exact,
                                             tol,
                                             data_name,
                                             fdr,
                                             noise_sd,
                                             bias_uns,
                                             bw_min,
                                             cp_min,
                                             boot_n,
                                             boot_sd,
                                             min_cell,
                                             params,
                                             plot,
                                             debug) {
  print("getting pre-adjustment gates")
  purrr::map_df(seq_along(ind_batch_list), function(i) {
    .debug(debug, "ind_batch_list", i) # nolint

    # print progress
    if (i %% 50 == 0 || i == length(ind_batch_list)) {
      print(paste0("batch ", i, " of ", length(ind_batch_list)))
    }
    .get_gate_batch_boot( # nolint
      data = data,
      ind_batch = ind_batch_list[[i]],
      ind_in_batch_gate = ind_in_batch_gate,
      ind_in_batch_uns = ind_in_batch_uns,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec,
      pop_gate = pop_gate,
      cut = cut,
      high = high,
      gate_combn = gate_combn,
      pop_man = pop_man,
      pop_man_match_exact = pop_man_match_exact,
      tol = tol,
      data_name = data_name,
      fdr = fdr,
      noise_sd = noise_sd,
      bias_uns = bias_uns,
      bw_min = bw_min,
      cp_min = cp_min,
      boot_n = boot_n,
      boot_sd = boot_sd,
      min_cell = min_cell,
      params = params,
      plot = plot,
      debug = debug
    ) |>
      dplyr::mutate(batch = i) |>
      dplyr::select(
        gate_name, gate_type, gate_combn, # nolint
        batch, ind, gate, everything() # nolint
      )
  })
}

.get_gate_obj_get_adj_gates <- function(gate_tbl,
                                        gate_tbl_params,
                                        tol_ctrl,
                                        tol_gate,
                                        gate_quant,
                                        pop_sub,
                                        data,
                                        params,
                                        cut,
                                        bw_min,
                                        path_project,
                                        debug,
                                        ind_batch_list,
                                        ind_in_batch_lab_vec,
                                        pop_gate,
                                        data_name,
                                        ind_in_batch_uns) {
  if (is.null(gate_tbl_params)) {
    .get_gate_obj_get_adj_gates_all( # nolint
      tol_ctrl = tol_ctrl,
      tol_gate = tol_gate,
      gate_tbl = gate_tbl,
      pop_sub = pop_sub,
      params = params,
      cut = cut,
      data = data,
      bw_min = bw_min,
      path_project = path_project,
      gate_quant = gate_quant,
      debug = debug,
      ind_batch_list = ind_batch_list,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec,
      pop_gate = pop_gate,
      data_name = data_name,
      ind_in_batch_uns = ind_in_batch_uns
    )
  } else {
    .get_gate_obj_gate_adj_gates_single(
      gate_tbl = gate_tbl,
      gate_tbl_params = gate_tbl_params,
      params = params,
      cut = cut,
      data = data,
      calc_cyt_pos_gates = TRUE,
      path_project = path_project,
      debug = debug,
      ind_batch_list = ind_batch_list,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec,
      pop_gate = pop_gate,
      data_name = data_name,
      ind_in_batch_uns = ind_in_batch_uns
    )
  }
}

.get_gate_obj_get_adj_gates_all <- function(tol_ctrl,
                                            tol_gate,
                                            gate_tbl,
                                            pop_sub,
                                            params,
                                            cut,
                                            data,
                                            bw_min,
                                            path_project,
                                            gate_quant,
                                            debug,
                                            ind_batch_list,
                                            ind_in_batch_lab_vec,
                                            pop_gate,
                                            data_name,
                                            ind_in_batch_uns) {
  if (!is.null(tol_ctrl)) {
    gate_tbl_ctrl <- gate_tbl |>
      dplyr::filter(gate_use == "ctrl") # nolint
  } else {
    gate_tbl_ctrl <- NULL
  }
  if (!is.null(tol_gate)) {
    gate_tbl_tg_gate <- gate_tbl |> dplyr::filter(gate_use == "tg_clust") # nolint
  } else {
    gate_tbl_tg_gate <- NULL
  }

  gate_tbl <- gate_tbl |> dplyr::filter(gate_use == "gate") # nolint
  gate_tbl <- gate_tbl |> dplyr::select(-gate_use) # nolint

  # =========================
  # Cluster-based gating
  # =========================

  if (!is.null(tol_gate)) {
    path_dir_stats <- .get_gate_stats( # nolint
      params = params,
      gate_tbl = gate_tbl |> dplyr::mutate(chnl = cut),
      gate_name = NULL,
      chnl = cut,
      data = data,
      filter_other_cyt_pos = FALSE,
      combn = FALSE,
      gate_type_single_pos_calc = "base",
      path_project = path_project,
      debug = debug,
      ind_batch_list = ind_batch_list,
      ind_in_batch_lab = ind_in_batch_lab_vec,
      pop_gate = pop_gate,
      data_name = data_name,
      ind_in_batch_uns = ind_in_batch_uns
    )
    gate_stats_tbl <- path_dir_stats |> .read_gate_stats() # nolint

    gate_tbl_cluster <- purrr::map_df(
      unique(gate_tbl$gate_name), function(gn) {
        gate_tbl_cluster <- .get_cp_cluster( # nolint
          gs = params$data,
          gate_tbl = gate_tbl |>
            dplyr::filter(gate_name == gn), # nolint
          gate_stats_tbl = gate_stats_tbl |>
            dplyr::filter(gate_name == gn),
          gate_tbl_ctrl = gate_tbl_tg_gate,
          chnl = cut,
          bw = bw_min,
          control = list(),
          filter_other_cyt_pos = FALSE,
          params = params,
          debug = debug
        )

        gate_tbl_cluster |>
          dplyr::select(ind, cp_join_tg_orig) |> # nolint
          dplyr::rename(gate = cp_join_tg_orig) |>
          dplyr::left_join(
            gate_tbl |>
              dplyr::filter(gate_name == gn) |> # nolint
              dplyr::select(
                gate_name, gate_type, gate_combn, # nolint
                boot, boot_ind, batch, ind # nolint
              ),
            by = c("ind")
          ) |>
          dplyr::select(
            gate_name, gate_type, gate_combn,
            batch, ind, gate, boot, boot_ind # nolint
          ) |>
          dplyr::mutate(
            gate_combn = paste0(gate_combn, "_clust"),
            gate_name = paste0(gate_type, "_", gate_combn)
          )
      }
    )
  }

  # =========================
  # Tail-gate controlled gating
  # =========================

  if (!is.null(tol_ctrl)) {
    if (is.null(tol_gate)) {
      path_dir_stats <- .get_gate_stats( # nolint
        params = params,
        gate_tbl = gate_tbl,
        pop_sub = pop_sub,
        data = data,
        path_project = path_project,
        debug = debug,
        ind_batch_list = ind_batch_list,
        ind_in_batch_lab = ind_in_batch_lab_vec,
        pop_gate = pop_gate,
        data_name = data_name,
        ind_in_batch_uns = ind_in_batch_uns
      )
      gate_stats_tbl <- path_dir_stats |> .read_gate_stats() # nolint
    }

    gate_tbl_2 <- .get_cp_adj_tbl( # nolint
      gate_stats_tbl = gate_stats_tbl,
      gate_quant = gate_quant,
      gate_tbl_ctrl = gate_tbl_ctrl
    )

    gate_tbl_2 <- gate_tbl_2 |>
      dplyr::select(-c(batch, boot)) # nolint
    gate_tbl_2 <- gate_tbl_2 |>
      dplyr::left_join(
        gate_tbl |>
          dplyr::select(batch, ind) |> # nolint
          dplyr::group_by(batch, ind) |>
          dplyr::slice(1) |>
          dplyr::ungroup(),
        by = c("ind")
      )

    gate_tbl <- gate_tbl_2
  }

  if (!is.null(tol_gate)) {
    gate_tbl <- gate_tbl |>
      dplyr::bind_rows(gate_tbl_cluster)
  }

  # Output
  # ------------------

  # output
  list(
    params = params,
    gate_tbl = gate_tbl
  )
}

.get_gate_obj_gate_adj_gates_single <- function(gate_tbl,
                                                params,
                                                gate_tbl_params,
                                                cut,
                                                data,
                                                calc_cyt_pos_gates,
                                                path_project,
                                                debug,
                                                ind_batch_list,
                                                ind_in_batch_lab_vec,
                                                pop_gate,
                                                data_name,
                                                ind_in_batch_uns) {
  # get gates
  gate_tbl_single <- gate_tbl

  # merge
  gate_tbl <- .get_gate_obj_gate_adj_gates_single_merge( # nolint
    gate_tbl_single = gate_tbl_single,
    gate_tbl_params = gate_tbl_params,
    cut = cut
  )

  # get stats table (if needed)
  gate_stats_tbl <- .get_gate_obj_gate_adj_gates_single_stats_tbl_get(
    gate_tbl = gate_tbl,
    params = params,
    cut = cut,
    data = data,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    path_project = path_project,
    debug = debug,
    ind_batch_list = ind_batch_list,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    pop_gate = pop_gate,
    data_name = data_name,
    ind_in_batch_uns = ind_in_batch_uns
  )

  gate_tbl_out <- .get_gate_obj_gate_adj_gates_single_out_get(
    gate_tbl = gate_tbl,
    gate_stats_tbl = gate_stats_tbl,
    gate_tbl_single = gate_tbl_single,
    params = params,
    debug = debug
  )

  list(
    params = params,
    gate_tbl = gate_tbl_out
  )
}

.get_gate_obj_gate_adj_gates_single_merge <- function(gate_tbl_single,
                                                      gate_tbl_params,
                                                      cut) {
  gate_tbl_params |>
    dplyr::left_join(
      gate_tbl_single |>
        dplyr::mutate(chnl = cut) |>
        dplyr::rename(gate_single = gate) |> # nolint
        dplyr::filter(gate_use == "gate") |> # nolint
        dplyr::select(
          chnl, gate_name, ind, gate_single # nolint
        ),
      by = c("chnl", "ind", "gate_name")
    ) |>
    dplyr::mutate(
      gate_type = purrr::map_chr(
        gate_name,
        function(gn) stringr::str_split(gn, "_")[[1]][1]
      ),
      gate_combn = gate_name |>
        stringr::str_remove("_adj") |>
        stringr::str_remove("_clust") |>
        stringr::str_remove(gate_type) |> # nolint
        stringr::str_remove("_")
    ) |>
    dplyr::select(
      chnl, marker, gate_name, # nolint
      gate_type, gate_combn, everything() # nolint
    )
}

.get_gate_obj_gate_adj_gates_single_stats_tbl_get <- function(gate_tbl,
                                                              params,
                                                              cut,
                                                              data,
                                                              calc_cyt_pos_gates, # nolint
                                                              path_project,
                                                              debug,
                                                              ind_batch_list,
                                                              ind_in_batch_lab_vec, # nolint
                                                              pop_gate,
                                                              data_name,
                                                              ind_in_batch_uns) { # nolint
  gate_name_vec <- unique(gate_tbl$gate_name)
  if (!.get_gate_obj_gate_adj_gates_single_stats_tbl_get_check(gate_name_vec)) {
    return(NULL)
  }


  gate_name_vec_clust <- gate_name_vec[
    stringr::str_detect(gate_name_vec, "_clust")
  ]
  gate_name_vec_adj <- gate_name_vec[
    stringr::str_detect(gate_name_vec, "_adj")
  ]

  .get_gate_stats( # nolint
    params = params,
    gate_tbl = gate_tbl |>
      dplyr::filter(
        gate_name %in% c(gate_name_vec_clust, gate_name_vec_adj) # nolint
      ),
    gate_name = NULL,
    chnl = cut,
    filter_other_cyt_pos = TRUE,
    gate_type_cyt_pos_filter = ifelse(calc_cyt_pos_gates, "cyt", "base"),
    data = data,
    gate_type_single_pos_filter = "base",
    gate_type_single_pos_calc = "base",
    combn = FALSE,
    path_project = path_project,
    debug = debug,
    ind_batch_list = ind_batch_list,
    ind_in_batch_lab = ind_in_batch_lab_vec,
    pop_gate = pop_gate,
    data_name = data_name,
    ind_in_batch_uns = ind_in_batch_uns
  )
}

.get_gate_obj_gate_adj_gates_single_stats_tbl_get_check <- function(gate_name_vec) { # nolint
  any_clust_ind <- any(
    purrr::map_lgl(
      gate_name_vec,
      function(gn) stringr::str_detect(gn, "_clust")
    )
  )
  any_adj_ind <- any(
    purrr::map_lgl(
      gate_name_vec,
      function(gn) stringr::str_detect(gn, "_adj")
    )
  )
  any_clust_ind || any_adj_ind
}

.get_gate_obj_gate_adj_gates_single_out_get <- function(gate_tbl,
                                                        gate_stats_tbl,
                                                        gate_tbl_single,
                                                        params,
                                                        debug) {
  # get tail-gate gates
  gate_tbl_ctrl_clust <- gate_tbl_single |>
    dplyr::filter(gate_use == "tg_clust") # nolint

  # get control gates
  gate_tbl_ctrl_ctrl <- gate_tbl_single |> # nolint
    dplyr::filter(gate_use == "ctrl") # nolint

  # get gate names
  gate_name_vec <- unique(gate_tbl$gate_name)
  purrr::map_df(gate_name_vec, function(gn) {
    .get_gate_obj_gate_adj_gates_single_out_get_gn(
      debug = debug,
      gn = gn,
      gate_tbl = gate_tbl,
      gate_stats_tbl = gate_stats_tbl,
      gate_tbl_ctrl_clust = gate_tbl_ctrl_clust,
      gate_tbl_ctrl_ctrl = gate_tbl_ctrl_ctrl,
      gate_tbl_single_gn = gate_tbl_single |>
        dplyr::filter(gate_name == gn), # nolint
      params = params
    )
  }) |>
    dplyr::mutate(gate_single = pmax(gate, gate_single)) # nolint
}

.get_gate_obj_gate_adj_gates_single_out_get_gn <- function(debug,
                                                           gn,
                                                           gate_tbl,
                                                           gate_stats_tbl,
                                                           gate_tbl_ctrl_clust,
                                                           gate_tbl_ctrl_ctrl,
                                                           gate_tbl_single_gn,
                                                           params) {
  .debug(debug, "gate_name_vec", gn) # nolint
  gate_tbl_gn <- gate_tbl |>
    dplyr::filter(gate_name == gn) # nolint

  adj_ind <- stringr::str_detect(gn, "_adj")

  clust_ind <- stringr::str_detect(gn, "_clust")
  if (!clust_ind && !adj_ind) {
    return(gate_tbl_gn |>
      dplyr::filter(chnl == params$cut) |> # nolint
      dplyr::select(
        chnl, marker, gate_name, gate_type, # nolint
        gate_combn, batch, ind, gate, gate_cyt, gate_single # nolint
      ))
  }

  # =========================
  # Tail-gate controlled gating
  # =========================

  if (adj_ind) {
    gate_stats_tbl_gn <- gate_stats_tbl |>
      dplyr::filter(gate_name == gn) # nolint
    gate_tbl_ctrl_ctrl <- gate_tbl_ctrl_ctrl |>
      dplyr::filter(gate_name == gn) # nolint

    gate_tbl_gn_2 <- .get_cp_adj_tbl( # nolint
      gate_stats_tbl = gate_stats_tbl_gn,
      gate_quant = params$gate_quant,
      gate_tbl_ctrl = gate_tbl_ctrl_ctrl
    )
    gate_tbl_gn_2 <- gate_tbl_gn_2 |> dplyr::select(-c(batch, boot)) # nolint
    gate_tbl_gn_2 <- gate_tbl_gn_2 |>
      dplyr::left_join(
        gate_tbl_single_gn |>
          dplyr::select(batch, ind) |> # nolint
          dplyr::group_by(batch, ind) |>
          dplyr::slice(1) |>
          dplyr::ungroup(),
        by = c("ind")
      ) |>
      dplyr::rename(gate_single = gate) # nolint

    gate_tbl_single_gn <- gate_tbl_gn_2
    return(gate_tbl_single_gn)
  }

  # =========================
  # Cluster-based gating
  # =========================

  if (clust_ind) {
    gate_stats_tbl_gn <- gate_stats_tbl |>
      dplyr::filter(gate_name == gn) # nolint
    gate_tbl_ctrl_clust_gn <- gate_tbl_ctrl_clust |>
      dplyr::filter(gate_name == gn) # nolint

    gate_tbl_cluster_gn <- .get_cp_cluster( # nolint
      gs = params$data,
      gate_tbl = gate_tbl_gn,
      gate_stats_tbl = gate_stats_tbl_gn,
      gate_tbl_ctrl = gate_tbl_ctrl_clust_gn,
      chnl = params$cut,
      bw = params$bw_min,
      control = list(),
      filter_other_cyt_pos = TRUE,
      params = params,
      debug = debug
    )

    gate_tbl_cluster_gn <- gate_tbl_cluster_gn |>
      dplyr::select(ind, cp_join_tg_orig) |> # nolint
      dplyr::rename(gate_single = cp_join_tg_orig) |>
      dplyr::left_join(
        gate_tbl |>
          dplyr::filter(
            gate_name == gn, # nolint
            chnl == params$cut # nolint
          ) |>
          dplyr::select(
            gate_name, gate_type, gate_combn, # nolint
            batch, ind, gate, gate_cyt # nolint
          ),
        by = c("ind")
      ) |>
      dplyr::mutate(
        chnl = params$cut,
        marker = params$chnl_lab[params$cut]
      ) |>
      dplyr::select(
        chnl, marker, gate_name, gate_type, gate_combn, # nolint
        batch, ind, gate, gate_cyt, gate_single # nolint
      ) |>
      dplyr::mutate(
        gate_combn = paste0(gate_combn, "_clust"),
        gate_name = paste0(gate_type, "_", gate_combn)
      )
  }

  gate_tbl_cluster_gn
}
