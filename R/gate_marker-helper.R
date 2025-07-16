.gate_marker_delete_old_gates <- function() {
  dir_save <- file.path(tempdir(), "stimgate")
  if (!dir.exists(dir_save)) {
    return(invisible(FALSE))
  }
  unlink(dir_save, recursive = TRUE)
  invisible(TRUE)
}

# Get gates for each sample within each batch
.gate_marker_pre_adj_gates_gate <- function(ind_batch_list,
                                            .data,
                                            pop_gate,
                                            chnl_cut,
                                            gate_combn,
                                            tol_clust,
                                            noise_sd,
                                            bias_uns,
                                            exc_min,
                                            bw_min,
                                            cp_min,
                                            min_cell,
                                            params,
                                            .debug,
                                            path_project) {
  print("getting pre-adjustment gates")
  purrr::map_df(seq_along(ind_batch_list), function(i) {
    .debug_msg(.debug, "ind_batch_list", i) # nolint

    # print progress
    if (i %% 50 == 0 || i == length(ind_batch_list)) {
      print(paste0("batch ", i, " of ", length(ind_batch_list)))
    }
    .gate_batch( # nolint
      .data = .data,
      ind_batch = ind_batch_list[[i]],
      pop_gate = pop_gate,
      chnl_cut = chnl_cut,
      gate_combn = gate_combn,
      tol_clust = tol_clust,
      noise_sd = noise_sd,
      bias_uns = bias_uns,
      exc_min = exc_min,
      bw_min = bw_min,
      cp_min = cp_min,
      min_cell = min_cell,
      params = params,
      batch = names(ind_batch_list)[i],
      .debug = .debug,
      path_project = path_project
    ) |>
      dplyr::select(
        gate_name, gate_type, gate_combn, # nolint
        batch, ind, gate, everything() # nolint
      )
  })
}

.gate_marker_get_adj_gates <- function(gate_tbl,
                                       gate_tbl_params,
                                       tol_clust,
                                       gate_quant,
                                       .data,
                                       params,
                                       chnl_cut,
                                       bw_min,
                                       path_project,
                                       .debug,
                                       ind_batch_list,
                                       pop_gate) {
  if (is.null(gate_tbl_params)) {
    .gate_marker_get_adj_gates_all( # nolint
      tol_clust = tol_clust,
      gate_tbl = gate_tbl,
      params = params,
      chnl_cut,
      .data = .data,
      bw_min = bw_min,
      path_project = path_project,
      gate_quant = gate_quant,
      .debug = .debug,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate
    )
  } else {
    .gate_marker_gate_adj_gates_single(
      gate_tbl = gate_tbl,
      gate_tbl_params = gate_tbl_params,
      params = params,
      chnl_cut,
      .data = .data,
      calc_cyt_pos_gates = TRUE,
      path_project = path_project,
      .debug = .debug,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate
    )
  }
}

.gate_marker_get_adj_gates_all <- function(tol_clust,
                                           gate_tbl,
                                           params,
                                           chnl_cut,
                                           .data,
                                           bw_min,
                                           path_project,
                                           gate_quant,
                                           .debug,
                                           ind_batch_list,
                                           pop_gate) {
  # START HERE!!!
  if (!is.null(tol_clust)) {
    gate_tbl_tg_gate <- gate_tbl |>
      dplyr::filter(gate_use == "tg_clust") # nolint
  } else {
    gate_tbl_tg_gate <- NULL
  }

  gate_tbl <- gate_tbl |> dplyr::filter(gate_use == "gate") # nolint
  gate_tbl <- gate_tbl |> dplyr::select(-gate_use) # nolint

  # =========================
  # Cluster-based gating
  # =========================

  if (!is.null(tol_clust)) {
    path_dir_stats <- .get_stats( # nolint
      params = params,
      gate_tbl = gate_tbl |> dplyr::mutate(chnl = chnl_cut),
      gate_name = NULL,
      chnl = chnl_cut,
      .data = .data,
      filter_other_cyt_pos = FALSE,
      combn = FALSE,
      gate_type_single_pos_calc = "base",
      path_project = path_project,
      .debug = .debug,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate
    )
    gate_stats_tbl <- path_dir_stats |>
      .read_gate_stats() # nolint

    gate_tbl_cluster <- purrr::map_df(
      unique(gate_tbl$gate_name), function(gn) {
        gate_tbl_cluster <- .get_cp_cluster( # nolint
          .data = .data,
          gate_tbl = gate_tbl |>
            dplyr::filter(gate_name == gn), # nolint
          gate_stats_tbl = gate_stats_tbl |>
            dplyr::filter(gate_name == gn),
          gate_tbl_ctrl = gate_tbl_tg_gate,
          chnl = chnl_cut,
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
                batch, ind # nolint
              ),
            by = c("ind")
          ) |>
          dplyr::select(
            gate_name, gate_type, gate_combn,
            batch, ind, gate # nolint
          ) |>
          dplyr::mutate(
            gate_combn = paste0(gate_combn, "_clust"),
            gate_name = paste0(gate_type, "_", gate_combn)
          )
      }
    )
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


.gate_marker_gate_adj_gates_single <- function(gate_tbl,
                                               params,
                                               gate_tbl_params,
                                               chnl_cut,
                                               .data,
                                               calc_cyt_pos_gates,
                                               path_project,
                                               .debug,
                                               ind_batch_list,
                                               pop_gate) {
  # get gates
  gate_tbl_single <- gate_tbl

  # merge
  gate_tbl <- .gate_marker_gate_adj_gates_single_merge( # nolint
    gate_tbl_single = gate_tbl_single,
    gate_tbl_params = gate_tbl_params,
    chnl_cut = chnl_cut
  )

  # get stats table (if needed)
  gate_stats_tbl <- .gate_marker_gate_adj_gates_single_stats_tbl_get(
    gate_tbl = gate_tbl,
    params = params,
    chnl_cut,
    .data = .data,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    path_project = path_project,
    .debug = .debug,
    ind_batch_list = ind_batch_list,
    pop_gate = pop_gate
  )

  gate_tbl_out <- .gate_marker_gate_adj_gates_single_out_get(
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

.gate_marker_gate_adj_gates_single_merge <- function(gate_tbl_single,
                                                     gate_tbl_params,
                                                     chnl_cut) {
  gate_tbl_params |>
    dplyr::left_join(
      gate_tbl_single |>
        dplyr::mutate(chnl = chnl_cut) |>
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
    ) |>
    dplyr::mutate(
      gate_single = ifelse(is.na(gate_single), gate, gate_single) # nolint
    )
}

.gate_marker_gate_adj_gates_single_stats_tbl_get <- function(gate_tbl,
                                                             params,
                                                             chnl_cut,
                                                             .data,
                                                             calc_cyt_pos_gates, # nolint
                                                             path_project,
                                                             .debug,
                                                             ind_batch_list) { # nolint
  gate_name_vec <- unique(gate_tbl$gate_name)
  if (!.gate_marker_gate_adj_gates_single_stats_tbl_get_check(gate_name_vec)) {
    return(NULL)
  }


  gate_name_vec_clust <- gate_name_vec[
    stringr::str_detect(gate_name_vec, "_clust")
  ]
  gate_name_vec_adj <- gate_name_vec[
    stringr::str_detect(gate_name_vec, "_adj")
  ]

  .get_stats( # nolint
    params = params,
    gate_tbl = gate_tbl |>
      dplyr::filter(
        gate_name %in% c(gate_name_vec_clust, gate_name_vec_adj) # nolint
      ),
    gate_name = NULL,
    chnl = chnl_cut,
    filter_other_cyt_pos = TRUE,
    gate_type_cyt_pos_filter = ifelse(calc_cyt_pos_gates, "cyt", "base"),
    .data = .data,
    gate_type_single_pos_filter = "base",
    gate_type_single_pos_calc = "base",
    combn = FALSE,
    path_project = path_project,
    .debug = .debug,
    ind_batch_list = ind_batch_list,
    pop_gate = pop_gate
  )
}

.gate_marker_gate_adj_gates_single_stats_tbl_get_check <- function(gate_name_vec) { # nolint
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

.gate_marker_gate_adj_gates_single_out_get <- function(gate_tbl,
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
    .gate_marker_gate_adj_gates_single_out_get_gn(
      .debug = .debug,
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

.gate_marker_gate_adj_gates_single_out_get_gn <- function(.debug,
                                                           gn,
                                                           gate_tbl,
                                                           gate_stats_tbl,
                                                           gate_tbl_ctrl_clust,
                                                           gate_tbl_ctrl_ctrl,
                                                           gate_tbl_single_gn,
                                                           params) {
  .debug_msg(.debug, "gate_name_vec", gn) # nolint
  gate_tbl_gn <- gate_tbl |>
    dplyr::filter(gate_name == gn) # nolint

  adj_ind <- stringr::str_detect(gn, "_adj")

  clust_ind <- stringr::str_detect(gn, "_clust")
  if (!clust_ind && !adj_ind) {
    return(gate_tbl_gn |>
      dplyr::filter(chnl == params$chnl_cut) |> # nolint
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
    gate_tbl_gn_2 <- gate_tbl_gn_2 |>
      dplyr::left_join(
        gate_tbl_single_gn,
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
      .data = .data,
      gate_tbl = gate_tbl_gn,
      gate_stats_tbl = gate_stats_tbl_gn,
      gate_tbl_ctrl = gate_tbl_ctrl_clust_gn,
      chnl = params$chnl_cut,
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
            chnl == params$chnl_cut # nolint
          ) |>
          dplyr::select(
            gate_name, gate_type, gate_combn, # nolint
            batch, ind, gate, gate_cyt # nolint
          ),
        by = c("ind")
      ) |>
      dplyr::mutate(
        chnl = params$chnl_cut,
        marker = params$chnl_lab[params$chnl_cut]
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
