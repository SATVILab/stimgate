#' @keywords internal
.gate_chnl_delete_old_gates <- function() {
  dir_save <- file.path(tempdir(), "stimgate")
  if (!dir.exists(dir_save)) {
    return(invisible(FALSE))
  }
  unlink(dir_save, recursive = TRUE)
  invisible(TRUE)
}

# Get gates for each sample within each batch
#' @keywords internal
.gate_chnl_pre_adj_gates_gate <- function(
  ind_batch_list,
  .data,
  chnl_settings,
  stage,
  path_project
) {
  message("getting pre-adjustment gates")
  purrr::map_df(seq_along(ind_batch_list), function(i) {
    .debug("ind_batch_list", i) # nolint

    # message progress
    if (i %% 50 == 0 || i == length(ind_batch_list)) {
      txt <- paste0("batch ", i, " of ", length(ind_batch_list))
      message(txt)
    }
    .gate_batch(
      # nolint
      .data = .data,
      ind_batch = ind_batch_list[[i]],
      batch = names(ind_batch_list)[i],
      chnl_settings = chnl_settings,
      stage = stage,
      path_project = path_project
    ) |>
      dplyr::select(
        gate_name,
        gate_type,
        gate_combn, # nolint
        batch,
        ind,
        gate,
        everything() # nolint
      )
  })
}

#' @keywords internal
.gate_chnl_get_adj_gates <- function(
  gate_tbl,
  gate_tbl_params,
  chnl_settings,
  .data,
  stage,
  path_project,
  ind_batch_list,
  calc_cyt_pos_gates
) {
  if (stage == "init") {
    .gate_chnl_get_adj_gates_all(
      # nolint
      gate_tbl = gate_tbl,
      .data = .data,
      path_project = path_project,
      stage = stage,
      ind_batch_list = ind_batch_list,
      chnl_settings = chnl_settings,
      calc_cyt_pos_gates = calc_cyt_pos_gates
    )
  } else if (stage == "single") {
    .gate_chnl_gate_adj_gates_single(
      gate_tbl = gate_tbl,
      gate_tbl_params = gate_tbl_params,
      .data = .data,
      calc_cyt_pos_gates = TRUE,
      path_project = path_project,
      ind_batch_list = ind_batch_list,
      stage = stage,
      chnl_settings = chnl_settings
    )
  } else {
    stop("stage not recognized")
  }
}

#' @keywords internal
.gate_chnl_get_adj_gates_all <- function(
  gate_tbl,
  .data,
  path_project,
  stage,
  ind_batch_list,
  chnl_settings,
  calc_cyt_pos_gates
) {
  if (!is.null(chnl_settings$tol_clust)) {
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

  if (!is.null(chnl_settings$tol_clust)) {
    path_dir_stats <- .get_stats(
      # nolint
      gate_tbl = gate_tbl |> dplyr::mutate(chnl = chnl_settings$chnl_cut),
      gate_name = NULL,
      chnl = chnl_settings$chnl_cut,
      .data = .data,
      filter_other_cyt_pos = FALSE,
      combn = FALSE,
      gate_type_single_pos_calc = "base",
      path_project = path_project,
      ind_batch_list = ind_batch_list,
      pop_gate = chnl_settings$pop_gate,
      tol_clust = chnl_settings$tol_clust
    )
    gate_stats_tbl <- path_dir_stats |>
      .read_gate_stats() # nolint

    gate_tbl_cluster <- purrr::map_df(
      unique(gate_tbl$gate_name),
      function(gn) {
        gate_tbl_cluster <- .get_cp_cluster(
          # nolint
          .data = .data,
          gate_tbl = gate_tbl |>
            dplyr::filter(gate_name == gn), # nolint
          gate_stats_tbl = gate_stats_tbl |>
            dplyr::filter(gate_name == gn),
          gate_tbl_ctrl = gate_tbl_tg_gate,
          chnl_settings = chnl_settings,
          filter_other_cyt_pos = FALSE,
          stage = stage,
          path_project = path_project,
          calc_cyt_pos_gates = calc_cyt_pos_gates,
          ind_batch_list = ind_batch_list
        )

        gate_tbl_cluster |>
          dplyr::select(ind, cp_join_tg_orig) |> # nolint
          dplyr::rename(gate = cp_join_tg_orig) |>
          dplyr::left_join(
            gate_tbl |>
              dplyr::filter(gate_name == gn) |> # nolint
              dplyr::select(
                gate_name,
                gate_type,
                gate_combn, # nolint
                batch,
                ind # nolint
              ),
            by = c("ind")
          ) |>
          dplyr::select(
            gate_name,
            gate_type,
            gate_combn,
            batch,
            ind,
            gate # nolint
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

  list(
    gate_tbl = gate_tbl
  )
}

#' @keywords internal
.gate_chnl_gate_adj_gates_single <- function(
  gate_tbl,
  gate_tbl_params,
  .data,
  calc_cyt_pos_gates,
  path_project,
  ind_batch_list,
  stage,
  chnl_settings
) {
  # get gates
  gate_tbl_single <- gate_tbl

  # merge
  gate_tbl <- .gate_chnl_gate_adj_gates_single_merge(
    # nolint
    gate_tbl_single = gate_tbl_single,
    gate_tbl_params = gate_tbl_params,
    chnl_cut = chnl_settings$chnl_cut
  )

  # get stats table (if needed)
  gate_stats_tbl <- .gate_chnl_gate_adj_gates_single_stats_tbl_get(
    gate_tbl = gate_tbl,
    .data = .data,
    chnl_settings = chnl_settings,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    path_project = path_project,
    ind_batch_list = ind_batch_list
  )

  gate_tbl_out <- .gate_chnl_gate_adj_gates_single_out_get(
    gate_tbl = gate_tbl,
    gate_stats_tbl = gate_stats_tbl,
    gate_tbl_single = gate_tbl_single,
    gate_quant = gate_quant,
    .data = .data,
    stage = stage,
    path_project = path_project,
    chnl_settings = chnl_settings
  )

  list(
    gate_tbl = gate_tbl_out
  )
}

#' @keywords internal
.gate_chnl_gate_adj_gates_single_merge <- function(
  gate_tbl_single,
  gate_tbl_params,
  chnl_cut
) {
  gate_tbl_params |>
    dplyr::left_join(
      gate_tbl_single |>
        dplyr::mutate(chnl = chnl_cut) |>
        dplyr::rename(gate_single = gate) |> # nolint
        dplyr::filter(gate_use == "gate") |> # nolint
        dplyr::select(
          chnl,
          gate_name,
          ind,
          gate_single # nolint
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
      chnl,
      marker,
      gate_name, # nolint
      gate_type,
      gate_combn,
      everything() # nolint
    ) |>
    dplyr::mutate(
      gate_single = ifelse(is.na(gate_single), gate, gate_single) # nolint
    )
}

#' @keywords internal
.gate_chnl_gate_adj_gates_single_stats_tbl_get <- function(
  gate_tbl,
  chnl_settings,
  .data,
  calc_cyt_pos_gates,
  path_project,
  ind_batch_list
) {
  gate_name_vec <- unique(gate_tbl$gate_name)
  if (!.gate_chnl_gate_adj_gates_single_stats_tbl_get_check(gate_name_vec)) {
    return(NULL)
  }

  gate_name_vec_clust <- gate_name_vec[
    stringr::str_detect(gate_name_vec, "_clust")
  ]
  gate_name_vec_adj <- gate_name_vec[
    stringr::str_detect(gate_name_vec, "_adj")
  ]

  .get_stats(
    # nolint
    gate_tbl = gate_tbl |>
      dplyr::filter(
        gate_name %in% c(gate_name_vec_clust, gate_name_vec_adj) # nolint
      ),
    gate_name = NULL,
    chnl = chnl_settings$chnl_cut,
    filter_other_cyt_pos = TRUE,
    gate_type_cyt_pos_filter = if (calc_cyt_pos_gates) "cyt" else "base",
    .data = .data,
    gate_type_single_pos_filter = "base",
    gate_type_single_pos_calc = "base",
    combn = FALSE,
    path_project = path_project,
    ind_batch_list = ind_batch_list,
    pop_gate = chnl_settings$chnl_cutpop_gate
  )
}

#' @keywords internal
.gate_chnl_gate_adj_gates_single_stats_tbl_get_check <- function(
  gate_name_vec
) {
  # nolint
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

#' @keywords internal
.gate_chnl_gate_adj_gates_single_out_get <- function(
  gate_tbl,
  gate_stats_tbl,
  gate_tbl_single,
  gate_quant,
  .data,
  stage,
  path_project,
  chnl_settings
) {
  # get tail-gate gates
  gate_tbl_ctrl_clust <- gate_tbl_single |>
    dplyr::filter(gate_use == "tg_clust") # nolint

  # get control gates
  gate_tbl_ctrl_ctrl <- gate_tbl_single |> # nolint
    dplyr::filter(gate_use == "ctrl") # nolint

  # get gate names
  gate_name_vec <- unique(gate_tbl$gate_name)
  purrr::map_df(gate_name_vec, function(gn) {
    .gate_chnl_gate_adj_gates_single_out_get_gn(
      gn = gn,
      gate_tbl = gate_tbl,
      gate_stats_tbl = gate_stats_tbl,
      gate_tbl_ctrl_clust = gate_tbl_ctrl_clust,
      gate_tbl_ctrl_ctrl = gate_tbl_ctrl_ctrl,
      gate_tbl_single_gn = gate_tbl_single |>
        dplyr::filter(gate_name == gn), # nolint
      gate_quant = gate_quant,
      .data = .data,
      stage = stage,
      path_project = path_project,
      chnl_settings = chnl_settings
    )
  }) |>
    dplyr::mutate(gate_single = pmax(gate, gate_single)) # nolint
}

#' @keywords internal
.gate_chnl_gate_adj_gates_single_out_get_gn <- function(
  gn,
  gate_tbl,
  gate_stats_tbl,
  gate_tbl_ctrl_clust,
  gate_tbl_ctrl_ctrl,
  gate_tbl_single_gn,
  gate_quant,
  .data,
  stage,
  path_project,
  chnl_settings
) {
  .debug("gate_name_vec", gn) # nolint
  gate_tbl_gn <- gate_tbl |>
    dplyr::filter(gate_name == gn) # nolint

  adj_ind <- stringr::str_detect(gn, "_adj")

  clust_ind <- stringr::str_detect(gn, "_clust")
  if (!clust_ind && !adj_ind) {
    return(
      gate_tbl_gn |>
        dplyr::filter(chnl == chnl_settings$chnl_cut) |> # nolint
        dplyr::select(
          chnl,
          marker,
          gate_name,
          gate_type, # nolint
          gate_combn,
          batch,
          ind,
          gate,
          gate_cyt,
          gate_single # nolint
        )
    )
  }

  # =========================
  # Tail-gate controlled gating
  # =========================

  if (adj_ind) {
    gate_stats_tbl_gn <- gate_stats_tbl |>
      dplyr::filter(gate_name == gn) # nolint
    gate_tbl_ctrl_ctrl <- gate_tbl_ctrl_ctrl |>
      dplyr::filter(gate_name == gn) # nolint

    gate_tbl_gn_2 <- .get_cp_adj_tbl(
      # nolint
      gate_stats_tbl = gate_stats_tbl_gn,
      gate_quant = gate_quant,
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

    gate_tbl_cluster_gn <- .get_cp_cluster(
      # nolint
      .data = .data,
      gate_tbl = gate_tbl_gn,
      gate_stats_tbl = gate_stats_tbl_gn,
      gate_tbl_ctrl = gate_tbl_ctrl_clust_gn,
      filter_other_cyt_pos = TRUE,
      stage = stage,
      path_project = path_project
    )

    gate_tbl_cluster_gn <- gate_tbl_cluster_gn |>
      dplyr::select(ind, cp_join_tg_orig) |> # nolint
      dplyr::rename(gate_single = cp_join_tg_orig) |>
      dplyr::left_join(
        gate_tbl |>
          dplyr::filter(
            gate_name == gn, # nolint
            chnl == chnl_settings$chnl_cut # nolint
          ) |>
          dplyr::select(
            gate_name,
            gate_type,
            gate_combn, # nolint
            batch,
            ind,
            gate,
            gate_cyt # nolint
          ),
        by = c("ind")
      ) |>
      dplyr::mutate(
        chnl = chnl_settings$chnl_cut,
        marker = chnl_settings$marker
      ) |>
      dplyr::select(
        chnl,
        marker,
        gate_name,
        gate_type,
        gate_combn, # nolint
        batch,
        ind,
        gate,
        gate_cyt,
        gate_single # nolint
      ) |>
      dplyr::mutate(
        gate_combn = paste0(gate_combn, "_clust"),
        gate_name = paste0(gate_type, "_", gate_combn)
      )
  }

  gate_tbl_cluster_gn
}

# Placeholder function for adjustment table - this may need proper implementation
#' @keywords internal
.get_cp_adj_tbl <- function(gate_stats_tbl, gate_quant, gate_tbl_ctrl) {
  unique_ind <- unique(gate_stats_tbl$ind)
  tibble::tibble(
    ind = unique_ind,
    gate = rep(0.5, length(unique_ind)) # placeholder gate values
  )
}
