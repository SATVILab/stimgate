# Get cutpoints using clustering approach
# Clusters thresholds from similar distributions to identify optimal cutpoints
.get_cp_cluster <- function(.data,
                            gate_tbl,
                            gate_stats_tbl,
                            gate_tbl_ctrl,
                            control = list(),
                            chnl,
                            bw,
                            params,
                            filter_other_cyt_pos,
                            .debug = FALSE) {
  .debug_msg(.debug, "Adjusting thresholds within clusters") # nolint
  # ==================================
  # PREPARATION
  # ==================================

  # control
  control <- .get_cp_cluster_control_update(control) # nolint

  # statistics
  gate_stats_tbl <- .get_cp_cluster_gate_stats_tbl_update( # nolint
    gate_stats_tbl, .debug
  )

  # cp
  cp_min <- .get_cp_cluster_cp_get_min(gate_tbl, gate_tbl_ctrl) # nolint
  max_cp <- .get_cp_cluster_cp_get_max(gate_tbl, gate_tbl_ctrl) # nolint

  prop_bs_by_cp_tbl_obj <- .get_cp_cluster_prop_bs_by_cp_tbl_obj( # nolint
    .data = .data,
    gate_tbl = gate_tbl,
    ind_batch_list = params$ind_batch_list,
    chnl_cut = params$chnl_cut,
    pop_gate = params$pop_gate,
    calc_cyt_pos_gates = params$calc_cyt_pos_gates,
    cp_min = cp_min,
    max_cp = max_cp,
    gate_stats_tbl = gate_stats_tbl,
    filter_other_cyt_pos = filter_other_cyt_pos,
    .debug = .debug
  )


  prop_bs_by_cp_tbl <- prop_bs_by_cp_tbl_obj[["prop_bs_by_cp_tbl"]]
  expr_max <- prop_bs_by_cp_tbl_obj[["expr_max"]]
  expr_min <- prop_bs_by_cp_tbl_obj[["expr_min"]]

  dens_tbl <- .get_cp_cluster_dens_tbl_get( # nolint
    ind_batch_list = params$ind_batch_list,
    .data = .data,
    filter_other_cyt_pos = filter_other_cyt_pos,
    calc_cyt_pos_gates = params$calc_cyt_pos_gates,
    chnl_cut = params$chnl_cut,
    expr_min = expr_min,
    expr_max = expr_max,
    pop_gate = params$pop_gate,
    gate_tbl = gate_tbl,
    control = control,
    bw = bw,
    .debug = .debug
  )

  n_clus <- .get_cp_cluster_n_clus( # nolint
    dens_tbl
  )
  clus_vec <- .get_cp_cluster_clus( # nolint
    dens_tbl, n_clus
  )

  dens_tbl <- dens_tbl |>
    dplyr::mutate(grp = as.character(clus_vec))

  dens_tbl_grp_lab <- dens_tbl |>
    dplyr::group_by(grp, ind) |> # nolint
    dplyr::slice(1) |>
    dplyr::ungroup()

  grp_ind_lab_vec <- dens_tbl_grp_lab$grp |>
    stats::setNames(dens_tbl_grp_lab$ind)

  if (FALSE) {
    .get_cp_cluster_plot_check_1(dens_tbl)
  }

  # bind onto the summary statistics the group
  # labels
  dens_tbl_join <- dens_tbl |>
    dplyr::group_by(batch, ind) |> # nolint
    dplyr::slice(1) |>
    dplyr::ungroup()
  dens_tbl_join <- dens_tbl_join[, c("ind", "grp")]
  prop_bs_by_cp_tbl <- prop_bs_by_cp_tbl |>
    dplyr::left_join(dens_tbl_join, by = "ind")

  # =====================================
  # CALCULATE THRESHOLDS
  # =====================================

  # length table,
  # making a column group that contains
  # the total number of groups
  # and are the group allocations
  # but this makes no sense, given the "grp" column
  # above?
  # fit smoothed models to proportion with se <= 1
  # --------------------------

  cp_grp_lab_vec <- .get_cp_cluster_cp_grp_lab_vec_get( # nolint
    prop_bs_by_cp_tbl = prop_bs_by_cp_tbl,
    expr_max = expr_max,
    expr_min = expr_min,
    .debug = .debug
  )

  gate_summ_stat_tbl <- .get_cp_cluster_gate_summ_stat_tbl_get( # nolint
    gate_tbl = gate_tbl,
    chnl_cut = params$chnl_cut,
    grp_ind_lab_vec = grp_ind_lab_vec,
    .debug = .debug
  )

  # calculate thresholds
  # ---------------------------

  # calculate cp_join for each ind

  cp_tbl <- .get_cp_cluster_cp_join_get( # nolint
    prop_bs_by_cp_tbl = prop_bs_by_cp_tbl,
    cp_grp_lab_vec = cp_grp_lab_vec,
    .debug = .debug
  )

  gate_tbl_chnl <- .get_cp_cluster_gate_tbl_chnl_get( # nolint
    gate_tbl, chnl
  )

  # add other tables with useful information
  # for creating new gates
  cp_tbl <- .get_cp_cluster_cp_tbl_add_info( # nolint
    cp_tbl = cp_tbl,
    gate_stats_tbl = gate_stats_tbl,
    gate_summ_stat_tbl = gate_summ_stat_tbl,
    gate_tbl_ctrl = gate_tbl_ctrl,
    gate_tbl_chnl = gate_tbl_chnl,
    .debug = .debug
  )

  cp_tbl <- .get_cp_cluster_cp_tbl_add_orig_quant_min( # nolint
    cp_tbl = cp_tbl,
    .debug = .debug
  )

  # calculate for each individual cp_lse (less than 0.01 standard errors)
  cp_tbl <- .get_cp_cluster_cp_join_lse_get( # nolint
    cp_tbl = cp_tbl,
    .debug = .debug
  )

  # add tail-gate-based thresholds

  cp_tbl <- .get_cp_cluster_cp_join_tg_get( # nolint
    cp_tbl = cp_tbl,
    .debug = .debug
  )

  # calculate cp where you don't go all the way to the new cp,
  # but only halfway from original cp
  # (if original cp higher)
  cp_tbl <- cp_tbl |>
    .get_cp_cluster_cp_lse_orig_mean(.debug = .debug) # nolint

  cp_tbl <- cp_tbl |>
    .get_cp_cluster_cp_join_tg_orig_mean(.debug = .debug) # nolint

  # calculate cp that is the minimum of lse_orig_mean and tg_orig
  cp_tbl <- cp_tbl |>
    .get_cp_cluster_cp_join_lse_orig_mean_tg(.debug = .debug) # nolint

  # filter at cp just above cp_join_lse_orig_mean_tg, in order
  # to get the prop_bs_cp_diff closest to it
  cp_tbl <- cp_tbl |>
    .get_cp_cluster_cp_filter_above_cp_join_lse_orig_mean_tg(.debug = .debug) # nolint

  # if no gate found above, then set it to base threshold OR
  # impute based on group. Not going to work when original threshold is NA.

  cp_tbl <- .get_cp_cluster_cp_impute_missing_batch( # nolint
    cp_tbl = cp_tbl,
    chnl_cut = params$chnl_cut,
    gate_tbl = gate_tbl,
    dens_tbl = dens_tbl,
    .debug = .debug
  )

  cp_tbl <- .get_cp_cluster_impute_missing_ind( # nolint
    cp_tbl = cp_tbl,
    chnl_cut = params$chnl_cut,
    gate_tbl = gate_tbl,
    dens_tbl = dens_tbl,
    .debug = .debug
  )

  cp_tbl <- .get_cp_cluster_impute_missing_final( # nolint
    cp_tbl = cp_tbl,
    chnl_cut = params$chnl_cut,
    gate_tbl = gate_tbl,
    dens_tbl = dens_tbl,
    .debug = .debug
  )

  cp_tbl <- .get_cp_cluster_impute_missing_final_batch( # nolint
    cp_tbl = cp_tbl,
    chnl_cut = params$chnl_cut,
    gate_tbl = gate_tbl,
    dens_tbl = dens_tbl,
    .debug = .debug
  )

  # =========================
  # OUTPUT
  # =========================

  cp_tbl
}
