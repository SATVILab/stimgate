#' @param gs GatingSet.
#' @param gate_stats_tbl dataframe.
#' Contains gates and bs_freq info. Must only contain
#' info using "base" threshold to use.
#' @param gate_tbl_ctrl dataframe. Control gates.
#' @param chnl character. Name of channel to get thresholds for.
#' @param bw numeric. Bandwidth for density estimates.
#' @param fcs_vec character. Specifies FCS file names.
#' @param control list. Named list of the following control parameters:
#' - min_threshold_frac -
#' Proportion of minimum threshold to gate on. Default is 0.8.
#' - min_threshold_quant -
#' Quantile of all thresholds to cluster with. Default is 0.1.
#' - bw - numeric. bandwidth for density estimates
#' - n_grp - integer vector of length two. Specifies the numbers of groups
#' to cluster densities into. Default is seq(6, 24, by = 2).
.get_cp_cluster <- function(gs,
                            gate_tbl,
                            gate_stats_tbl,
                            gate_tbl_ctrl,
                            control = list(),
                            chnl,
                            bw,
                            params,
                            filter_other_cyt_pos,
                            debug = FALSE) {
  .debug(debug, "Adjusting thresholds within clusters") # nolint
  # ==================================
  # PREPARATION
  # ==================================

  # control
  control <- .get_cp_cluster_control_update( # nolint
    control
  )

  # statistics
  gate_stats_tbl <- .get_cp_cluster_gate_stats_tbl_update( # nolint
    gate_stats_tbl,
    debug = debug
  )

  # cp
  cp_min <- .get_cp_cluster_cp_get_min( # nolint
    gate_tbl, gate_tbl_ctrl
  )
  max_cp <- .get_cp_cluster_cp_get_max( # nolint
    gate_tbl, gate_tbl_ctrl
  )

  prop_bs_by_cp_tbl_obj <- .get_cp_cluster_prop_bs_by_cp_tbl_obj( # nolint
    gs = gs,
    gate_tbl = gate_tbl,
    ind_batch_list = params$ind_batch_list,
    ind_in_batch_lab_vec = params$ind_in_batch_lab_vec,
    ind_in_batch_uns = params$ind_in_batch_uns,
    high = params$high,
    cut = params$cut,
    pop_gate = params$pop_gate,
    data_name = params$data_name,
    calc_cyt_pos_gates = params$calc_cyt_pos_gates,
    cp_min = cp_min,
    max_cp = max_cp,
    gate_stats_tbl = gate_stats_tbl,
    filter_other_cyt_pos = filter_other_cyt_pos
  )


  prop_bs_by_cp_tbl <- prop_bs_by_cp_tbl_obj[["prop_bs_by_cp_tbl"]]
  expr_max <- prop_bs_by_cp_tbl_obj[["expr_max"]]
  expr_min <- prop_bs_by_cp_tbl_obj[["expr_min"]]

  dens_tbl <- .get_cp_cluster_dens_tbl_get( # nolint
    ind_batch_list = params$ind_batch_list,
    gs = gs,
    ind_in_batch_lab_vec = params$ind_in_batch_lab_vec,
    ind_in_batch_uns = params$ind_in_batch_uns,
    high = params$high,
    data_name = params$data_name,
    filter_other_cyt_pos = filter_other_cyt_pos,
    calc_cyt_pos_gates = params$calc_cyt_pos_gates,
    cut = params$cut,
    expr_min = expr_min,
    expr_max = expr_max,
    pop_gate = params$pop_gate,
    gate_tbl = gate_tbl,
    control = control,
    bw = bw,
    debug = debug
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
    dplyr::group_by(batch_sh, stim, ind) |> # nolint
    dplyr::slice(1) |>
    dplyr::ungroup()
  dens_tbl_join <- dens_tbl_join[, c("ind", "grp")]
  prop_bs_by_cp_tbl <- prop_bs_by_cp_tbl |>
    dplyr::left_join(dens_tbl_join, by = "ind")

  # =====================================
  # CALCULATE THRESHOLDS
  # =====================================

  if (FALSE) {
    .get_cp_cluster_plot_check_1lse(prop_bs_by_cp_tbl)
  }

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
    debug = debug
  )

  gate_summ_stat_tbl <- .get_cp_cluster_gate_summ_stat_tbl_get( # nolint
    gate_tbl = gate_tbl,
    cut = params$cut,
    grp_ind_lab_vec = grp_ind_lab_vec,
    debug = debug
  )

  # calculate thresholds
  # ---------------------------

  # calculate cp_join for each ind

  cp_tbl <- .get_cp_cluster_cp_join_get( # nolint
    prop_bs_by_cp_tbl = prop_bs_by_cp_tbl,
    cp_grp_lab_vec = cp_grp_lab_vec,
    debug = debug
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
    debug = debug
  )

  cp_tbl <- .get_cp_cluster_cp_tbl_add_orig_quant_min( # nolint
    cp_tbl = cp_tbl,
    debug = debug
  )

  # calculate for each individual cp_lse (less than 0.01 standard errors)
  cp_tbl <- .get_cp_cluster_cp_join_lse_get( # nolint
    cp_tbl = cp_tbl,
    debug = debug
  )

  # add tail-gate-based thresholds

  cp_tbl <- .get_cp_cluster_cp_join_tg_get( # nolint
    cp_tbl = cp_tbl,
    debug = debug
  )

  # calculate cp where you don't go all the way to the new cp,
  # but only halfway from original cp
  # (if original cp higher)
  cp_tbl <- cp_tbl |>
    .get_cp_cluster_cp_lse_orig_mean(debug = debug) # nolint

  cp_tbl <- cp_tbl |>
    .get_cp_cluster_cp_join_tg_orig_mean(debug = debug) # nolint

  # calculate cp that is the minimum of lse_orig_mean and tg_orig
  cp_tbl <- cp_tbl |>
    .get_cp_cluster_cp_join_lse_orig_mean_tg(debug = debug) # nolint

  # filter at cp just above cp_join_lse_orig_mean_tg, in order
  # to get the prop_bs_cp_diff closest to it
  cp_tbl <- cp_tbl |>
    .get_cp_cluster_cp_filter_above_cp_join_lse_orig_mean_tg(debug = debug) # nolint

  # if no gate found above, then set it to base threshold OR
  # impute based on group. Not going to work when original threshold is NA.

  cp_tbl <- .get_cp_cluster_cp_impute_missing_batch( # nolint
    cp_tbl = cp_tbl,
    cut = params$cut,
    gate_tbl = gate_tbl,
    dens_tbl = dens_tbl,
    debug = debug
  )

  cp_tbl <- .get_cp_cluster_impute_missing_ind( # nolint
    cp_tbl = cp_tbl,
    cut = params$cut,
    gate_tbl = gate_tbl,
    dens_tbl = dens_tbl,
    debug = debug
  )

  cp_tbl <- .get_cp_cluster_impute_missing_final( # nolint
    cp_tbl = cp_tbl,
    cut = params$cut,
    gate_tbl = gate_tbl,
    dens_tbl = dens_tbl,
    debug = debug
  )

  cp_tbl <- .get_cp_cluster_impute_missing_final_batch( # nolint
    cp_tbl = cp_tbl,
    cut = params$cut,
    gate_tbl = gate_tbl,
    dens_tbl = dens_tbl,
    debug = debug
  )

  # ================================
  # PLOT THRESHOLDS
  # ================================

  if (FALSE) {
    .get_cp_cluster_plot_thresholds(
      cp_tbl = cp_tbl,
      params = params,
      marker_list = marker_list,
      chnl = chnl,
      gs = gs,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec
    )
  }

  # =========================
  # OUTPUT
  # =========================

  cp_tbl
}
