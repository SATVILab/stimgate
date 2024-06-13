.get_gate_stats <- function(params = NULL,
                            gate_tbl = NULL,
                            chnl = NULL,
                            filter_other_cyt_pos = FALSE,
                            gate_name = NULL,
                            combn = TRUE,
                            gate_type_cyt_pos_filter = "base",
                            gate_type_single_pos_filter = "base",
                            gate_type_cyt_pos_calc,
                            gate_type_single_pos_calc,
                            debug = FALSE,
                            pop_gate,
                            chnl_lab = NULL,
                            data,
                            ind_in_batch_lab,
                            ind_in_batch_gate,
                            fcs,
                            data_name,
                            ind_in_batch_uns,
                            save = FALSE,
                            ind_batch_list,
                            save_gate_tbl = FALSE,
                            stim_lab = NULL,
                            sampleid_lab = NULL,
                            path_project) {
  # prep
  # ---------------

  chnl_lab <- .get_gate_stats_chnl_lab_get( # nolint
    chnl_lab = chnl_lab,
    data = data,
    chnl = chnl
  )

  params <- .get_gate_stats_params_get( # nolint
    params = params,
    chnl = chnl,
    pop_gate = pop_gate,
    chnl_lab = chnl_lab,
    ind_in_batch_lab = ind_in_batch_lab,
    ind_in_batch_gate = ind_in_batch_gate,
    data_name = data_name,
    fcs = fcs,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_batch_list = ind_batch_list,
    data = data
  )
  gate_tbl <- .get_gate_stats_gate_tbl_get( # nolint
    gate_tbl = gate_tbl,
    chnl_lab = chnl_lab,
    path_project = path_project,
    params = params,
    gate_name = gate_name
  )
  chnl <- .get_gate_stats_chnl_get( # nolint
    chnl = chnl,
    gate_tbl = gate_tbl
  )
  gate_name <- .get_gate_stats_gate_name_get( # nolint
    gate_name = gate_name,
    gate_tbl = gate_tbl
  )

  if ((!filter_other_cyt_pos) && combn) {
    n_chnl <- length(chnl)
    combn_mat_list <-
      .get_gate_stats_combn_mat_list_get( # nolint
        n_chnl = n_chnl,
        n_pos = 2 # not sure if it should be 2,
        # but it wasn't really set before
      )
    cyt_combn_vec_list <-
      .get_gate_stats_cyt_combn_vec_list_get( # nolint
        combn_mat_list = combn_mat_list,
        chnl = chnl
      )
  } else {
    combn_mat_list <- NULL
    cyt_combn_vec_list <- NULL
  }

  .get_gate_stats_gate_tbl_save( # nolint
    gate_tbl = gate_tbl,
    path_project = path_project,
    params = params,
    chnl = chnl,
    save = save_gate_tbl,
    sampleid_lab = sampleid_lab,
    stim_lab = stim_lab
  )

  stat_tbl <- .get_gate_stats_overall( # nolint
    ind_batch_list = ind_batch_list,
    gate_tbl = gate_tbl,
    chnl = chnl,
    combn = combn,
    cyt_combn_vec_list = cyt_combn_vec_list,
    gate_name = gate_name,
    gate_type_cyt_pos_calc = gate_type_cyt_pos_calc,
    gate_type_single_pos_calc = gate_type_single_pos_calc,
    gate_type_cyt_pos_filter = gate_type_cyt_pos_filter,
    gate_type_single_pos_filter = gate_type_single_pos_filter,
    pop_gate = pop_gate,
    data = data,
    chnl_lab = chnl_lab,
    cut = params$cut,
    ind_in_batch_lab = ind_in_batch_lab,
    data_name = data_name,
    ind_in_batch_uns = ind_in_batch_uns,
    sampleid_lab = sampleid_lab,
    stim_lab = stim_lab,
    debug = debug,
    filter_other_cyt_pos = filter_other_cyt_pos,
    combn_mat_list = combn_mat_list
  )

  # save it
  .get_gate_stats_save( # nolint
    stat_tbl = stat_tbl,
    path_project = path_project,
    params = params,
    save = save,
    chnl = chnl
  )
}

.read_gate_stats <- function(path_dir_stats) {
  path_dir_stats |>
    file.path("gate_stats.rds") |>
    readRDS()
}
