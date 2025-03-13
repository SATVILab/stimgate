.get_stats <- function(params = NULL,
                       gate_tbl = NULL,
                       chnl = NULL,
                       filter_other_cyt_pos = FALSE,
                       combn = TRUE,
                       gate_type_cyt_pos_filter = "base",
                       gate_type_single_pos_filter = "base",
                       gate_type_cyt_pos_calc,
                       gate_type_single_pos_calc,
                       debug = FALSE,
                       pop_gate,
                       chnl_lab = NULL,
                       .data,
                       save = FALSE,
                       ind_batch_list,
                       save_gate_tbl = FALSE,
                       gate_name = NULL,
                       tol_clust = NULL,
                       path_project) {
  # prep
  # ---------------

  chnl_lab <- .get_stats_chnl_lab_get( # nolint
    chnl_lab = chnl_lab,
    .data = .data,
    chnl = chnl
  )

  params <- .get_stats_params_get( # nolint
    params = params,
    chnl = chnl,
    pop_gate = pop_gate,
    chnl_lab = chnl_lab,
    ind_batch_list = ind_batch_list,
    .data = .data
  )
  gate_tbl <- .get_stats_gate_tbl_get( # nolint
    gate_tbl = gate_tbl,
    chnl_lab = chnl_lab,
    path_project = path_project,
    params = params,
    gate_name = gate_name,
    tol_clust = tol_clust
  )
  chnl <- .get_stats_chnl_get( # nolint
    chnl = chnl,
    gate_tbl = gate_tbl
  )
  gate_name <- .get_stats_gate_name_get( # nolint
    gate_name = gate_name,
    gate_tbl = gate_tbl
  )

  if ((!filter_other_cyt_pos) && combn) {
    n_chnl <- length(chnl)
    combn_mat_list <-
      .get_stats_combn_mat_list_get( # nolint
        n_chnl = n_chnl,
        n_pos = 2 # not sure if it should be 2,
        # but it wasn't really set before
      )
    cyt_combn_vec_list <-
      .get_stats_cyt_combn_vec_list_get( # nolint
        combn_mat_list = combn_mat_list,
        chnl = chnl
      )
  } else {
    combn_mat_list <- NULL
    cyt_combn_vec_list <- NULL
  }

  .get_stats_gate_tbl_save( # nolint
    gate_tbl = gate_tbl,
    path_project = path_project,
    params = params,
    chnl = chnl,
    save = save_gate_tbl
  )

  stat_tbl <- .get_stats_overall( # nolint
    ind_batch_list = ind_batch_list,
    gate_tbl = gate_tbl,
    chnl = chnl,
    combn = combn,
    cyt_combn_vec_list = cyt_combn_vec_list,
    gate_type_cyt_pos_calc = gate_type_cyt_pos_calc,
    gate_type_single_pos_calc = gate_type_single_pos_calc,
    gate_type_cyt_pos_filter = gate_type_cyt_pos_filter,
    gate_type_single_pos_filter = gate_type_single_pos_filter,
    pop_gate = pop_gate,
    .data = .data,
    chnl_lab = chnl_lab,
    chnl_cut = params$chnl_cut,
    debug = debug,
    filter_other_cyt_pos = filter_other_cyt_pos,
    combn_mat_list = combn_mat_list,
    gate_name = gate_name
  )

  # save it
  .stats_save( # nolint
    stat_tbl = stat_tbl,
    path_project = path_project,
    params = params,
    save = save,
    chnl = chnl
  )
}

.read_gate_stats <- function(stats_save_output) {
  if (inherits(stats_save_output, "data.frame")) {
    return(stats_save_output)
  }
  if (!inherits(stats_save_output, "character")) {
    stop(
      "stats_save_output must be a character string if not a data.frame."
    )
  }
  path_stats <- file.path(stats_save_output, "gate_stats.rds")
  gate_stats_tbl <- readRDS(path_stats)
  if ("ind" %in% colnames(gate_stats_tbl)) {
    gate_stats_tbl[, "ind"] <- as.character(gate_stats_tbl[["ind"]])
  }
  if ("batch" %in% colnames(gate_stats_tbl)) {
    gate_stats_tbl[, "batch"] <- as.character(gate_stats_tbl[["batch"]])
  }
  gate_stats_tbl
}

#' @title Get gating statistics
#' @params path_project character. Path to the project directory.
#' @return A data frame with gating statistics.
#' @export
get_stats <- function(path_project) {
  path_stats_partial <- file.path(path_project, "gate_stats")
  if (file.exists(paste0(path_stats_partial, ".rds"))) {
    readRDS(paste0(path_stats_partial, ".rds"))
  } else if (file.exists(paste0(path_stats_partial, ".csv"))) {
    read.csv(paste0(path_stats_partial, ".csv"))
  } else {
    stop(
      "No stats file found"
    )
  }
}