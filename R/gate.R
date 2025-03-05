#' @title Get gates and stats and plots based thereon
#'
#' @description
#'
#' @param .data GatingSet. GatingSet from which to draw the data.
#' @param pop_gate character vector. Populations for which separate gates are
#'   to be calculated.
#' @param path_project character. Path to project directory. Results are saved
#'   here.
#' @param batch_size numeric. Number of samples per batch. Batches are
#'   calculated assuming that numbers 1 to \code{batch_size} are batch one,
#'   (\code{batch_siz} + 1) to 2 * \code{batch_size} are batch two, and so on.
#' @param ind_in_batch_uns integer. Index within a given batch corresponding to
#'   unstim. For example, if \code{ind_in_batch_uns} is four, then the fourth
#'   sample within each batch is regarded as unstim.
#' @param ind_in_batch_gate integer vector. Indices within a given batch for
#'   which gates need to be calculated.
#'
#' # gating parameters
#' @param marker list. List where each element specifies the parameters to be
#'   used to gate a given marker (including the name of the marker itself).
#'   List elements are as follows:
#'   \describe{
#'     \item{cut}{character. Name of channel to get gate for.}
#'     \item{tol}{numeric. Tolerance value for the
#'       \code{cytoUtils:::.cytokine_cutpoint} method. Default is 0.5e-8 for
#'       CyTOF and flow.}
#'     \item{gate_combn}{named list. Named list where names are one of 'grp',
#'       'mean', 'median', 'trim20', 'min' or 'max', and elements are character
#'       vectors of 'scp', 'dcp', 'tg', 'midp', 'uns#' and 'uns#r' (where # are
#'       FDRs expressed as percentages). Each element therefore specifies the
#'       method of combining gates from individual samples within a group for a
#'       subset of the automated gating methods. \code{'grp'} means to gate
#'       jointly, whereas all of the others do what they sound like. If not
#'       specified (i.e. \code{NULL}), then all gates are performed
#'       individually on each sample. Default is \code{NULL}.
#'   }

#' @importFrom flowCore exprs<- parameters<-
#' @export
stimgate_gate <- function(.data,
                          path_project,
                          pop_gate,
                          marker,
                          batch_list,
                          bias_uns = NULL,
                          cp_min = NULL,
                          bw_min = NULL,
                          min_cell = 10,
                          max_pos_prob_x,
                          gate_quant = c(0.25, 0.75),
                          tol_ctrl = NULL,
                          tol_gate = NULL,
                          calc_cyt_pos_gates = TRUE,
                          calc_single_pos_gates = TRUE,
                          tol_gate_single = NULL,
                          debug = FALSE) {
  force(.data)

  # get unspecified levels in marker elements
  marker <- .complete_marker_list( # nolint
    marker = marker,
    .data = .data,
    pop_gate = pop_gate,
    cut = cut,
    debug = debug,
    bias_uns = bias_uns,
    bw_min = bw_min,
    cp_min = cp_min,
    ind_batch_list = batch_list
  )

  # inital gates
  .gate_init(
    pop_gate = pop_gate,
    marker = marker,
    .data = .data,
    ind_batch_list = batch_list,
    path_project = path_project,
    noise_sd = NULL,
    max_pos_prob_x = max_pos_prob_x,
    gate_quant = gate_quant,
    tol_ctrl = tol_ctrl,
    tol_gate = tol_gate,
    tol_gate_single = tol_gate_single,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    debug = debug
  )

  # cytokine-positive gates
  gate_tbl <- .gate_cyt_pos( # nolint
    marker_list = marker,
    pop_gate = pop_gate,
    .data = .data,
    calc_cyt_pos = calc_cyt_pos_gates,
    debug = debug,
    path_project = path_project
  )

  # single-positive gates
  .gate_single(
    pop_gate = pop_gate,
    marker = marker,
    .data = .data,
    ind_batch_list = batch_list,
    path_project = path_project,
    noise_sd = NULL,
    max_pos_prob_x = max_pos_prob_x,
    gate_quant = gate_quant,
    tol_ctrl = tol_ctrl,
    tol_gate = tol_gate,
    tol_gate_single = tol_gate_single,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    debug = debug,
    gate_tbl = gate_tbl
  )

  print("")
  print("")
  print("")
  print("getting cyt combn frequencies")

  path_dir_stats <- .gate_stats(
    .data = .data,
    params = NULL,
    gate_tbl = NULL,
    filter_other_cyt_pos = FALSE,
    combn = TRUE,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    debug = debug,
    save = TRUE,
    pop_gate = pop_gate,
    marker = marker,
    ind_batch_list = batch_list,
    path_project = path_project,
    save_gate_tbl = TRUE
  )

  path_dir_stats
}

.gate_init <- function(pop_gate,
                       marker,
                       .data,
                       ind_batch_list,
                       path_project,
                       noise_sd,
                       max_pos_prob_x,
                       gate_quant,
                       tol_ctrl,
                       tol_gate,
                       tol_gate_single,
                       calc_cyt_pos_gates,
                       debug) {
  # loop over populations
  print("----")
  print("getting base gates")
  print("----")
  print("")
  # loop over markers
  purrr::walk(marker, function(marker_curr) {
    print(paste0("chnl: ", marker_curr$cut))
    # get gates for each sample within each batch

    gate_obj <- .gate_marker( # nolint
      .data = .data,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate,
      cut = marker_curr$cut,
      gate_combn = marker_curr$gate_combn,
      tol = marker_curr$tol,
      noise_sd = NULL,
      bias_uns = marker_curr$bias_uns,
      bw_min = marker_curr$bw_min,
      cp_min = marker_curr$cp_min,
      min_cell = marker_curr$min_cell,
      max_pos_prob_x = max_pos_prob_x,
      gate_quant = gate_quant,
      tol_ctrl = tol_ctrl,
      tol_gate = tol_gate,
      tol_gate_single = tol_gate_single,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      path_project = path_project,
      debug = debug
    )

    saveRDS(
      gate_obj$gate_tbl,
      file = file.path(
        path_project, marker_curr$cut, "gate_tbl_init.rds"
      )
    )
  })
}

.gate_single <- function(pop_gate,
                         marker,
                         .data,
                         ind_batch_list,
                         path_project,
                         noise_sd,
                         max_pos_prob_x,
                         gate_quant,
                         tol_ctrl,
                         tol_gate,
                         tol_gate_single,
                         calc_cyt_pos_gates,
                         calc_single_pos_gates,
                         debug,
                         gate_tbl) {
  # loop over populations
  print("")
  print("")
  print("----")
  print("getting single+ gates")
  print("----")
  print("")
  if (!calc_single_pos_gates) {
    .debug(debug, "Not gating single-pos gates") # nolint
    return(invisible(TRUE))
  } else {
    .debug(debug, "Gating single-pos gates") # nolint
  }
  # loop over markers
  purrr::walk(marker, function(marker_curr) {
    print(paste0("chnl: ", marker_curr$cut))
    # get gates for each sample within each batch

    gate_obj <- .gate_marker( # nolint
      .data = .data,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate_curr,
      cut = marker_curr$cut,
      gate_combn = marker_curr$gate_combn,
      tol = marker_curr$tol,
      noise_sd = NULL,
      bias_uns = marker_curr$bias_uns,
      bw_min = marker_curr$bw_min,
      cp_min = marker_curr$cp_min,
      min_cell = marker_curr$min_cell,
      max_pos_prob_x = max_pos_prob_x,
      gate_quant = gate_quant,
      tol_ctrl = tol_ctrl,
      tol_gate = tol_gate,
      tol_gate_single = tol_gate_single,
      gate_tbl = gate_tbl,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      path_project = path_project
    )

    saveRDS(
      gate_obj$gate_tbl,
      file = file.path(path_project, marker_curr$cut, "gate_tbl.rds")
    )

  })
}

.gate_stats <- function(.data,
                        params = NULL,
                        gate_tbl = NULL,
                        filter_other_cyt_pos = FALSE,
                        gate_name_stats,
                        combn = TRUE,
                        calc_cyt_pos_gates,
                        calc_single_pos_gates,
                        debug,
                        save = TRUE,
                        pop_gate,
                        marker,
                        ind_batch_list,
                        path_project,
                        save_gate_tbl = TRUE) {
  force(.data)
  .get_gate_stats( # nolint
    params = params,
    gate_tbl = gate_tbl,
    filter_other_cyt_pos = filter_other_cyt_pos,
    combn = combn,
    gate_type_cyt_pos_filter =
      if (calc_cyt_pos_gates) "cyt" else "base",
    gate_type_single_pos_filter =
      if (calc_single_pos_gates) "single" else "base",
    gate_type_cyt_pos_calc =
      if (calc_cyt_pos_gates) "cyt" else "base",
    gate_type_single_pos_calc =
      if (calc_single_pos_gates) "single" else "base",
    debug = debug,
    save = save,
    pop_gate = pop_gate,
    chnl = purrr::map_chr(marker, function(x) x$cut),
    ind_batch_list = ind_batch_list,
    .data = .data,
    save_gate_tbl = save_gate_tbl,
    path_project = path_project
  )
}
