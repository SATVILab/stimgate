#' @title Get gates and stats and plots based thereon
#'
#' @description Main function to identify cytokine-positive cells by gating and generate associated statistics and plots.
#'
#' @param path_project character. Path to project directory. Results are saved
#'   here.
#' @param .data GatingSet. GatingSet from which to draw the data.
#' @param batch_list list. List of indices grouped by batch.
#' @param marker list. List where each element specifies the parameters to be
#'   used to gate a given marker including the name of the marker itself.
#' @param pop_gate character vector. Populations for which separate gates are
#'   to be calculated. Default is "root".
#' @param bias_uns numeric. Bias for unstimulated data. Default is NULL.
#' @param bias_uns_factor numeric. Factor for unstimulated bias. Default is 1.
#' @param exc_min logical. Whether to exclude minimum values. Default is TRUE.
#' @param cp_min numeric. Minimum cutpoint value. Default is NULL.
#' @param bw_min numeric. Minimum bandwidth. Default is NULL.
#' @param min_cell numeric. Minimum number of cells required. Default is 1e2.
#' @param max_pos_prob_x numeric. Maximum positive probability x-value. Default is Inf.
#' @param gate_quant numeric vector. Quantiles for gating. Default is c(0.25, 0.75).
#' @param tol_clust numeric. Tolerance for clustering. Default is 1e-7.
#' @param gate_combn character. Method for combining gates. Default is "min".
#' @param marker_settings list. Additional marker settings. Default is NULL.
#' @param calc_cyt_pos_gates logical. Whether to calculate cytokine-positive gates. Default is TRUE.
#' @param calc_single_pos_gates logical. Whether to calculate single-positive gates. Default is FALSE.
#' @param debug logical. Whether to enable debug output. Default is FALSE.
#' @examples
#' \dontrun{
#'   # Basic usage
#'   result <- stimgate_gate(
#'     path_project = "/path/to/project",
#'     .data = gs,
#'     batch_list = list(batch1 = 1:10, batch2 = 11:20),
#'     marker = list(
#'       list(cut = "IL2", tol = 0.5e-8),
#'       list(cut = "TNFa", tol = 0.5e-8)
#'     )
#'   )
#' }
#' @importFrom flowCore exprs<- parameters<-
#' @importFrom stats approx as.formula binomial density glm kmeans median optim predict quantile rnorm sd
#' @importFrom utils read.csv write.csv
#' @importFrom cluster clusGap maxSE
#' @importFrom ggplot2 ggplot aes geom_line geom_smooth geom_vline geom_hline
#' @importFrom dplyr everything
#' @export
stimgate_gate <- function(path_project,
                          .data,
                          batch_list,
                          marker,
                          pop_gate = "root",
                          bias_uns = NULL,
                          bias_uns_factor = 1,
                          exc_min = TRUE,
                          cp_min = NULL,
                          bw_min = NULL,
                          min_cell = 1e2,
                          max_pos_prob_x = Inf,
                          gate_quant = c(0.25, 0.75),
                          tol_clust = 1e-7,
                          gate_combn = "min",
                          marker_settings = NULL,
                          calc_cyt_pos_gates = TRUE,
                          calc_single_pos_gates = FALSE,
                          debug = FALSE) {
  force(.data)
  # capture and force-evaluate the .debug flag into a local .debug object
  .debug <- debug

  if (is.null(names(batch_list))) {
    batch_list <- batch_list |>
      stats::setNames(paste0("batch_", seq_along(batch_list)))
  }

  # get unspecified levels in marker elements
  marker <- .complete_marker_list( # nolint
    marker = marker,
    bias_uns = bias_uns,
    bias_uns_factor = bias_uns_factor,
    exc_min = exc_min,
    .data = .data,
    pop_gate = pop_gate,
    ind_batch_list = batch_list,
    bw_min = bw_min,
    cp_min = cp_min,
    min_cell = min_cell,
    tol_clust = tol_clust,
    max_pos_prob_x = max_pos_prob_x,
    gate_combn = gate_combn,
    marker_settings = marker_settings,
    path_project = path_project,
    .debug = .debug
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
    tol_clust = tol_clust,
    tol_gate_single = tol_clust * 1e-1,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    .debug = .debug
  )

  # cytokine-positive gates
  gate_tbl <- .gate_cyt_pos( # nolint
    marker_list = marker,
    ind_batch_list = batch_list,
    pop_gate = pop_gate,
    .data = .data,
    calc_cyt_pos = calc_cyt_pos_gates,
    .debug = .debug,
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
    tol_clust = tol_clust,
    tol_gate_single = tol_gate_single,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    .debug = .debug,
    gate_tbl = gate_tbl
  )

  message("")
  message("")
  message("")
  message("getting cyt combn frequencies")

  path_dir_stats <- .gate_stats(
    .data = .data,
    params = NULL,
    gate_tbl = NULL,
    filter_other_cyt_pos = FALSE,
    combn = TRUE,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    .debug = .debug,
    save = TRUE,
    pop_gate = pop_gate,
    marker = marker,
    ind_batch_list = batch_list,
    path_project = path_project,
    tol_clust = tol_clust,
    save_gate_tbl = TRUE
  )

  path_project
}

.gate_init <- function(pop_gate,
                       marker,
                       .data,
                       ind_batch_list,
                       path_project,
                       noise_sd,
                       max_pos_prob_x,
                       gate_quant,
                       tol_clust,
                       tol_gate_single,
                       calc_cyt_pos_gates,
                       .debug) {
  # loop over populations
  message("----")
  message("getting base gates")
  message("----")
  message("")
  # loop over markers
  purrr::walk(marker, function(marker_curr) {
    txt <- paste0("chnl: ", marker_curr$chnl_cut)
    message(txt)
    # get gates for each sample within each batch

    gate_obj <- .gate_marker( # nolint
      .data = .data,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate,
      chnl_cut = marker_curr$chnl_cut,
      gate_combn = marker_curr$gate_combn,
      noise_sd = NULL,
      bias_uns = marker_curr$bias_uns,
      exc_min = marker_curr$exc_min,
      bw_min = marker_curr$bw_min,
      cp_min = marker_curr$cp_min,
      min_cell = marker_curr$min_cell,
      max_pos_prob_x = marker_curr$max_pos_prob_x,
      gate_quant = gate_quant,
      tol_clust = tol_clust,
      tol_gate_single = tol_gate_single,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      path_project = path_project,
      .debug = .debug
    )

    path_dir_save <- file.path(path_project, marker_curr$chnl_cut)
    dir.create(path_dir_save, recursive = TRUE, showWarnings = TRUE)

    saveRDS(
      gate_obj$gate_tbl,
      file = file.path(path_dir_save, "gate_tbl_init.rds"
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
                         tol_clust,
                         tol_gate_single,
                         calc_cyt_pos_gates,
                         calc_single_pos_gates,
                         .debug,
                         gate_tbl) {
  # loop over populations
  message("")
  message("")
  message("----")
  message("getting single+ gates")
  message("----")
  message("")
  if (!calc_single_pos_gates) {
    .debug_msg(.debug, "Not gating single-pos gates") # nolint
    purrr::walk(marker, function(marker_curr) {
      saveRDS(
        gate_tbl |>
          dplyr::filter(chnl == marker_curr$chnl_cut) |>
          dplyr::mutate(gate_single = gate),
        file.path(path_project, marker_curr$chnl_cut, "gate_tbl.rds")
      )
    })
    return(invisible(TRUE))
  } else {
    .debug_msg(.debug, "Gating single-pos gates") # nolint
  }
  # loop over markers
  purrr::walk(marker, function(marker_curr) {
    txt <- paste0("chnl: ", marker_curr$chnl_cut)
    message(txt)
    # get gates for each sample within each batch

    gate_obj <- .gate_marker( # nolint
      .data = .data,
      ind_batch_list = ind_batch_list,
      pop_gate = pop_gate_curr,
      chnl_cut = marker_curr$chnl_cut,
      gate_combn = marker_curr$gate_combn,
      tol = marker_curr$tol,
      noise_sd = NULL,
      bias_uns = marker_curr$bias_uns,
      exc_min = marker_curr$exc_min,
      bw_min = marker_curr$bw_min,
      cp_min = marker_curr$cp_min,
      min_cell = marker_curr$min_cell,
      max_pos_prob_x =  marker_curr$max_pos_prob_x,
      gate_quant = gate_quant,
      tol_clust = tol_clust,
      tol_gate_single = tol_gate_single,
      gate_tbl = gate_tbl,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      path_project = path_project,
      .debug = .debug
    )

    saveRDS(
      gate_obj$gate_tbl,
      file = file.path(path_project, marker_curr$chnl_cut, "gate_tbl.rds")
    )

  })
}

.gate_stats <- function(.data,
                        params = NULL,
                        gate_tbl = NULL,
                        filter_other_cyt_pos = FALSE,
                        combn = TRUE,
                        calc_cyt_pos_gates,
                        calc_single_pos_gates,
                        .debug,
                        save = TRUE,
                        pop_gate,
                        marker,
                        ind_batch_list,
                        path_project,
                        tol_clust,
                        save_gate_tbl = TRUE) {
  force(.data)
  .get_stats( # nolint
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
    .debug = .debug,
    save = save,
    pop_gate = pop_gate,
    chnl = purrr::map_chr(marker, function(x) x$chnl_cut),
    ind_batch_list = ind_batch_list,
    .data = .data,
    save_gate_tbl = save_gate_tbl,
    path_project = path_project,
    tol_clust = tol_clust
  )
}
