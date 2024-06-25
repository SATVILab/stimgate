#' @title Get gates and stats and plots based thereon
#'
#' @description
#'
#' @param data GatingSet. GatingSet from which to draw the data.
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
#' @param in_in_batch_gate integer vector. Indices within a given batch for
#'   which gates need to be calculated.
#'
#' # gating parameters
#' @param marker list. List where each element specifies the parameters to be
#'   used to gate a given marker (including the name of the marker itself).
#'   List elements are as follows:
#'   \describe{
#'     \item{cut}{character. Name of channel to get gate for.}
#'     \item{high}{named numeric vector. Names are channel names and values are
#'       thresholds above which a value for a given cell in this channel is
#'       deemed 'high'.}
#'     \item{fdr}{numeric vector. False discovery rates at which the
#'       unstim-based gate is to be calculated. Default is \code{c(0.3, 0.2,
#'       0.1)}.}
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
#'     \item{pop_man}{character vector. Population(s) in GatingHierarchy for
#'       manual count. If \code{pop_man_match_exact} is true, then it is
#'       appended to pop_gate, after the addition of a forward slash (so do not
#'       begin with a forward slash). If \code{pop_man_match_exact} is
#'       \code{FALSE}, then this is passed through as is, and all pops in the
#'       GatingHierarchy are checked for a match with any of these characters.}
#'     \item{pop_man_match_exact}{logical. If \code{TRUE}, then the population
#'       specified above by pop_man is taken to be exact (after appending to
#'       pop_gate, in the manner specified above). If \code{FALSE}, then the
#'       population specified above by pop_man is checked for any match with no
#'       appending. Default is \code{TRUE}.}
#'   }
#' # plotting parameters
#' @param pop_sub character vector. Sub-populations of \code{pop_gate} to plot,
#'   in addition to pop_gate. These populations are appended to pop_gate, after
#'   the addition of a forward slash (so do not begin with a forward slash).

#' @importFrom flowCore exprs<- parameters<-
#' @import rlang
#' @export
stimgate_gate <- function(data,
                          data_name,
                          path_project,
                          pop_gate,
                          pop_sub = NULL,
                          marker,
                          batch_size = NULL,
                          ind_in_batch_gate = NULL,
                          ind_in_batch_uns = NULL,
                          ind_in_batch_lab_vec = NULL,
                          plot = TRUE,
                          bias_uns = NULL,
                          cp_min = NULL,
                          bw_min = NULL,
                          boot_n = NULL,
                          boot_sd = NULL,
                          min_cell = 10,
                          perm_n = NULL,
                          ind_skip = NULL,
                          fcs = NULL,
                          max_pos_prob_x,
                          gate_quant = c(0.25, 0.75),
                          tol_ctrl = NULL,
                          tol_gate = NULL,
                          calc_cyt_pos_gates = TRUE,
                          calc_single_pos_gates = TRUE,
                          gate_name_plot = NULL,
                          gate_name_stats = NULL,
                          tol_gate_single = NULL,
                          debug_stats = FALSE,
                          debug = FALSE,
                          sampleid_lab = NULL,
                          stim_lab = NULL,
                          stats_only = FALSE,
                          plots_only = FALSE) {
  force(data)
  # prep
  debug_stats <- debug_stats || debug

  # get list of batch indices
  ind_batch_list <- .get_ind_batch_list( # nolint
    data = data, batch_size = batch_size,
    fcs = fcs,
    ind_skip = ind_skip
  )

  # get gating indices
  if (is.null(ind_in_batch_gate)) ind_in_batch_gate <- seq_len(batch_size)

  # get unspecified levels in marker elements
  marker <- .complete_marker_list( # nolint
    marker = marker, data_name = data_name,
    pop_gate = pop_gate, cut = cut,
    debug = debug, bias_uns = bias_uns,
    bw_min = bw_min, cp_min = cp_min,
    ind_batch_list = ind_batch_list,
    ind_in_batch_gate = ind_in_batch_gate,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec
  )

  if (stats_only) {
    path_dir_stats <- .get_gate_stats( # nolint
      params = NULL, gate_tbl = NULL,
      filter_other_cyt_pos = FALSE,
      gate_name = gate_name_stats, combn = TRUE,
      gate_type_cyt_pos_calc =
        if (calc_cyt_pos_gates) "cyt" else "base",
      gate_type_single_pos_calc =
        if (calc_single_pos_gates) "single" else "base",
      debug = debug_stats,
      save = TRUE,
      pop_gate = pop_gate,
      chnl = purrr::map_chr(marker, function(x) x$cut),
      ind_in_batch_lab = ind_in_batch_lab_vec,
      ind_in_batch_gate = ind_in_batch_gate,
      data_name = data_name,
      fcs = fcs,
      ind_in_batch_uns = ind_in_batch_uns,
      ind_batch_list = ind_batch_list,
      data = data,
      save_gate_tbl = FALSE,
      sampleid_lab = sampleid_lab,
      stim_lab = stim_lab,
      path_project = path_project
    )
    stats_tbl <- path_dir_stats |> .read_gate_stats() # nolint
    if (!plots_only) {
      return(invisible(TRUE))
    }
  }

  if (plots_only) {
    plot_cp_all( # nolint
      data = data, gate_name = gate_name_plot,
      chnl_base = c("Nd146Di", "Gd156Di", "Ho165Di"),
      params = NULL,
      pop_gate = pop_gate,
      chnl = purrr::map_chr(marker, function(x) x$cut),
      ind_in_batch_lab = ind_in_batch_lab_vec,
      ind_in_batch_gate = ind_in_batch_gate,
      data_name = data_name,
      fcs = fcs,
      ind_in_batch_uns = ind_in_batch_uns,
      ind_batch_list = ind_batch_list,
      path_project = path_project
    )
    return(invisible(TRUE))
  }

  # inital gates
  # ------------------------------

  .gate_init(
    pop_gate = pop_gate,
    marker = marker,
    data = data,
    ind_batch_list = ind_batch_list,
    ind_in_batch_gate = ind_in_batch_gate,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    data_name = data_name,
    path_project = path_project,
    noise_sd = NULL,
    bias_uns = marker$bias_uns,
    bw_min = marker$bw_min,
    cp_min = marker$cp_min,
    plot = plot,
    max_pos_prob_x = max_pos_prob_x,
    gate_quant = gate_quant,
    tol_ctrl = tol_ctrl,
    tol_gate = tol_gate,
    tol_gate_single = tol_gate_single,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    debug = debug,
    pop_sub = pop_sub
  )

  # cytokine-positive gates
  # -----------------------------

  gate_tbl <- .gate_cyt_pos( # nolint
    marker_list = marker,
    pop_gate = pop_gate,
    data = data,
    gate_name = NULL,
    bw_min = bw_min,
    data_name = data_name,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    ind_in_batch_gate = ind_in_batch_gate,
    calc_cyt_pos = calc_cyt_pos_gates,
    debug = debug,
    path_project = path_project
  )

  # single-positive gates
  # -----------------------------

  .gate_single(
    pop_gate = pop_gate,
    marker = marker,
    data = data,
    ind_batch_list = ind_batch_list,
    ind_in_batch_gate = ind_in_batch_gate,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    data_name = data_name,
    path_project = path_project,
    noise_sd = NULL,
    bias_uns = marker$bias_uns,
    bw_min = marker$bw_min,
    cp_min = marker$cp_min,
    plot = plot,
    max_pos_prob_x = max_pos_prob_x,
    gate_quant = gate_quant,
    tol_ctrl = tol_ctrl,
    tol_gate = tol_gate,
    tol_gate_single = tol_gate_single,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    debug = debug,
    pop_sub = pop_sub,
    gate_tbl = gate_tbl,
    gate_name_plot = gate_name_plot,
    gate_name_stats = gate_name_stats,
    stats_only = stats_only,
    plots_only = plots_only
  )

  print("")
  print("")
  print("")
  print("getting cyt combn frequencies")

  path_dir_stats <- .gate_stats(
    data = data,
    params = NULL,
    gate_tbl = NULL,
    filter_other_cyt_pos = FALSE,
    gate_name_stats = gate_name_stats,
    combn = TRUE,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    calc_single_pos_gates = calc_single_pos_gates,
    debug = debug_stats,
    save = TRUE,
    pop_gate = pop_gate,
    marker = marker,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    ind_in_batch_gate = ind_in_batch_gate,
    data_name = data_name,
    fcs = fcs,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_batch_list = ind_batch_list,
    sampleid_lab = sampleid_lab,
    stim_lab = stim_lab,
    path_project = path_project,
    save_gate_tbl = TRUE
  )

  if (plot) {
    print("")
    print("")
    print("")
    print("plotting all 2d plots with final gates")

    plot_cp_all( # nolint
      data = data, gate_name = gate_name_plot,
      chnl_base = c("Nd146Di", "Gd156Di", "Ho165Di"),
      params = NULL,
      pop_gate = pop_gate,
      chnl = purrr::map_chr(marker, function(x) x$cut),
      ind_in_batch_lab = ind_in_batch_lab_vec,
      ind_in_batch_gate = ind_in_batch_gate,
      data_name = data_name,
      fcs = fcs,
      ind_in_batch_uns = ind_in_batch_uns,
      ind_batch_list = ind_batch_list,
      path_project = path_project
    )
  }

  path_dir_stats
}

.gate_init <- function(pop_gate,
                       marker,
                       data,
                       ind_batch_list,
                       ind_in_batch_gate,
                       ind_in_batch_uns,
                       ind_in_batch_lab_vec,
                       data_name,
                       path_project,
                       noise_sd,
                       bias_uns,
                       bw_min,
                       cp_min,
                       plot,
                       max_pos_prob_x,
                       gate_quant,
                       tol_ctrl,
                       tol_gate,
                       tol_gate_single,
                       calc_cyt_pos_gates,
                       debug,
                       pop_sub) {
  # loop over populations
  print("----")
  print("getting base gates")
  print("----")
  print("")
  purrr::walk(pop_gate, function(pop_gate_curr) {
    print(paste0("pop: ", pop_gate_curr))
    # loop over markers
    purrr::walk(marker, function(marker_curr) {
      print(paste0("chnl: ", marker_curr$cut))
      # get gates for each sample within each batch


      gate_obj <- .get_gate_obj( # nolint
        data = data,
        ind_batch_list = ind_batch_list,
        ind_in_batch_gate = ind_in_batch_gate,
        ind_in_batch_uns = ind_in_batch_uns,
        ind_in_batch_lab_vec = ind_in_batch_lab_vec,
        pop_gate = pop_gate_curr,
        data_name = data_name,
        cut = marker_curr$cut,
        high = marker_curr$high,
        gate_combn = marker_curr$gate_combn,
        pop_man_sub = marker_curr$pop_man_sub,
        pop_man_match_exact = marker_curr$pop_man_match_exact,
        tol = marker_curr$tol,
        fdr = marker_curr$fdr,
        noise_sd = NULL,
        bias_uns = marker$bias_uns,
        bw_min = marker$bw_min,
        cp_min = marker$cp_min,
        min_cell = marker_curr$min_cell,
        plot = plot,
        max_pos_prob_x = max_pos_prob_x,
        gate_quant = gate_quant,
        tol_ctrl = tol_ctrl,
        tol_gate = tol_gate,
        tol_gate_single = tol_gate_single,
        pop_sub = pop_sub,
        calc_cyt_pos_gates = calc_cyt_pos_gates,
        path_project = path_project,
        debug = debug
      )

      dir_base <- stimgate_dir_base_create( # nolint
        dir_base_init = path_project,
        params = gate_obj$params
      )
      saveRDS(
        gate_obj$gate_tbl,
        file = file.path(dir_base, "gate_tbl_init.rds")
      )
    })
  })
}

.gate_single <- function(pop_gate,
                         marker,
                         data,
                         ind_batch_list,
                         ind_in_batch_gate,
                         ind_in_batch_uns,
                         ind_in_batch_lab_vec,
                         data_name,
                         path_project,
                         noise_sd,
                         bias_uns,
                         bw_min,
                         cp_min,
                         plot,
                         max_pos_prob_x,
                         gate_quant,
                         tol_ctrl,
                         tol_gate,
                         tol_gate_single,
                         calc_cyt_pos_gates,
                         calc_single_pos_gates,
                         debug,
                         pop_sub,
                         gate_tbl,
                         gate_name_plot,
                         gate_name_stats,
                         stats_only,
                         plots_only) {
  # loop over populations
  print("")
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
  purrr::walk(pop_gate, function(pop_gate_curr) {
    print(paste0("pop: ", pop_gate_curr))
    # loop over markers
    purrr::walk(marker, function(marker_curr) {
      print(paste0("chnl: ", marker_curr$cut))
      # get gates for each sample within each batch

      gate_obj <- .get_gate_obj( # nolint
        data = data,
        ind_batch_list = ind_batch_list,
        ind_in_batch_gate = ind_in_batch_gate,
        ind_in_batch_uns = ind_in_batch_uns,
        ind_in_batch_lab_vec = ind_in_batch_lab_vec,
        pop_gate = pop_gate_curr,
        data_name = data_name,
        cut = marker_curr$cut,
        high = marker_curr$high,
        gate_combn = marker_curr$gate_combn,
        pop_man_sub = marker_curr$pop_man_sub,
        pop_man_match_exact = marker_curr$pop_man_match_exact,
        tol = marker_curr$tol,
        fdr = marker_curr$fdr,
        noise_sd = NULL,
        bias_uns = marker$bias_uns,
        bw_min = marker$bw_min,
        cp_min = marker$cp_min,
        min_cell = marker_curr$min_cell,
        plot = plot,
        max_pos_prob_x = max_pos_prob_x,
        gate_quant = gate_quant,
        tol_ctrl = tol_ctrl,
        tol_gate = tol_gate,
        tol_gate_single = tol_gate_single,
        gate_tbl = gate_tbl,
        calc_cyt_pos_gates = calc_cyt_pos_gates,
        path_project = path_project
      )

      dir_base <- stimgate_dir_base_create( # nolint
        dir_base_init = path_project,
        params = gate_obj$params
      )

      saveRDS(
        gate_obj$gate_tbl,
        file = file.path(dir_base, "gate_tbl.rds")
      )

      # get summary html
      if (plot && FALSE) {
        print("getting plots")
        gate_tbl <- gate_obj[["gate_tbl"]]
        if (!is.null(gate_name_plot)) {
          gate_tbl <- gate_tbl |>
            dplyr::filter(gate_name %in% gate_name_plot) # nolint
        }
        plot_cp( # nolint
          gate_tbl = gate_tbl,
          params = gate_obj[["params"]],
          pop_sub = pop_sub,
          plot = plot
        )
        print("saving results as rmd")
        render_gate( # nolint
          params_knit = gate_obj[["params"]],
          pop_sub = pop_sub,
          path_project = path_project
        )
      }
    })
  })
}

.gate_stats <- function(data,
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
                        ind_in_batch_lab_vec,
                        ind_in_batch_gate,
                        data_name,
                        fcs,
                        ind_in_batch_uns,
                        ind_batch_list,
                        sampleid_lab,
                        stim_lab,
                        path_project,
                        save_gate_tbl = TRUE) {
  force(data)
  .get_gate_stats( # nolint
    params = params,
    gate_tbl = gate_tbl,
    filter_other_cyt_pos = filter_other_cyt_pos,
    gate_name = gate_name_stats,
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
    ind_in_batch_lab = ind_in_batch_lab_vec,
    ind_in_batch_gate = ind_in_batch_gate,
    data_name = data_name,
    fcs = fcs,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_batch_list = ind_batch_list,
    data = data,
    save_gate_tbl = save_gate_tbl,
    sampleid_lab = sampleid_lab,
    stim_lab = stim_lab,
    path_project = path_project
  )
}
