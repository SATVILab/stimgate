#' @title Check that parameters for each marker for which a gate is required are complete
#'
#' @inheritParams gate # marker
#' @data_name 'gs_cytof' or `gs_proto`. Name of dataset to be gated in R environment.
#'
#' @return

.complete_marker_list <- function(marker,
                                  bias_uns,
                                  bias_uns_factor,
                                  exc_min,
                                  .data,
                                  pop_gate,
                                  ind_batch_list,
                                  bw_min,
                                  cp_min,
                                  min_cell,
                                  tol_clust,
                                  max_pos_prob_x,
                                  gate_combn,
                                  marker_settings,
                                  path_project,
                                  debug) {
  marker_settings_common <- list(
    bias_uns = bias_uns, bias_uns_factor = bias_uns_factor,
    exc_min = exc_min, cp_min = cp_min, bw_min = bw_min,
    min_cell = min_cell, tol_clust = tol_clust, gate_combn = gate_combn,
    max_pos_prob_x = max_pos_prob_x
  )
  marker_list <- purrr::map(marker, function(marker_curr) {
   .complete_marker_list_ind(
      marker = marker_curr,
      marker_settings_common = marker_settings_common,
      marker_settings_spec = list(chnl_cut = marker_curr) |>
        append(marker_settings[[marker_curr]]),
      .data = .data,
      pop_gate = pop_gate,
      .debug = .debug,
      ind_batch_list = ind_batch_list
    )
  })

  .complete_marker_list_save(
    marker_list = marker_list,
    path_project = path_project
  )

  marker_list
}

.complete_marker_list_ind <- function(marker_settings_common,
                                      marker_settings_spec,
                                      marker,
                                      .data,
                                      pop_gate,
                                      .debug,
                                      ind_batch_list) {

  marker_settings <- .complete_marker_list_add_common(
    marker_settings_common = marker_settings_common,
    marker_settings = marker_settings_spec
  )

  marker_settings$bias_uns <- .complete_marker_list_bias_uns(
    bias_uns = marker_settings$bias_uns,
    bias_uns_factor = marker_settings$bias_uns_factor,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = marker,
    .debug = .debug,
    ind_batch_list = ind_batch_list
  )

  marker_settings$bw_min <- .complete_marker_list_min_bw(
    bw_min = marker_settings$bw_min,
    bias_uns = max(marker_settings$bias_uns)
  )

  marker_settings$cp_min <- .complete_marker_list_cp_min(
    cp_min = marker_settings$cp_min,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = marker,
    .debug = .debug,
    ind_batch_list = ind_batch_list
  )

  marker_settings
}

.complete_marker_list_add_common <- function(marker_settings_common,
                                             marker_settings) {
  marker_settings |>
    append(marker_settings_common[
      setdiff(names(marker_settings_common), names(marker_settings))
    ])
}

.complete_marker_list_bias_uns <- function(bias_uns,
                                           bias_uns_factor,
                                           .data,
                                           pop_gate,
                                           chnl_cut,
                                           .debug,
                                           ind_batch_list) {
  if (!is.null(bias_uns)) {
    return(bias_uns)
  }
  .debug_msg(.debug, "calculating bias_uns automatically") # nolint
  mean_range <- .complete_marker_list_bias_uns_get_mean_range(
    ind_batch_list = ind_batch_list,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = chnl_cut
  )
  (mean_range / 12 * bias_uns_factor) |> signif(3)
}

.complete_marker_list_bias_uns_get_mean_range <- function(ind_batch_list,
                                                          .data,
                                                          pop_gate,
                                                          chnl_cut) {
  purrr::map(
    seq_len(min(2, length(ind_batch_list))),
    function(i) {
      ex_list <- .get_ex_list( # nolint
        .data = .data,
        ind_batch = ind_batch_list[[i]],
        pop = pop_gate,
        chnl_cut,
        batch = names(ind_batch_list)[i]
      )
      purrr::map_dbl(ex_list, function(ex) {
        abs(
          diff(quantile(.get_cut(ex)[.get_cut(ex) > min(.get_cut(ex))], c(0.99, 0.01)),
            na.rm = TRUE
          )
        )[[1]]
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
}


.complete_marker_list_min_bw <- function(bw_min, bias_uns) {
  if (!is.null(bw_min)) {
    return(bw_min)
  }
  bias_uns * 2.25
}

.complete_marker_list_cp_min <- function(cp_min,
                                         .data,
                                         pop_gate,
                                         chnl_cut,
                                         .debug,
                                         ind_batch_list) {
  if (!is.null(cp_min)) {
    return(cp_min)
  }
  .debug_msg(.debug, "calculating cp_min automatically") # nolint
  purrr::map(
    seq_len(min(2, length(ind_batch_list))),
    function(i) {
      ex_list <- .get_ex_list( # nolint
        .data = .data,
        ind_batch = ind_batch_list[[i]],
        pop = pop_gate,
        chnl_cut,
        batch = names(ind_batch_list)[i]
      )
      purrr::map_dbl(ex_list, function(ex) {
        median(.get_cut(ex)[.get_cut(ex) > min(.get_cut(ex))], na.rm = TRUE)[[1]]
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
}



#' @title Get all cp type names
#'
#' @inheritParams get_cp #fdr
#'
#' @return Character vector, where each element is name of a
#' cutpoint, and all elements together represent names of all cutpoints.
.get_full_cp_type_vec <- function(fdr) {
  # Get cutpoint names for unstim-based cuts
  # cp_name_vec_uns <- .get_cp_uns_name_vec(fdr)
  # cp_name_vec_uns_root <- paste0(cp_name_vec_uns, 'root')

  # output all cutpoint names
  c(
    "man", "tg", "dcp", "midp", "scp",
    "uns", "unsr", "loc"
  )
}

.complete_marker_list_save <- function(marker_list, path_project) {
  path_save <- file.path(path_project, "marker_list.rds")
  if (file.exists(path_save)) {
    invisible(file.remove(path_save))
  }
  if (!dir.exists(path_project)) {
    dir.create(path_project, recursive = TRUE)
  }
  saveRDS(
    marker_list,
    file = path_save
  )
}