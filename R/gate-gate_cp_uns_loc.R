#' @title Calculate the local fdr-based cut
.get_cp_uns_loc <- function(ex_list,
                            ind_gate,
                            ind_uns,
                            gate_combn,
                            pop_root = NULL,
                            data,
                            bias_uns = 0,
                            noise_sd = NULL,
                            min_bw = 80,
                            cp_min,
                            min_cell,
                            params,
                            plot,
                            path_project,
                            debug = FALSE){

  # get cutpoints for each level of bias
  .get_cp_uns_loc_bias(
    ex_list = ex_list,
    ind_gate = ind_gate,
    ind_uns = ind_uns,
    data = data,
    bias_uns = bias_uns,
    noise_sd = NULL,
    cp_min = cp_min,
    gate_combn = gate_combn,
    min_bw = min_bw,
    min_cell = min_cell,
    params = params,
    plot = plot,
    path_project = path_project,
    debug = debug
  )

}

#' @title Get the unstim-based local fdr-method cutpoint for each level of bias
.get_cp_uns_loc_bias <- function(ex_list,
                                 ind_gate, ind_uns,
                                 data,
                                 bias_uns,
                                 noise_sd,
                                 cp_min,
                                 gate_combn,
                                 min_bw,
                                 min_cell,
                                 params,
                                 plot,
                                 path_project,
                                 debug){

  # get ecdf of uns
  purrr::map(bias_uns, function(bias){
    .debug(debug, "bias_uns", bias)
    #print('getting loc fdr gate for a given bias')

    # get ecdf of uns sample
    # -------------------------------------
    cut_tbl_uns <- .get_cut_list(
      ex_list = ex_list,
      ind = ind_uns,
      exc_min = TRUE,
      bias = bias,
      debug = debug,
      noise_sd = NULL)[[1]]

    # get gates for given level of bias across gate combination methods
    # --------------------------------------
    cp_uns_gate_combn_obj <- .get_cp_uns_loc_gate_combn(
      ex_list =  ex_list,
      cut_uns = cut_tbl_uns,
      ind_uns = ind_uns,
      ind_gate = ind_gate,
      exc_min = TRUE,
      gate_combn = gate_combn,
      cp_min = cp_min,
      min_bw = min_bw,
      min_cell = min_cell,
      params = params,
      plot = plot,
      bias = bias,
      path_project = path_project,
      debug = debug
    )
    # extract and save plots
    # ---------------------------------------

    # save to temp directory
    if (plot && FALSE){
      # add bias label
      cp_uns_plot_list <- cp_uns_gate_combn_obj[['p_list']]
      for(i in seq_along(cp_uns_plot_list)){
        for(j in seq_along(cp_uns_plot_list[[i]])){
          names(cp_uns_plot_list[[i]][[j]]) <- stringr::str_replace(
            names(cp_uns_plot_list[[i]][[j]]),
            "loc",
            paste0("locb", bias)
          )
        }
      }

      dir_save <- file.path(
        tempdir(),
        params$data_name,
        paste0("cp_locb", bias, "_plots")
      )
      if(!dir.exists(dir_save)) {
        dir.create(dir_save, recursive = TRUE)
      }
      saveRDS(
        cp_uns_plot_list,
        file.path(dir_save, paste0(
          names(cp_uns_plot_list),
          ".rds"
        ))
      )
    }

    # extract and add bias label to gates
    # ---------------------------------------
    cp_uns_gate_combn_list <- cp_uns_gate_combn_obj[['cp_uns']]
    names(cp_uns_gate_combn_list) <-
      paste0(names(cp_uns_gate_combn_list), "b", bias)
    cp_uns_gate_combn_list

    #print('done getting loc fd
    # gate across gate combinations for a given bias')

    cp_uns_gate_combn_list

  }) |>
    purrr::flatten()
}

#' @title Get the unstim-based local fdr-method cutpoint for a given bias across gate combination methods
.get_cp_uns_loc_gate_combn <- function(ex_list,
                                       cut_uns,
                                       ind_uns,
                                       ind_gate,
                                       exc_min = TRUE,
                                       gate_combn,
                                       cp_min,
                                       min_bw,
                                       min_cell,
                                       params,
                                       plot,
                                       bias,
                                       path_project,
                                       debug = FALSE){
  .debug(debug, "getting gate_combn")

  #print('getting loc fdr gate across gate combinations')

  # get cutpoints for prejoin gate combination method
  if ('prejoin' %in% gate_combn) {
    .debug(debug, "prejoin")

    # get marker expression for stim samples,
    # join and then sort into descending order
    cut_list_stim <- .get_cut_list(
      ex_list = ex_list,
      ind = setdiff(ind_gate, ind_uns),
      exc_min = TRUE,
      bias = 0,
      noise_sd = NULL,
      debug = debug
    ) |>
      unlist() |>
        sort() |>
        rev() |>
        list()

    # get cutpoints for gate combn method for a range of fdr's
    cp_uns_list_prejoin <- .get_cp_uns_loc_sample(
      cut_stim = cut_list_stim,
      cut_uns = cut_uns,
      cp_min = cp_min,
      ind_uns = ind_uns,
      ind_gate = ind_gate,
      min_bw = min_bw,
      min_cell = min_cell,
      params = params,
      plot = plot,
      bias = bias,
      path_project = path_project
    ) |>
      purrr::map(function(x) list('prejoin' = x))


  } else cp_uns_list_prejoin <- list()

  # get cutpoint if group method is not only prejoin
  non_prejoin_combn_vec <- setdiff(gate_combn, 'prejoin')

  if (length(non_prejoin_combn_vec) > 0){
    .debug(debug, "non-prejoin")

    # get marker expression for stim samples, and sort into descending order
    cut_list_stim <- .get_cut_list(
      ex_list = ex_list,
      ind = setdiff(ind_gate, ind_uns),
      exc_min = TRUE,
      bias = 0,
      noise_sd = NULL
      ) |>
      purrr::map(function(x) x |> dplyr::arrange(desc(expr)))

    cp_uns_list_nonjoin <- .get_cp_uns_loc_sample(
      cut_stim = cut_list_stim,
      cut_uns = cut_uns,
      cp_min = cp_min,
      ind_uns = ind_uns,
      ind_gate = ind_gate,
      min_bw = min_bw,
      min_cell = min_cell,
      params = params,
      plot = plot,
      bias = bias,
      path_project = path_project
    )

    # get list of plots organised
    # ---------------------------
    .debug(debug, "Organising plots")

    indices_name_vec <- paste0(
      names(cp_uns_list_nonjoin[["loc"]]),
      collapse = "_"
    )
    sample_name_vec <- purrr::map_chr(
      names(cp_uns_list_nonjoin[['p_list']]),
      function(ind){
        paste0(ex_list[[ind]]$batch_sh[1],
              "_", ex_list[[ind]]$stim[1])
    })

    p_list_sample_level <- stats::setNames(
      cp_uns_list_nonjoin[["p_list"]],
      sample_name_vec
    )

    p_list <- stats::setNames(
      list(p_list_sample_level), indices_name_vec
    )


    # get list of cutpoints combined in the appropriate way
    # ---------------------------

    .debug(debug, "Combining cutpoints")
    cp_uns_list_nonjoin <- .combine_cp(
      cp = cp_uns_list_nonjoin[["loc"]],
      gate_combn = non_prejoin_combn_vec
    )
  } else cp_uns_list_nonjoin <- list()

  .debug(debug, "Getting unstim cutpoints")

  combined_list <- cp_uns_list_prejoin |>
    append(cp_uns_list_nonjoin)
  cp_uns_list <- purrr::map(
    unique(names(combined_list)),
    function(x){
      .debug(debug, "cutpoint name", paste0(x, collapse = "-"))
      cp_uns_list_prejoin[[x]] |>
        append(cp_uns_list_nonjoin[[x]])
    }) |>
      stats::setNames(unique(names(combined_list)))

  #print('done getting loc fdr gates across gate combinations')
  .debug(debug, "done getting gate_combn")

  list(
    "cp_uns" = list("loc" = cp_uns_list),
    "p_list" = p_list
  )
}

#' @title Get cutpoint for a range of samples given the q-value and fdr
#'
#' @description Calculate the cutpoint for each sample in a batch at a given FDR.
#'
#' @param cut_stim list. List where each element are the marker expression readings
#' of the marker to be cut on for the cells in a sample. Note that the i-th element
#' in \code{cut_stim} must correspond to the i-th element in \code{q_list}, i.e.
#' must be related to the same marker in same cell population from the same blood sample and stimulation.
#' @param fdr numeric. A value between 0 and 1 specifying the false discovery rate
#' the sample should be cut at.
#'
#' @return Numeric vector. A cutpoint for each sample.
.get_cp_uns_loc_sample <- function(cut_stim,
                                   cut_uns,
                                   cp_min,
                                   min_bw,
                                   ind_uns,
                                   ind_gate,
                                   min_cell,
                                   params,
                                   plot = TRUE,
                                   bias,
                                   path_project,
                                   debug = FALSE){
  .debug(debug, "getting loc gate at sample level")
  #print('getting loc fdr gates across sample')

  # get cutpoints for each sample
  cp_uns_loc_obj_list <- purrr::map(seq_along(cut_stim), function(i) {
    .debug(debug, "sample", i)

    # return early if there are too few cells
    too_few_cells_lgl <- .get_cp_uns_loc_sample_check_cell_number(
      cut_stim = cut_stim[[i]], min_cell = min_cell, cut_uns = cut_uns
    )
    if (too_few_cells_lgl) {
      return(.get_cp_uns_loc_sample_out_cell_number(debug))
    }

    # remove any cytokine-positive cells from unstim using gates from
    # sample for which single-positive gates are required
    cut_uns <- .get_cp_uns_loc_sample_uns_rm_cyt_pos(
      debug = debug,
      ex_uns = params$ex_uns,
      gate_tbl = params$gate_tbl,
      cut_stim = cut_stim[[i]],
      gate_name = params$gate_name_curr,
      cut = params$cut,
      calc_cyt_pos_gates = params$calc_cyt_pos_gates,
      bias = bias
    )

    cp_uns_loc_obj <- .get_cp_uns_loc_ind(
      cut_stim = cut_stim[[i]],
      cut_uns = cut_uns,
      min_bw = min_bw,
      cp_min = cp_min,
      min_cell = min_cell,
      params = params,
      plot = plot,
      bias = bias,
      path_project = path_project,
      debug = debug
    )
  })

  # name sample
  # ------------------

  .debug(debug, "Possibly re-using calculated cutpoints")

  # extract vector of cutpoints
  cp_vec <- purrr::map_dbl(
    cp_uns_loc_obj_list,
    function(cp_uns_loc_obj) {
      cp_uns_loc_obj$cp
    }
  )

  # stimulation indices
  ind_stim <- setdiff(ind_gate, ind_uns)

  cp_vec <- purrr::map_dbl(
    cp_uns_loc_obj_list,
    function(cp_uns_loc_obj) {
      cp_uns_loc_obj$cp
    }
  )

  # Repeat cutpoint if it was prejoined
  if (length(cp_vec) == 1 && ind_stim != 1) {
    cp_vec <- stats::setNames(
      rep(cp_vec, length(ind_stim)), ind_stim
    )
  } else {
    # name gate indices if not prejoined
    cp_vec <- stats::setNames(cp_vec, ind_stim)
  }

  # add unstim if unstim gate required
  if(ind_uns %in% ind_gate){
    if(!all(purrr::map_lgl(cp_vec, is.na))){
      cp_vec <- c(cp_vec, stats::setNames(mean(cp_vec, na.rm = TRUE), ind_uns))
    } else cp_vec <- c(cp_vec, NA)
  }

  p_list <- purrr::map(cp_uns_loc_obj_list, function(x) x$p_list) |>
    stats::setNames(ind_stim)

  #print('done getting loc fdr gates for this batch')
  .debug(debug, "done getting loc gate at sample level")

  # collate plots
  list("loc" = cp_vec,
       "p_list" = p_list)
}

.get_cp_uns_loc_sample_check_cell_number <- function(cut_stim, min_cell, cut_uns) {
  nrow(cut_stim[[i]]) < min_cell || nrow(cut_uns) < min_cell
}
.get_cp_uns_loc_sample_out_cell_number <- function(debug) {
  .debug(debug, "Too few cells")
  p_list <- .get_cp_uns_loc_p_list_empty()
  list(p = NA, p_list = p_list)
}

.get_cp_uns_loc_sample_uns_rm_cyt_pos <- function(debug,
                                                  ex_uns,
                                                  gate_tbl,
                                                  cut_stim,
                                                  gate_name,
                                                  cut,
                                                  calc_cyt_pos_gates,
                                                  bias) {
  if (!is.null(gate_tbl)) {
    return(cut_uns)
  }
  .debug(debug, "Removing cytokine-positive cells from unstim")

  # first filter
  gate_tbl_gn_ind <- gate_tbl |>
    # filter to use gates from cut_stim
    dplyr::filter(
      ind == cut_stim$ind[1],
      .data$gate_name == .env$gate_name
    )

  pos_ind_vec_but_single_pos_curr <-
    .get_pos_ind_but_single_pos_for_one_cyt(
      ex = ex_uns,
      gate_tbl = gate_tbl_gn_ind,
      chnl_single_exc = cut,
      chnl = NULL,
      gate_type_cyt_pos = ifelse(calc_cyt_pos_gates, "cyt", "base"),
      gate_type_single_pos = "base"
    )

  ex_uns <- ex_uns[
    !pos_ind_vec_but_single_pos_curr, , drop = FALSE
  ]

  # filter unstim based on this
  .get_cut_list(
    ex_list = stats::setNames(list(ex_uns), ex_uns$ind[1]),
    ind = ex_uns$ind[1],
    exc_min = TRUE,
    bias = bias,
    noise_sd = NULL
  )[[1]]
}

.get_cp_uns_loc_p_list_empty <- function() {
  lapply(1:3, function(x) ggplot()) |>
    stats::setNames(c("p_loc_dens", "p_loc_prob", "p_loc_ctb"))
}

.get_cp_uns_loc_ind_orig <- function(cut_stim, cut_uns) {
  max_dens_x <- .get_cp_uns_loc_ind_max_dens_x(cut_stim)
  list(stim = cut_stim, uns = cut_uns, max_x = max_dens_x)
}
.get_cp_uns_loc_ind <- function(cut_uns,
                                cut_stim,
                                min_bw,
                                cp_min,
                                min_cell,
                                params,
                                plot = TRUE,
                                prob_min = 0.1,
                                bias,
                                path_project,
                                debug = FALSE) {

  .debug(debug, "getting loc gate for single sample")

  # estimate densities for stim and unstim over stim range
  if (.get_cp_uns_loc_check_early(cut_stim, min_cell, cp_min)) {
    return(.get_cp_uns_loc_ind_check_out(
      cp_min, cut_stim, debug, "Too few cells"
    ))
  }

  orig_list <- .get_cp_uns_loc_ind_orig(cut_stim, cut_uns)

  # stop expr being higher than max_x to prevent really far away values creating modes
  cut_stim <- get_cp_uns_loc_set_max_expr(cut_stim, orig_list$max_x)
  cut_uns <- get_cp_uns_loc_set_max_expr(cut_uns, orig_list$max_x)

  # get smoothed probabilities
  data_mod <- .get_cp_uns_loc_get_prob(
    orig_list, cut_stim, cut_uns, debug, min_bw
  )

  # get threshold
  cp <- .get_cp_uns_loc_get_cp(data_mod, orig_list$stim)

  list(cp = cp, p_list = .get_cp_uns_loc_p_list_empty())
}

# initial checks
# ---------------------
.get_cp_uns_loc_ind_check_n_cell <- function(cut_stim, min_cell) {
  nrow(cut_stim) < min_cell
}

.get_cp_uns_loc_ind_max_dens_x <- function(cut_stim) {
  max(cut_stim$expr) - 0.05 * (diff(range(cut_stim$expr)))
}

.get_cp_uns_loc_ind_check_max_x <- function(cut_stim, cp_min) {
  .get_cp_uns_loc_ind_max_dens_x(cut_stim) <= cp_min
}

.get_cp_uns_loc_check_early <- function(cut_stim, min_cell, cp_min) {
  .get_cp_uns_loc_ind_check_n_cell(cut_stim, min_cell) ||
    .get_cp_uns_loc_ind_check_max_x(cut_stim, cp_min)
}

.get_cp_uns_loc_ind_check_out <- function(cp_min,
                                                cut_stim,
                                                debug,
                                                msg) {
  .debug(debug, msg)
  list(
    cp = get_cp_uns_loc_ind_cp_non_loc(cp_min, cut_stim),
    p_list = .get_cp_uns_loc_p_list_empty()
  )
}

.get_cp_uns_loc_ind_cp_non_loc <- function(cp_min, cut_stim) {
  max(cp_min, cut_stim$expr + (max(cut_stim$expr) - min(cut_stim$expr)) / 5)
}

.get_cp_uns_loc_set_max_expr <- function(.data, max_x) {
  .data |> dplyr::mutate(expr = pmin(.data$expr, orig_list$max_x))
}

# get dens_tbl_raw
# -----------------------
.get_cp_uns_loc_get_dens_raw <- function(cut_stim,
                                         cut_uns,
                                         debug,
                                         min_bw) {
  .debug(debug, "Calculating densities")

  dens_list <- .get_cp_uns_loc_get_dens_raw_densities(
    cut_stim, cut_uns, debug, min_bw
  )

  # put raw densities into table
  get_cp_uns_loc_get_dens_raw_tabulate(dens_list = dens_list)
}

.get_cp_uns_loc_get_dens_raw_densities <- function(cut_stim,
                                                   cut_uns,
                                                   debug,
                                                   min_bw) {
  bw <- .get_cp_uns_loc_get_dens_raw_densities_bw(cut_stim, cut_uns, min_bw)
  dens_stim <- .get_cp_uns_loc_get_dens_raw_densities_stim(cut_stim, bw)
  dens_uns <- .get_cp_uns_loc_get_dens_raw_densities_uns(cut_uns, dens_stim, bw)
  list(stim = dens_stim, uns = dens_uns, bw = bw)
}

.get_cp_uns_loc_get_dens_raw_densities_bw_init <- function(.data, min_bw) {
  bw_calc <- try(density(.data, bw  = "SJ")$bw)
  if(inherits(bw_calc, "try-error")) min_bw else bw_calc
}

.get_cp_uns_loc_get_dens_raw_densities_bw <- function(cut_stim,
                                                      cut_uns,
                                                      min_bw) {
  bw_stim <- .get_cp_uns_loc_get_dens_raw_densities_bw_init(cut_stim$expr, min_bw)
  bw_uns <- .get_cp_uns_loc_get_dens_raw_densities_bw_init(cut_uns$expr, min_bw)
  max(bw_stim, min_bw, bw_uns)
}

.get_cp_uns_loc_get_dens_raw_densities_stim <- function(cut_stim, bw) {
  density(cut_stim$expr, bw = bw)
}
.get_cp_uns_loc_get_dens_raw_densities_uns <- function(cut_uns, dens_stim, bw) {
  density(cut_uns$expr, from = min(dens_stim$x), to = max(dens_stim$x), bw = bw)
}

.get_cp_uns_loc_get_dens_raw_tabulate <- function(stim_x,
                                                  stim_y,
                                                  uns_x,
                                                  uns_y) {
  dens_tbl_raw_stim <- tibble::tibble(x_stim = stim_x, y_stim = stim_y)
  dens_tbl_raw_wide <- .get_cp_uns_loc_get_dens_raw_tabulate_uns_interp(
    dens_tbl_raw_stim, uns_x, uns_y
  )
  .get_cp_uns_loc_get_dens_raw_tabulate_format(dens_tbl_raw_wide)
}

.get_cp_uns_loc_get_dens_raw_tabulate_uns_interp <- function(.data, uns_x, uns_y) {
  .data |>
    dplyr::mutate(y_uns = purrr::map_dbl(x_stim,
      function(marker) .interp(val = marker, x = uns_x, y = uns_y)
    ))
}

.get_cp_uns_loc_get_dens_raw_tabulate_format <- function(.data) {
  .data |>
    tidyr::pivot_longer(y_stim:y_uns, names_to = "stim", values_to = "dens") |>
    dplyr::mutate(stim = ifelse(stim == "y_stim", "yes", "no"))
}

#
.get_cp_uns_loc_get_prob_tbl <- function(dens_tbl_raw,
                                         debug,
                                         cp_min,
                                         ex_stim) {
  .debug(debug, "Normalising probabilities")

  # calculate raw and
  # normed probability based on densities for density measurements
  prob_tbl <- get_cp_uns_loc_get_prob_tbl_init(dens_tbl_raw, cp_min)

  # choose probabilities to right of largest peak and
  # with sufficient evidence of a response ito probabilities
  # to be worth taking time smoothing over
  prob_tbl_pos <- .get_cp_uns_loc_prob_tbl_filter(ex_stim, prob_tbl_init, debug)
  list(all = prob_tbl, pos = prob_tbl_pos)
}

.get_cp_uns_loc_get_prob_tbl_init <- function(dens_tbl_raw, cp_min) {
  dens_tbl_raw |>
    tidyr::pivot_wider(
      id_cols = x_stim,
      names_from = stim,
      values_from = dens
      ) |>
    dplyr::mutate(
      prob_stim = 1 - no/yes,
      prob_stim = ifelse(yes == 0 & no == 0, 0, prob_stim),
      prob_stim_norm = pmin(1, prob_stim),
      prob_stim_norm = pmax(0, prob_stim_norm)
    ) |>
    dplyr::filter(x_stim > cp_min)
}

.get_cp_uns_loc_prob_tbl_filter <- function(ex_stim,
                                            prob_tbl,
                                            debug) {
  .debug(debug, "Filtering before smoothing")
  # I don't know why this is ex_stim orig, which
  # excludes the minimum. Shouldn't it be
  # ex_stim? (2024 June 14)
  # changing it to ex_stim_orig.

  # find highest peak (assumed to be the left-most one)
  density_exc_min <- density(ex_stim)
  dens_tbl <- tibble::tibble(x = density_exc_min$x, y = density_exc_min$y)
  peak <- dens_tbl |> dplyr::filter(y == max(y)) |> dplyr::pull("x")

  prob_tbl <- prob_tbl |>
    dplyr::filter(x_stim > peak + 0.02 * diff(range(x_stim)))

  # get range of values for which we'd want to calculate
  # probability:
  # - those cells for which the probability
  # of responding from their position onwards
  # is 0.025 or more
  prob_tbl |>
    dplyr::mutate(
      ge10 = prob_stim_norm >= 0.025, ge10 = cumsum(ge10) > 0
    ) |>
    dplyr::filter(ge10) |>
    dplyr::select(-ge10)
}

.get_cp_uns_loc_get_min_prob_x <- function(prob_tbl_pos) {
  min(prob_tbl_pos$x_stim)
}

.get_cp_uns_loc_check_response <- function(prob_tbl_pos, ex_stim_orig) {
  nrow(prob_tbl_pos) == 0 ||
    max(ex_stim_orig) < .get_cp_uns_loc_get_min_prob_x(prob_tbl_pos)
}

.get_cp_uns_loc_get_data_mod <- function(ex_stim_orig, ex_stim, ex_uns, prob_tbl_list) {
  if (.get_cp_uns_loc_check_response(prob_tbl_list$pos, ex_stim_orig)) {
    return(.get_cp_uns_loc_ind_check_early_out(
      cp_min, ex_stim_orig, debug, "No responding cells"
    ))
  }
  margin <- get_cp_uns_loc_get_data_mod_margin(ex_stim, ex_uns)

  ex_stim_orig |>
    dplyr::filter(expr >
        (min(.get_cp_uns_loc_get_min_prob_x(prob_tbl_list$pos) - margin)))|>
    dplyr:::mutate(prob_smooth = expr)
}

get_cp_uns_loc_get_data_mod_margin <- function(ex_stim, ex_uns) {
  abs(max(diff(ex_stim$expr), diff(ex_uns$expr))) * 0.05
}

# smooth
# ---------------------
.get_cp_uns_loc_get_prob_smooth <- function(data_mod) {
  
  # enough cells to bother smoothing
  if (!.get_cp_uns_loc_get_prob_smooth_check_n_cell(data_mod)) {
    return(.get_cp_uns_loc_get_prob_smooth_check_n_cell_out(data_mod))
  }
  
  # get predictions after smoothing
  pred_vec <- get_cp_uns_loc_get_prob_smooth_actual(data_mod, debug)
  data_mod |> dplyr::mutate(pred = pred_vec)
}

.get_cp_uns_loc_get_prob_smooth_check_n_cell <- function(data_mod) {
  nrow(data_mod) < 10
}

.get_cp_uns_loc_get_prob_smooth_check_n_cell_out <- function(data_mod) {
  data_mod |> dplyr::mutate(pred = prob_smooth - 1e-4)
}

.get_cp_uns_loc_get_prob_smooth_actual <- function(data_mod, debug) {
  fit_1 <- .get_cp_uns_loc_get_prob_smooth_actual_first(data_mod, debug)
  .get_cp_uns_loc_get_prob_smooth_actual_first_response(
    fit_1, data_mod, debug
  )
}

.get_cp_uns_loc_get_prob_smooth_actual_first <- function(data_mod, debug) {
  .debug(debug, "Smoothing I")
  try(
    scam::scam(
      prob_smooth ~ s(expr, bs = "mpi"),
      family = "binomial",
      data = data_mod |>
        dplyr::mutate(
          prob_smooth = pmin(prob_smooth, 0.999),
          prob_smooth = pmax(prob_smooth, 0.001)
        ),
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        devtol.fit = 0.5,
        steptol.fit = 1e-1,
        maxHalf = 5,
        bfgs = list(steptol.bfgs = 1e-1),
        maxit = 1e1
      )
    ),
    silent = TRUE
  )
}

.get_cp_uns_loc_get_prob_smooth_actual_first_response <- function(fit, 
                                                                  data_mod,
                                                                  debug) {
  # return predictions if success
  if (.get_cp_uns_loc_get_prob_smooth_actual_check(fit, data_mod)) {
    .debug(debug, "Smoothed")
    return(.get_cp_uns_loc_get_prob_smooth_actual_response_success(
      fit, data_mod
    )$pred)
  }
  # fit again if not a success
  .get_cp_uns_loc_get_prob_smooth_actual_first_response_failure(
    debug, data_mod
  )
}

.get_cp_uns_loc_get_prob_smooth_actual_check <- function(fit, data_mod) {
  if (inherits(fit, "try-error")) {
    return(FALSE)
  }
  out_list <- .get_cp_uns_loc_get_prob_smooth_actual_response_success(
      fit, data_mod
    )
  !(all(pred_vec > 0.99) ||  mean_abs_error > 0.3)
}

.get_cp_uns_loc_get_prob_smooth_actual_response_success <- function(fit, # nolint
                                                                          data_mod) { # nolint
    pred_vec <- predict(fit, type = 'response')
    mean_abs_error <- mean(abs(pred_vec - data_mod$prob_smooth))
    list("pred" = pred_vec, "mean_abs_error" = mean_abs_error)
  }

.get_cp_uns_loc_get_prob_smooth_actual_first_response_failure <- function(debug,# nolint
                                                                          data_mod) {
    fit_2 <- .get_cp_uns_loc_get_prob_smooth_actual_second(data_mod, debug)
    if (.get_cp_uns_loc_get_prob_smooth_actual_check(fit_2, data_mod)) {
      .debug(debug, "Smoothed")
      return(.get_cp_uns_loc_get_prob_smooth_actual_response_success(
          fit_2, data_mod
        )$pred)
    }
    .get_cp_uns_loc_get_prob_smooth_actual_third(data_mod, debug)
  }

.get_cp_uns_loc_get_prob_smooth_actual_second <- function(data_mod, debug) {
  .debug(debug, "Smoothing II")
  try(
    scam::scam(
      prob_smooth ~ s(expr, bs = "micv"),
      family = "binomial",
      data = data_mod,
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        devtol.fit = 0.01
        )),
        silent = TRUE
    )
}

.get_cp_uns_loc_get_prob_smooth_actual_third <- function(data_mod) {
  .debug(debug, "Failed to smooth")
  data_mod$prob_smooth - 0.0001
}

# get probabilities
.get_cp_uns_loc_get_prob <- function(orig_list,
                                     cut_stim,
                                     cut_uns,
                                     debug,
                                     min_bw) {

  # get raw densities
  dens_tbl_raw <- .get_cp_uns_loc_get_dens_raw(
    cut_stim, cut_uns, debug, min_bw
  )

  # get probabilities
  prob_tbl_list <- .get_cp_uns_loc_get_prob_tbl(
    dens_tbl_raw, debug, cp_min, ex_stim
  )

  # get data to smooth over
  data_mod <- .get_cp_uns_loc_get_data_mod(
    orig_list$stim, ex_stim, ex_uns, prob_tbl_list
  )

  # smooth
  .get_cp_uns_loc_get_prob_smooth(data_mod)
}

# get cp
.get_cp_uns_loc_get_cp <- function(data_mod,
                                   ex_stim_orig,
                                   debug) {
  if (!inherits(data_mod, "data.frame")) {
    return(data_mod)
  }

  data_threshold <- .get_cp_uns_loc_get_cp_data_threshold(
    data_mod, ex_stim_orig
  )
  cp <- .get_cp_uns_loc_get_cp_actual(data_threshold)
  .debug(debug, "Completed loc gate for single sample")
  cp
}

.get_cp_uns_loc_get_cp_data_threshold <- function(data_mod,
                                                  ex_stim_orig) {
  data_count <- .get_cp_uns_loc_get_cp_data_threshold_count(data_mod)
  prob_bs_est <- .get_cp_uns_loc_get_cp_data_threshold_prop_bs_est(data_count)
  .get_cp_uns_loc_get_cp_data_threshold_actual(
    data_count, prob_bs_est, ex_stim_orig
  )
}

.get_cp_uns_loc_get_cp_data_threshold_count <- function(data_mod) {
  data_mod |>
    dplyr::filter(expr >= min(.data$expr)) |>
    dplyr::arrange(expr) |>
    dplyr::mutate(n_row = seq_len(dplyr::n())) |>
    dplyr::filter(cumsum(pred > prob_smooth) != n_row) |>
    dplyr::select(-n_row)
}

.get_cp_uns_loc_get_cp_data_threshold_prop_bs_est <- function(data_count) {
  sum(data_count$pred) / nrow(orig_list$stim)
}

.get_cp_uns_loc_get_cp_data_threshold_actual <- function(data_count,
                                                         prop_bs_est) {
  data_count |>
    dplyr::arrange(desc(expr)) |>
    dplyr::mutate(count_stim = seq_len(dplyr::n())) |>
    dplyr::mutate(
      prop_stim = count_stim/nrow(orig_list$stim),
      prop_uns = purrr::map_dbl(expr, function(x){
        sum((cut_uns$expr - bias) >= x) / nrow(orig_list$uns)
      }),
      prop_bs = prop_stim - prop_uns,
      prop_bs_diff = prop_bs - prop_bs_est
    )
}

.get_cp_uns_loc_get_cp_actual <- function(data_threshold) {
  data_threshold |>
    dplyr::filter(abs(prop_bs_diff) ==  min(abs(prop_bs_diff))) |>
    dplyr::slice(1) |>
    dplyr::pull(expr)
}

#' @title Plot gating plots for local fdr method
#'
#' @description
#' Plot the following three plots:
#' \describe{
#'   \item{p_loc_dens}{Density of the stim and unstim (pos and neg) distributions.}
#'   \item{p_loc_prob}{Raw, normalised and smoothed probability of positivty by marker value.}
#'   \item{p_loc_ctb}{Contribution to total count by binned marker values}.
#' }
#'
#' @param dens_tbl dataframe. Dataframe equal to \code{dens_tbl_raw}
#' in the \code{.get_cp_uns_loc_ind} function.
#' @inheritParams plot_cp
#' @param prob_mod object of class \code{glm}. Models probability of positivity
#' as a function of \code{x_stim}. Created within \code{.get_cp_uns_loc_ind}
#' function.
#' @param prob_tbl dataframe. Dataframe contains columns x_stim,
#' prob_stim, prob_stim_norm and pred, where x_stim are marker values
#' and prob_stim, prob_stim_norm and pred are the raw, normalised and
#' smoothed probabilities of posititivity, respectively. Dataframe is
#' produced appropriately within \code{.get_cp_uns_loc_ind} function.
#' @param cut_stim numeric vector. Each element is an observed value of the marker
#' in the stim condition.
#' @param min_x_pos_prob numeric. Minimum value for \code{cut_stim} such that
#' the probability of positivity is greater than the value specified by the
#' \code{prob_min} parameter in the \code{.get_cp_uns_loc_ind} function.
#'
#' @return A list with elements named p_loc_dens, p_loc_prob and
#' p_loc_ctb, where each element is a ggplot2 plot.
#'
#' @examples
#' For example of use, set \code{debugonce(:::.get_cp_uns_loc_ind)}
#' before running \code{::gate}. Step through the debuggec function
#' until near the end.
.plot_cp_loc_fdr <- function(dens_tbl,
                             params,
                             prob_tbl,
                             pred_tbl,
                             cut_stim,
                             sample,
                             min_x_pos_prob,
                             path_project){

  dir_base <- stimgate_dir_base_create(
    params = params, dir_base_init = path_project
  )
  dir_base <- file.path(dir_base, "gating_plots", "gate")
  dir_len <- stringr::str_length(dir_base) + 10
  underscore_loc <- stringr::str_locate(sample, "_")[1,"start"][[1]]
  subject_visittype <- stringr::str_sub(sample, end = underscore_loc - 1)
  stim <- stringr::str_sub(sample, underscore_loc + 1)
  subjectid_visittype_len <- 255 - dir_len - 4 - 4
  subjectid_visittype_len <- min(
    stringr::str_length(subject_visittype), subjectid_visittype_len
  )
  subject_visittype <- stringr::str_sub(
    subject_visittype,
    end = subjectid_visittype_len
  )
  fn <- paste0(subject_visittype, "-", stim, ".png")

  # Plot densities
  # -------------------------------
  p_loc_dens <- ggplot(dens_tbl |>
                         dplyr::filter(x_stim > min(min_x_pos_prob, 100)),
                       aes(x = x_stim, y = dens, linetype = stim)) +
    cowplot::theme_cowplot(font_size = 20) +
    geom_line() +
    scale_linetype_manual(values = c("yes" = "solid",
                                     "no" = "dotted"),
                          labels = c("yes" = "stim",
                                     "no" = "unstim")) +
    labs(x = params$chnl_lab[params$cut], y = "Density") +
    theme(legend.title = element_blank())

  dir_save <- file.path(dir_base, "p_loc_dens")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  cowplot::ggsave2(filename = file.path(file.path(dir_save, fn)),
                   height = 10, width = 10)

  # Plot probabilities of being stim-induced
  # -------------------------------

  plot_tbl_prob <- prob_tbl |>
    tidyr::pivot_longer(prob_stim:prob_stim_norm,
                 names_to = "prob_type",
                 values_to = "prob") |>
    dplyr::select(x_stim, prob_type, prob)

  plot_tbl_pred <- pred_tbl |>
    dplyr::rename(x_stim = cut_stim,
           prob = pred) |>
    dplyr::mutate(prob_type = 'pred') |>
    dplyr::select(x_stim, prob_type, prob )

  plot_tbl <- plot_tbl_prob |>
    dplyr::bind_rows(plot_tbl_pred)
  p_loc_prob <- ggplot(plot_tbl,
                       aes(x = x_stim, y = prob, linetype = prob_type)) +
    cowplot::theme_cowplot(font_size = 20) +
    cowplot::background_grid(major = 'y', minor = 'y') +
    geom_hline(yintercept = 0) +
    geom_line(size = 2) +
    scale_linetype_manual(values = c("pred" = "solid",
                                     "prob_stim_norm" = "dotted",
                                     "prob_stim" = "dashed"),
                          labels = c("pred" = "Smoothed - Model",
                                     "prob_stim" = "Raw",
                                     "prob_stim_norm" = "Smoothed - Raw")) +
    lims(y = c(-1, 1)) +
    theme(legend.title = element_blank()) +
    labs(x = params$chnl_lab[params$cut],
         y = "Probability of being stimulation-induced")

  dir_save <- file.path(dir_base, "p_loc_prob")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  cowplot::ggsave2(filename = file.path(file.path(dir_save, fn)),
                   height = 10, width = 10)

  # Plot contribution to count by bin
  # -------------------------------

  # plot binned contributions to total count by range
  # new_pred_tbl <- tibble::tibble(x_stim = cut_stim[cut_stim >= min_x_pos_prob])
  # new_pred_vec <- predict(prob_mod, newdata = new_pred_tbl, type = 'response')
  # new_pred_tbl <- new_pred_tbl |> dplyr::mutate(pred = new_pred_vec)

  #  n_bin <- hist(new_pred_tbl$x_stim, plot = FALSE)$mids |> length
  n_bin <- hist(pred_tbl$cut_stim, plot = FALSE)$mids |> length()
  if(n_bin == 1){
    bin_tbl <- pred_tbl |>
      dplyr::mutate(bin = paste0("(", cut_stim, ",", cut_stim, "]"),
             ctb = pred)
  } else{
    bin_tbl <- pred_tbl |>
      dplyr::mutate(bin = cut(.data$cut_stim, breaks = n_bin)) |>
      dplyr::group_by(bin) |>
      dplyr::summarise(ctb = sum(pred))
  }

  #' @title Get the midpoint of a bin returned by base::cut
  .get_bin_mid <- function(bin){
    purrr::map_dbl(bin, function(bin_curr){
      comma_loc <- stringr::str_locate(bin_curr, ",")[1,"start"][[1]]
      num_1 <- stringr::str_sub(bin_curr, 2, comma_loc - 1) |>
        as.numeric()
      num_2 <- stringr::str_sub(
        bin_curr, comma_loc + 1, stringr::str_length(bin_curr) - 1
        ) |>
        as.numeric()
      mean(c(num_1, num_2))
    })
  }

  bin_tbl <- bin_tbl |>
    dplyr::mutate(bin_mid = .get_bin_mid(bin))

  # plot of contribution per bin
  p_loc_ctb <- ggplot(bin_tbl,
                      aes(x = bin_mid, y = ctb)) +
    geom_bar(stat = 'identity',
             col = 'gray25',
             fill = 'gray85')

  dir_save <- file.path(dir_base, "p_loc_ctb")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  cowplot::ggsave2(filename = file.path(file.path(dir_save, fn)),
                   height = 10, width = 10)

  #list(p_loc_dens = p_loc_dens, p_loc_prob = p_loc_prob,
  #     p_loc_ctb = p_loc_ctb)
  list()

}

.get_prob_mod <- function(prob_tbl, nrow_level){

  prob_mod <- try(scam::scam(prob_stim_norm ~ s(x_stim, bs = "mpi"),
                         family = "binomial",
                         data = prob_tbl |>
                           dplyr::filter(prob_stim_norm > 0.25)))

  if(class(prob_mod) == 'try-error'){
    prob_mod <- try(scam::scam(prob_stim_norm ~ s(x_stim, bs = "micx"),
                               family = "binomial",
                               data = prob_tbl))
  } else return(prob_mod)

  if(class(prob_mod) == 'try-error'){
    prob_mod <- try(mgcv::gam(prob_stim_norm ~ s(x_stim),
                               family = "binomial",
                               data = prob_tbl))
  }

  prob_mod
}

.get_cell_specific_response_probs <- function(prob_tbl, nrow_level, params,
                                              prob_min, cut_stim){

  # =====================
  # Fit model
  # =====================

  # fit a model to get smoothed probabilities
  prob_mod <- .get_prob_mod(prob_tbl, nrow_level = nrow_level)

  # if model fails, then return it
  if('ultimate error' %in% class(prob_mod)){

    out <- list(pred_tbl = NA,
                prob_mod = NA)

    # out <- FALSE
    # class(out) <- 'ultimate error'
    return(out)
  }

  # =====================
  # Get predictions
  # =====================

  # get fitted values
  pred_vec <- fitted.values(prob_mod)

  # add them to table
  pred_tbl <- prob_tbl |>
    dplyr::mutate(pred = pred_vec)

  out <- list(pred_tbl = pred_tbl,
              prob_mod = prob_mod)

  out
}


