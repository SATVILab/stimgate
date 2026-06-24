# Calculate local FDR-based cutpoint
#' @keywords internal
.get_cp_uns_loc <- function(
  ex_list,
  .data,
  chnl_settings,
  stage,
  path_project
) {
  # get cutpoints for each level of bias
  .get_cp_uns_loc_bias(
    ex_list = ex_list,
    .data = .data,
    chnl_settings = chnl_settings,
    stage = stage,
    path_project = path_project
  )
}

# Get unstim-based local FDR cutpoint for each bias level
#' @keywords internal
.get_cp_uns_loc_bias <- function(
  ex_list,
  .data,
  chnl_settings,
  stage,
  path_project
) {
  # get ecdf of uns
  purrr::map(chnl_settings$bias_uns, function(bias) {
    .debug("bias_uns", bias) # nolint

    ex_list_prep <- .prepare_data_with_bias_and_noise(
      ex_list = ex_list,
      bias = bias,
      noise_sd = NULL,
      exc_min = chnl_settings$exc_min
    )

    # get gates for given level of bias across gate combination methods
    # --------------------------------------
    cp_uns_gate_combn_obj <- .get_cp_uns_loc_gate_combn(
      ex_list_orig = ex_list_prep[["ex_list_orig"]],
      ex_list_no_min = ex_list_prep[["ex_list_no_min"]],
      ex_tbl_uns_bias = ex_list_prep[["ex_tbl_uns_bias"]],
      chnl_settings = chnl_settings,
      bias = bias,
      path_project = path_project,
      stage = stage
    )

    # extract and add bias label to gates
    # ---------------------------------------
    .get_cp_uns_loc_gate_label(cp_uns_gate_combn_obj, bias)
  }) |>
    purrr::flatten()
}

#' @keywords internal
.get_cp_uns_loc_gate_label <- function(cp_uns_gate_combn_obj, bias) {
  cp_uns_gate_combn_list <- cp_uns_gate_combn_obj[["cp_uns"]]
  names(cp_uns_gate_combn_list) <-
    paste0(names(cp_uns_gate_combn_list))
  cp_uns_gate_combn_list
}

#' @keywords internal
.prepare_data_with_bias_and_noise <- function(
  ex_list,
  bias,
  exc_min,
  noise_sd
) {
  # rename `cut` column` to `expr`
  # -------------------------------------
  ex_list_orig <- .prepare_ex_list_with_bias_and_noise(
    # exc_min should always be FALSE here, as we're trying
    # to keep all the original data
    ex_list = ex_list,
    ind = names(ex_list),
    exc_min = FALSE,
    bias = 0,
    noise_sd = NULL
  ) |>
    .arrange_samples_by_desc_expr()

  # separate stim and uns samples, rename
  # `cut` column to `expr` and exclude min
  # values
  # -------------------------------------
  ex_list_no_min <- .prepare_ex_list_with_bias_and_noise(
    ex_list = ex_list,
    ind = names(ex_list),
    exc_min = exc_min,
    bias = 0,
    noise_sd = NULL
  ) |>
    .arrange_samples_by_desc_expr()

  # adjust expression in unstim,
  # applying bias, excluding the min val
  # and /or adding noise
  # -------------------------------------
  ex_list_bias <- .prepare_ex_list_with_bias_and_noise(
    ex_list = ex_list,
    ind = names(ex_list)[length(ex_list)],
    exc_min = exc_min,
    bias = bias,
    noise_sd = NULL
  ) |>
    .arrange_samples_by_desc_expr()
  ex_tbl_uns_bias <- ex_list_bias[[1]]

  list(
    "ex_list_orig" = ex_list_orig,
    "ex_list_no_min" = ex_list_no_min,
    "ex_tbl_uns_bias" = ex_tbl_uns_bias
  )
}

# Get the unstim-based local fdr-method
# cutpoint for a given bias across gate combination methods
#' @keywords internal
.get_cp_uns_loc_gate_combn <- function(
  ex_list_orig,
  ex_list_no_min,
  ex_tbl_uns_bias,
  chnl_settings,
  bias,
  path_project,
  stage
) {
  .debug("getting gate_combn") # nolint

  # get cutpoints for prejoin gate combination method
  cp_uns_list_prejoin <- .get_cp_uns_loc_gate_combn_prejoin(
    ex_list_no_min = ex_list_no_min,
    ex_list_orig = ex_list_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    chnl_settings = chnl_settings,
    bias = bias,
    path_project = path_project,
    stage = stage
  )

  # get cutpoint for non-prejoin grouping methods
  cp_uns_list_prejoin_non <- .get_cp_uns_loc_gate_combn_prejoin_non(
    non_prejoin_combn = setdiff(chnl_settings$gate_combn, "prejoin"),
    ex_list_no_min_stim = ex_list_no_min[-length(ex_list_no_min)],
    ex_list_orig = ex_list_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    chnl_settings = chnl_settings,
    bias = bias,
    path_project = path_project,
    stage = stage
  )

  # merge above two lists
  cp_uns_list <- .get_cp_uns_loc_gate_combn_merge(
    cp_uns_list_prejoin,
    cp_uns_list_prejoin_non,
    stage
  )

  list(
    "cp_uns" = list("loc" = cp_uns_list),
    "p_list" = list()
  )
}

#' @keywords internal
.get_cp_uns_loc_gate_combn_merge <- function(
  cp_uns_list_prejoin,
  cp_uns_list_prejoin_non,
  stage
) {
  .debug("done getting gate_combn") # nolint

  combined_list <- cp_uns_list_prejoin |>
    append(cp_uns_list_prejoin_non)
  purrr::map(
    unique(names(combined_list)),
    function(x) {
      .debug("cutpoint name", paste0(x, collapse = "-")) # nolint
      cp_uns_list_prejoin[[x]] |>
        append(cp_uns_list_prejoin_non[[x]])
    }
  ) |>
    stats::setNames(unique(names(combined_list)))
}

# --------------------------------
# gate using prejoined .data
# --------------------------------
#' @keywords internal
.get_cp_uns_loc_gate_combn_prejoin <- function(
  ex_list_no_min,
  ex_list_orig,
  ex_tbl_uns_bias,
  chnl_settings,
  bias,
  path_project,
  stage
) {
  if (!"prejoin" %in% chnl_settings$gate_combn) {
    return(.get_cp_uns_loc_gate_combn_prejoin_not())
  }
  .get_cp_uns_loc_gate_combn_prejoin_actual(
    ex_list_no_min = ex_list_no_min,
    ex_list_orig = ex_list_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    chnl_settings = chnl_settings,
    bias = bias,
    path_project = path_project,
    stage = stage
  )
}

#' @keywords internal
.get_cp_uns_loc_gate_combn_prejoin_not <- function() {
  list("cp" = list(), "p_list" = list())
}

#' @keywords internal
.prepare_data_for_prejoin <- function(ex_list_orig, ex_list_no_min) {
  ex_list_no_min <- .prepare_data_for_prejoin_ind(
    ex_list = ex_list_no_min
  )
  ex_list_orig <- .prepare_data_for_prejoin_ind(
    ex_list = ex_list_orig
  )

  list(
    "ex_list_orig" = ex_list_orig,
    "ex_list_no_min" = ex_list_no_min
  )
}

#' @keywords internal
.prepare_data_for_prejoin_ind <- function(ex_list) {
  ind_uns <- names(ex_list)[length(ex_list)]
  ex_tbl_stim <- ex_list[seq_len(length(ex_list) - 1)] |>
    dplyr::bind_rows()
  ex_tbl_stim <- ex_tbl_stim[order(.get_cut(ex_tbl_stim)), ]
  list(
    ex_tbl_stim,
    ex_list[[length(ex_list)]]
  ) |>
    stats::setNames(c(names(ex_list)[1], ind_uns))
}

#' @keywords internal
.get_cp_uns_loc_gate_combn_prejoin_actual <- function(
  ex_list_no_min,
  ex_list_orig,
  ex_tbl_uns_bias,
  chnl_settings,
  bias,
  path_project,
  stage
) {
  .debug("prejoin") # nolint

  # get marker expression for stim samples,
  # join and then sort into descending order
  ex_list_prejoin <- .prepare_data_for_prejoin(
    ex_list_orig = ex_list_orig,
    ex_list_no_min = ex_list_no_min
  )

  # get cutpoints for gate combn method for a range of fdr's
  .get_cp_uns_loc_sample(
    ex_list_orig = ex_list_prejoin[["ex_list_orig"]],
    ex_list_no_min_stim = ex_list_prejoin[["ex_list_no_min"]],
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    chnl_settings = chnl_settings,
    bias = bias,
    stage = stage,
    path_project = path_project
  ) |>
    purrr::map(function(x) list("prejoin" = x))
}

# --------------------------------
# gate each sample individually
# --------------------------------
#' @keywords internal
.get_cp_uns_loc_gate_combn_prejoin_non <- function(
  non_prejoin_combn,
  ex_list_no_min_stim,
  ex_list_orig,
  ex_tbl_uns_bias,
  chnl_settings,
  bias,
  path_project,
  stage
) {
  if (length(non_prejoin_combn) == 0L) {
    return(
      .get_cp_uns_loc_gate_combn_prejoin_non_not()
    )
  }
  .get_cp_uns_loc_gate_combn_prejoin_non_actual(
    ex_list_no_min_stim = ex_list_no_min_stim,
    ex_list_orig = ex_list_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    chnl_settings = chnl_settings,
    bias = bias,
    non_prejoin_combn = non_prejoin_combn,
    path_project = path_project,
    stage = stage
  )
}

#' @keywords internal
.get_cp_uns_loc_gate_combn_prejoin_non_not <- function() {
  list("cp" = list(), "p_list" = list())
}

#' @keywords internal
.get_cp_uns_loc_gate_combn_prejoin_non_actual <- function(
  ex_list_no_min_stim,
  ex_list_orig,
  ex_tbl_uns_bias,
  chnl_settings,
  bias,
  path_project,
  stage,
  non_prejoin_combn
) {
  .debug("non-prejoin") # nolint
  cp_uns_list_nonjoin <- .get_cp_uns_loc_sample(
    ex_list_orig = ex_list_orig,
    ex_list_no_min_stim = ex_list_no_min_stim,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    chnl_settings = chnl_settings,
    bias = bias,
    stage = stage,
    path_project = path_project
  )

  cp_uns_list_nonjoin <- .get_cp_uns_loc_gate_combn_prejoin_non_actual_combn(
    stage,
    cp_uns_list_nonjoin,
    non_prejoin_combn
  )
  list("cp" = cp_uns_list_nonjoin, "p_list" = list())
}

#' @keywords internal
.arrange_samples_by_desc_expr <- function(ex_list) {
  # arrange in descending order of expression
  ex_list |>
    purrr::map(function(x) x[order(.get_cut(x)), ]) # nolint
}

#' @keywords internal
.get_cp_uns_loc_gate_combn_prejoin_non_actual_combn <- function(
  stage,
  cp_uns_list_nonjoin, # nolint
  non_prejoin_combn_vec
) {
  .debug("Combining cutpoints") # nolint
  .combine_cp(
    cp = cp_uns_list_nonjoin[["loc"]],
    gate_combn = non_prejoin_combn_vec
  )
}

# ------------------------------------------
# get cutpoints for a range of samples, and then individual samples
# ------------------------------------------

# Get cutpoint for a range of samples given the q-value and fdr
#' @keywords internal
.get_cp_uns_loc_sample <- function(
  ex_list_orig,
  ex_list_no_min_stim,
  ex_tbl_uns_bias,
  chnl_settings,
  bias,
  path_project,
  stage
) {
  .debug("getting loc gate at sample level") # nolint
  force(path_project)

  # get cutpoints for each sample
  cp_uns_loc_obj_list <- purrr::map(
    seq_along(ex_list_no_min_stim),
    function(i) {
      .debug("sample", i) # nolint

      ex_tbl_no_min_stim <- ex_list_no_min_stim[[i]]
      ind <- .get_ind(ex_tbl_no_min_stim)
      .debug("ind", ind) # nolint
      ex_tbl_uns_orig <- ex_list_orig[[length(ex_list_orig)]]
      ex_tbl_stim_orig <- ex_list_orig[[i]]
      chnl <- chnl_settings$chnl_cut %||% .get_cp_uns_loc_get_chnl(ex_tbl_no_min_stim)
      stage_chnl <- file.path(stage, chnl)
      .int_save(
        ind,
        stage_chnl,
        path_project,
        ex_tbl_no_min_stim,
        ex_tbl_uns_orig,
        ex_tbl_stim_orig
      )

      # return early if there are too few cells
      too_few_cells_lgl <- .get_cp_uns_loc_sample_check_cell_number(
        ex_tbl_stim_no_min = ex_tbl_no_min_stim,
        min_cell = chnl_settings$min_cell,
        ex_tbl_uns_bias = ex_tbl_uns_bias
      )
      if (too_few_cells_lgl) {
        obj_out <- .get_cp_uns_loc_sample_too_few(
          stage = stage,
          path_project = path_project,
          ex_tbl_no_min_stim = ex_tbl_no_min_stim,
          ex_tbl_uns_bias = ex_tbl_uns_bias,
          cp_min = chnl_settings$cp_min
        )
        return(obj_out)
      }

      # remove any cytokine-positive cells from unstim using gates from
      # sample for which single-positive gates are required
      ex_tbl_uns_bias <- .get_cp_uns_loc_sample_uns_rm_cyt_pos(
        ex_tbl_uns_orig = ex_tbl_uns_orig,
        chnl_settings = chnl_settings,
        ex_tbl_stim_no_min = ex_tbl_no_min_stim,
        bias = bias,
        ex_tbl_uns_bias = ex_tbl_uns_bias,
        stage = stage
      )
      .int_save(ind, stage_chnl, path_project, ex_tbl_uns_bias)

      .get_cp_uns_loc_ind(
        ex_tbl_uns_bias = ex_tbl_uns_bias,
        ex_tbl_stim_no_min = ex_tbl_no_min_stim,
        chnl_settings = chnl_settings,
        ex_tbl_stim_orig = ex_tbl_stim_orig,
        ex_tbl_uns_orig = ex_tbl_uns_orig,
        bias = bias,
        path_project = path_project,
        stage = stage
      )
    }
  ) |>
    stats::setNames(names(ex_list_no_min_stim))

  .path_project <- path_project
  chnl <- .get_cp_uns_loc_get_chnl(ex_list_no_min_stim[[1]])
  .get_cp_uns_loc_output(
    cp_uns_loc_obj_list = cp_uns_loc_obj_list,
    ind_uns = names(ex_list_orig)[length(ex_list_orig)],
    ind_stim = names(ex_list_no_min_stim),
    stage = stage,
    path_project = .path_project,
    chnl = chnl
  )
}

#' @keywords internal
.get_cp_uns_loc_sample_check_cell_number <- function(
  ex_tbl_stim_no_min,
  min_cell,
  ex_tbl_uns_bias
) {
  nrow(ex_tbl_stim_no_min) < min_cell ||
    nrow(ex_tbl_uns_bias) < min_cell
}

#' @keywords internal
.get_cp_uns_loc_sample_too_few <- function(
  stage,
  path_project,
  ex_tbl_no_min_stim,
  ex_tbl_uns_bias,
  cp_min
) {
  chnl <- .get_cp_uns_loc_get_chnl(ex_tbl_no_min_stim)
  stage_chnl <- file.path(stage, chnl)
  .int_save_nm(
    "too_few_cells_sample_fn",
    NULL,
    .get_ind(ex_tbl_no_min_stim),
    stage_chnl,
    path_project
  ) # nolint
  obj_out <- .get_cp_uns_loc_ind_check_out(
    cp_min = cp_min,
    ex_tbl_stim_no_min = ex_tbl_no_min_stim,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    stage = stage,
    msg = "Too few cells"
  )
  .int_save_nm(
    "cp_ind",
    obj_out$cp,
    .get_ind(ex_tbl_no_min_stim),
    stage_chnl,
    path_project
  )
  obj_out
}

#' @keywords internal
.get_cp_uns_loc_sample_uns_rm_cyt_pos <- function(
  ex_tbl_uns_orig,
  chnl_settings,
  ex_tbl_stim_no_min,
  bias,
  ex_tbl_uns_bias,
  stage
) {
  if (stage == "init") {
    return(ex_tbl_uns_bias)
  }
  .debug("Removing cytokine-positive cells from unstim") # nolint

  # first filter
  gate_tbl_gn_ind <- chnl_settings$gate_tbl |>
    dplyr::filter(
      ind == ex_tbl_stim_no_min$ind[1], # nolint
      .data$gate_name == .env$chnl_settings$gate_name_curr # nolint
    )

  pos_ind_vec_but_single_pos_curr <-
    .get_pos_ind_but_single_pos_for_one_cyt(
      ex = ex_tbl_uns_orig,
      gate_tbl = gate_tbl_gn_ind,
      chnl_single_exc = chnl_settings$chnl_cut,
      chnl = NULL,
      gate_type_cyt_pos = ifelse(chnl_settings$calc_cyt_pos_gates, "cyt", "base"),
      gate_type_single_pos = "base"
    )

  ex_tbl_uns_orig <- ex_tbl_uns_orig[
    !pos_ind_vec_but_single_pos_curr,
    ,
    drop = FALSE
  ]

  # re-apply bias, noise and exclude minimum after removing cytokine-positive cells
  .prepare_ex_list_with_bias_and_noise(
    ex_list = stats::setNames(
      list(ex_tbl_uns_orig),
      attr(ex_tbl_uns_orig, "ind_uns")
    ),
    ind = attr(ex_tbl_uns_orig, "ind_uns"),
    exc_min = chnl_settings$exc_min,
    bias = bias,
    noise_sd = NULL
  )[[1]]
}

#' @keywords internal
.get_cp_uns_loc_p_list_empty <- function() {
  lapply(seq_len(3), function(x) list()) |>
    stats::setNames(c("p_loc_dens", "p_loc_prob", "p_loc_ctb"))
}

#' @keywords internal
.get_cp_uns_loc_get_chnl <- function(ex_tbl) {
  if (is.null(ex_tbl) || !is.data.frame(ex_tbl)) {
    return("unknown_chnl")
  }
  chnl <- attr(ex_tbl, "chnl_cut")
  if (!is.null(chnl) && nzchar(chnl)) {
    return(chnl)
  }
  cut_vec <- .get_cut(ex_tbl)
  cn_vec <- colnames(ex_tbl)
  if (is.null(cn_vec) || length(cn_vec) == 0L) {
    return("unknown_chnl")
  }
  match_ind <- vapply(
    cn_vec,
    function(nm) identical(ex_tbl[[nm]], cut_vec),
    logical(1)
  )
  if (any(match_ind)) {
    return(cn_vec[match_ind][1])
  }
  cn_vec[1]
}

#' @keywords internal
.get_cp_uns_loc_ind <- function(
  ex_tbl_uns_bias,
  ex_tbl_stim_no_min,
  chnl_settings,
  ex_tbl_stim_orig,
  ex_tbl_uns_orig,
  plot = TRUE,
  prob_min = 0.1,
  bias,
  path_project,
  stage
) {
  .debug("getting loc gate for single sample") # nolint
  ind <- .get_ind(ex_tbl_stim_no_min)
  chnl <- chnl_settings$chnl_cut %||% .get_cp_uns_loc_get_chnl(ex_tbl_stim_no_min)
  stage_chnl <- file.path(stage, chnl)
  .debug("ind", ind) # nolint

  # estimate densities for stim and unstim over stim range
  if (.get_cp_uns_loc_check_early(ex_tbl_stim_no_min, chnl_settings$min_cell, chnl_settings$cp_min)) {
    obj_out <- .get_cp_uns_loc_ind_too_few(
      stage = stage,
      path_project = path_project,
      ex_tbl_no_min_stim = ex_tbl_stim_no_min,
      ex_tbl_uns_bias = ex_tbl_uns_bias,
      cp_min = chnl_settings$cp_min
    )
    return(obj_out)
  }

  # stop expr being higher than max_x to prevent really far away values creating modes
  ex_tbl_stim_threshold <- .get_cp_uns_loc_set_max_expr(
    ex_tbl_stim_no_min,
    .get_cp_uns_loc_ind_max_dens_x(ex_tbl_stim_no_min)
  )
  ex_tbl_uns_threshold <- .get_cp_uns_loc_set_max_expr(
    ex_tbl_uns_bias,
    .get_cp_uns_loc_ind_max_dens_x(ex_tbl_stim_no_min)
  )
  .int_save(
    .get_ind(ex_tbl_stim_no_min),
    stage_chnl,
    path_project,
    ex_tbl_stim_threshold,
    ex_tbl_uns_threshold
  )

  # get smoothed probabilities
  data_mod <- .get_cp_uns_loc_get_prob(
    ex_tbl_stim_no_min = ex_tbl_stim_no_min,
    ex_tbl_stim_threshold = ex_tbl_stim_threshold,
    ex_tbl_uns_threshold = ex_tbl_uns_threshold,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    stage = stage,
    bw_min = chnl_settings$bw_min,
    cp_min = chnl_settings$cp_min + bias,
    ex_tbl_uns_orig = ex_tbl_uns_orig,
    path_project = path_project
  )
  .int_save(
    .get_ind(ex_tbl_stim_no_min),
    stage_chnl,
    path_project,
    data_mod
  )

  # get threshold
  .get_cp_uns_loc_get_cp(
    data_mod = data_mod,
    ex_tbl_stim_no_min = ex_tbl_stim_no_min,
    ex_tbl_stim_orig = ex_tbl_stim_orig,
    ex_tbl_uns_orig = ex_tbl_uns_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    bias = bias,
    cp_min = chnl_settings$cp_min,
    stage = stage,
    path_project = path_project
  )
}

.get_cp_uns_loc_ind_too_few <- function(
  stage,
  path_project,
  ex_tbl_no_min_stim,
  ex_tbl_uns_bias,
  cp_min
) {
  chnl <- .get_cp_uns_loc_get_chnl(ex_tbl_no_min_stim)
  stage_chnl <- file.path(stage, chnl)
  .int_save_nm(
    "too_few_cells_ind_fn",
    NULL,
    .get_ind(ex_tbl_no_min_stim),
    stage_chnl,
    path_project
  ) # nolint
  obj_out <- .get_cp_uns_loc_ind_check_out(
    cp_min = cp_min,
    ex_tbl_stim_no_min = ex_tbl_no_min_stim,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    stage = stage,
    msg = "Too few cells"
  )
  .int_save_nm(
    "cp_ind",
    obj_out$cp,
    .get_ind(ex_tbl_no_min_stim), # BUG FIX: Get IND, don't pass whole tbl
    stage_chnl,
    path_project
  )
  return(obj_out)
}

# initial checks
# ---------------------
#' @keywords internal
.get_cp_uns_loc_ind_check_n_cell <- function(ex_tbl_stim_no_min, min_cell) {
  nrow(ex_tbl_stim_no_min) < min_cell
}

#' @keywords internal
.get_cp_uns_loc_ind_max_dens_x <- function(ex_tbl_stim_no_min) {
  max(.get_cut(ex_tbl_stim_no_min)) -
    0.05 * (diff(range(.get_cut(ex_tbl_stim_no_min))))
}

#' @keywords internal
.get_cp_uns_loc_ind_check_max_x <- function(ex_tbl_stim_no_min, cp_min) {
  .get_cp_uns_loc_ind_max_dens_x(ex_tbl_stim_no_min) <= cp_min
}

#' @keywords internal
.get_cp_uns_loc_check_early <- function(ex_tbl_stim_no_min, min_cell, cp_min) {
  .get_cp_uns_loc_ind_check_n_cell(ex_tbl_stim_no_min, min_cell) ||
    .get_cp_uns_loc_ind_check_max_x(ex_tbl_stim_no_min, cp_min)
}

#' @keywords internal
.get_cp_uns_loc_ind_check_out <- function(
  cp_min,
  ex_tbl_stim_no_min,
  ex_tbl_uns_bias,
  stage,
  msg
) {
  .debug(msg) # nolint
  list(
    cp = .get_cp_uns_loc_ind_cp_non_loc(
      cp_min = cp_min,
      ex_tbl_stim_no_min = ex_tbl_stim_no_min,
      ex_tbl_uns_bias = ex_tbl_uns_bias
    ),
    p_list = .get_cp_uns_loc_p_list_empty()
  )
}

#' @keywords internal
.get_cp_uns_loc_ind_cp_non_loc <- function(
  cp_min,
  ex_tbl_stim_no_min,
  ex_tbl_uns_bias
) {
  range_vec_stim <- range(.get_cut(ex_tbl_stim_no_min))
  range_vec_uns <- range(.get_cut(ex_tbl_uns_bias))
  range_len <- max(diff(range_vec_stim), diff(range_vec_uns))
  max(
    cp_min,
    range_vec_stim[[2]] + range_len / 5,
    range_vec_uns[[2]] + range_len / 3
  )
}

#' @keywords internal
.get_cp_uns_loc_set_max_expr <- function(.data, max_x) {
  .data[, attr(.data, "chnl_cut")] <- pmin(.get_cut(.data), max_x)
  .data
}

# get probabilities
#' @keywords internal
.get_cp_uns_loc_get_prob <- function(
  ex_tbl_stim_no_min,
  ex_tbl_stim_threshold,
  ex_tbl_uns_threshold,
  ex_tbl_uns_bias,
  bw_min,
  cp_min,
  ex_tbl_uns_orig,
  stage,
  path_project
) {
  ind <- .get_ind(ex_tbl_stim_no_min)
  chnl <- .get_cp_uns_loc_get_chnl(ex_tbl_stim_no_min)
  stage_chnl <- file.path(stage, chnl)
  
  # get raw densities
  dens_tbl_raw <- .get_cp_uns_loc_get_dens_raw(
    ex_tbl_stim_threshold = ex_tbl_stim_threshold,
    ex_tbl_uns_threshold = ex_tbl_uns_threshold,
    stage = stage,
    path_project = path_project,
    bw_min = bw_min
  )
  .int_save(ind, stage_chnl, path_project, dens_tbl_raw)

  # get probabilities
  prob_tbl_list <- .get_cp_uns_loc_get_prob_tbl(
    dens_tbl_raw = dens_tbl_raw,
    stage = stage,
    cp_min = cp_min,
    ex_vec_stim_threshold = .get_cut(ex_tbl_stim_threshold),
    ex_vec_uns_threshold = .get_cut(ex_tbl_uns_threshold)
  )
  .int_save(ind, stage_chnl, path_project, prob_tbl_list)

  # get .data to smooth over
  data_mod <- .get_cp_uns_loc_get_data_mod(
    ex_tbl_stim_threshold = ex_tbl_stim_threshold,
    ex_tbl_stim_no_min = ex_tbl_stim_no_min,
    ex_tbl_uns_threshold = ex_tbl_uns_threshold,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    prob_tbl_list = prob_tbl_list,
    cp_min = cp_min,
    stage = stage
  )
  .int_save(ind, stage_chnl, path_project, data_mod)

  # smooth
  .get_cp_uns_loc_get_prob_smooth(data_mod, stage, path_project, chnl)
}

# get dens_tbl_raw
# -----------------------
#' @keywords internal
.get_cp_uns_loc_get_dens_raw <- function(
  ex_tbl_stim_threshold,
  ex_tbl_uns_threshold,
  stage,
  path_project,
  bw_min
) {
  .debug("Calculating densities") # nolint

  dens_list <- .get_cp_uns_loc_get_dens_raw_densities(
    ex_tbl_stim_threshold = ex_tbl_stim_threshold,
    ex_tbl_uns_threshold = ex_tbl_uns_threshold,
    stage = stage,
    path_project = path_project,
    bw_min = bw_min
  )

  # put raw densities into table
  .get_cp_uns_loc_get_dens_raw_tabulate(
    stim_x = dens_list$stim$x,
    stim_y = dens_list$stim$y,
    uns_x = dens_list$uns$x,
    uns_y = dens_list$uns$y
  )
}

#' @keywords internal
.get_cp_uns_loc_get_dens_raw_densities <- function(
  ex_tbl_stim_threshold,
  ex_tbl_uns_threshold,
  stage,
  path_project,
  bw_min
) {
  bw <- .get_cp_uns_loc_get_dens_raw_densities_bw(
    ex_tbl_stim_threshold = ex_tbl_stim_threshold,
    ex_tbl_uns_threshold = ex_tbl_uns_threshold,
    bw_min = bw_min
  )
  chnl <- .get_cp_uns_loc_get_chnl(ex_tbl_stim_threshold)
  stage_chnl <- file.path(stage, chnl)
  .int_save_nm(
    "bw_cp_uns_loc",
    bw,
    .get_ind(ex_tbl_stim_threshold),
    stage_chnl,
    path_project
  )
  dens_stim <- .get_cp_uns_loc_get_dens_raw_densities_stim(
    ex_tbl_stim_threshold = ex_tbl_stim_threshold,
    bw = bw
  )
  dens_uns <- .get_cp_uns_loc_get_dens_raw_densities_uns(
    ex_tbl_uns_threshold = ex_tbl_uns_threshold,
    dens_stim = dens_stim,
    bw = bw
  )
  list(stim = dens_stim, uns = dens_uns, bw = bw)
}

#' @keywords internal
.get_cp_uns_loc_get_dens_raw_densities_bw <- function(
  ex_tbl_stim_threshold,
  ex_tbl_uns_threshold,
  bw_min
) {
  bw_stim <- .get_cp_uns_loc_get_dens_raw_densities_bw_init(
    .data = .get_cut(ex_tbl_stim_threshold),
    bw_min = bw_min
  )
  bw_uns <- .get_cp_uns_loc_get_dens_raw_densities_bw_init(
    .data = .get_cut(ex_tbl_uns_threshold),
    bw_min = bw_min
  )
  bw_init <- max(bw_stim, bw_min, bw_uns)
  min(bw_min * 1.5, bw_init)
}

#' @keywords internal
.get_cp_uns_loc_get_dens_raw_densities_bw_init <- function(.data, bw_min) {
  bw_calc <- try(density(.data, bw = "SJ")$bw)
  if (inherits(bw_calc, "try-error")) bw_min else bw_calc
}

#' @keywords internal
.get_cp_uns_loc_get_dens_raw_densities_stim <- function(
  ex_tbl_stim_threshold,
  bw
) {
  dens_obj <- density(.get_cut(ex_tbl_stim_threshold), bw = bw)
  if (is.null(attr(ex_tbl_stim_threshold, "prob_g_min"))) {
    return(dens_obj)
  }
  dens_obj$y <- dens_obj$y * attr(ex_tbl_stim_threshold, "prob_g_min")
  dens_obj
}

#' @keywords internal
.get_cp_uns_loc_get_dens_raw_densities_uns <- function(
  ex_tbl_uns_threshold,
  dens_stim,
  bw
) {
  dens_obj <- density(
    .get_cut(ex_tbl_uns_threshold),
    from = min(dens_stim$x),
    to = max(dens_stim$x),
    bw = bw
  )
  if (is.null(attr(ex_tbl_uns_threshold, "prob_g_min"))) {
    return(dens_obj)
  }
  dens_obj$y <- dens_obj$y * attr(ex_tbl_uns_threshold, "prob_g_min")
  dens_obj
}

#' @keywords internal
.get_cp_uns_loc_get_dens_raw_tabulate <- function(
  stim_x,
  stim_y,
  uns_x,
  uns_y
) {
  dens_tbl_raw_stim <- tibble::tibble(x_stim = stim_x, y_stim = stim_y)
  dens_tbl_raw_wide <- .get_cp_uns_loc_get_dens_raw_tabulate_uns_interp(
    .data = dens_tbl_raw_stim,
    uns_x = uns_x,
    uns_y = uns_y
  )
  .get_cp_uns_loc_get_dens_raw_tabulate_format(dens_tbl_raw_wide)
}

#' @keywords internal
.get_cp_uns_loc_get_dens_raw_tabulate_uns_interp <- function(
  .data,
  uns_x,
  uns_y
) {
  .data |>
    dplyr::mutate(
      y_uns = purrr::map_dbl(
        x_stim, # nolint
        function(marker) .interp(val = marker, x = uns_x, y = uns_y) # nolint
      )
    )
}

#' @keywords internal
.get_cp_uns_loc_get_dens_raw_tabulate_format <- function(.data) {
  .data |>
    tidyr::pivot_longer(y_stim:y_uns, names_to = "stim", values_to = "dens") |> # nolint
    dplyr::mutate(stim = ifelse(stim == "y_stim", "yes", "no")) # nolint
}

# get probabilities
# -------------------
#' @keywords internal
.get_cp_uns_loc_get_prob_tbl <- function(
  dens_tbl_raw,
  stage,
  cp_min,
  ex_vec_stim_threshold,
  ex_vec_uns_threshold
) {
  .debug("Normalising probabilities") # nolint

  prob_tbl <- .get_cp_uns_loc_get_prob_tbl_init(dens_tbl_raw, cp_min)

  prob_tbl_pos <- .get_cp_uns_loc_prob_tbl_filter(
    ex_vec_stim_threshold = ex_vec_stim_threshold,
    ex_vec_uns_threshold = ex_vec_uns_threshold,
    prob_tbl = prob_tbl,
    stage = stage
  )
  list(all = prob_tbl, pos = prob_tbl_pos)
}

#' @keywords internal
.get_cp_uns_loc_get_prob_tbl_init <- function(dens_tbl_raw, cp_min) {
  dens_tbl_raw |>
    tidyr::pivot_wider(
      id_cols = x_stim,
      names_from = stim,
      values_from = dens # nolint
    ) |>
    dplyr::mutate(
      prob_stim = 1 - no / yes, # nolint
      prob_stim = ifelse(yes == 0 & no == 0, 0, prob_stim), # nolint
      prob_stim_norm = pmin(1, prob_stim), # nolint
      prob_stim_norm = pmax(0, prob_stim_norm) # nolint
    ) |>
    dplyr::filter(x_stim > cp_min)
}

#' @keywords internal
.get_cp_uns_loc_prob_tbl_filter <- function(
  ex_vec_stim_threshold,
  ex_vec_uns_threshold,
  prob_tbl,
  stage
) {
  .debug("Filtering before smoothing") # nolint

  density_exc_min_stim <- density(ex_vec_stim_threshold)
  dens_tbl_stim <- tibble::tibble(
    x = density_exc_min_stim$x,
    y = density_exc_min_stim$y
  )
  peak_stim <- dens_tbl_stim |>
    dplyr::filter(y == max(y)) |> # nolint
    dplyr::pull("x") # nolint
  density_exc_min_uns <- density(ex_vec_uns_threshold)
  dens_tbl_uns <- tibble::tibble(
    x = density_exc_min_uns$x,
    y = density_exc_min_uns$y
  )
  peak_uns <- dens_tbl_uns |>
    dplyr::filter(y == max(y)) |> # nolint
    dplyr::pull("x") # nolint
  peak <- max(peak_stim, peak_uns)

  window_width <- 0.15 * diff(quantile(prob_tbl$x_stim, c(0.05, 0.5)))

  prob_tbl <- prob_tbl |>
    dplyr::filter(x_stim > peak + window_width) # nolint

  prob_tbl |>
    dplyr::filter(
      cumsum(prob_stim_norm >= 0.025) > 0
    )

  prob_tbl <- prob_tbl |>
    dplyr::mutate(
      minor_response_ind = prob_stim_norm >= 0.025,
      moderate_response_ind = prob_stim_norm >= 0.075,
      n_remaining = dplyr::n() - seq_len(dplyr::n()) + 1
    )
  prob_tbl |>
    dplyr::filter(
      cumsum(minor_response_ind) > 0 # nolint
    ) |>
    dplyr::mutate(
      prob_larger_count = purrr::map_int(x_stim, function(x) {
        sum(prob_tbl$moderate_response_ind[prob_tbl$x_stim >= x])
      }),
      prob_larger_prop = prob_larger_count / n_remaining # nolint
    ) |>
    dplyr::filter(
      prob_larger_prop > 0.25 # nolint
    ) |>
    dplyr::select(
      -c(
        prob_larger_prop,
        minor_response_ind,
        moderate_response_ind,
        n_remaining,
        prob_larger_count
      )
    )
}

#' @keywords internal
.get_cp_uns_loc_get_min_prob_x <- function(prob_tbl_pos) {
  min(prob_tbl_pos$x_stim)
}

#' @keywords internal
.get_cp_uns_loc_check_response <- function(prob_tbl_pos, ex_tbl_stim_orig) {
  nrow(prob_tbl_pos) == 0 ||
    max(.get_cut(ex_tbl_stim_orig)) <
      .get_cp_uns_loc_get_min_prob_x(prob_tbl_pos)
}

#' @keywords internal
.get_cp_uns_loc_get_data_mod <- function(
  ex_tbl_stim_threshold,
  ex_tbl_stim_no_min,
  ex_tbl_uns_threshold,
  ex_tbl_uns_bias,
  prob_tbl_list,
  cp_min,
  stage
) {
  if (.get_cp_uns_loc_check_response(prob_tbl_list$pos, ex_tbl_stim_no_min)) {
    return(.get_cp_uns_loc_ind_check_out(
      cp_min = cp_min,
      ex_tbl_stim_no_min = ex_tbl_stim_no_min,
      ex_tbl_uns_bias = ex_tbl_uns_bias,
      stage = stage,
      msg = "No responding cells" # nolint
    ))
  }
  margin <- get_cp_uns_loc_get_data_mod_margin(
    ex_tbl_stim_no_min = ex_tbl_stim_no_min,
    ex_tbl_uns_no_min = ex_tbl_uns_threshold
  )

  data_mod <- ex_tbl_stim_threshold
  data_mod <- data_mod[
    .get_cut(data_mod) >=
      (min(.get_cp_uns_loc_get_min_prob_x(prob_tbl_list$pos) - margin)),
  ]
  if (nrow(data_mod) == 0L) {
    return(.get_cp_uns_loc_ind_check_out(
      cp_min = cp_min,
      ex_tbl_stim_no_min = ex_tbl_stim_no_min,
      ex_tbl_uns_bias = ex_tbl_uns_bias,
      stage = stage,
      msg = "No responding cells" # nolint
    ))
  }
  prob_vec <- approx(
    x = prob_tbl_list$pos$x_stim,
    y = prob_tbl_list$pos$prob_stim_norm,
    xout = data_mod[[1]],
    method = "linear",
    f = 0.5,
    rule = 2
  )$y

  data_mod |>
    dplyr::mutate(prob_smooth = prob_vec)
}

get_cp_uns_loc_get_data_mod_margin <- function(
  ex_tbl_stim_no_min,
  ex_tbl_uns_no_min
) {
  abs(max(
    diff(.get_cut(ex_tbl_stim_no_min)),
    diff(.get_cut(ex_tbl_uns_no_min))
  )) *
    0.05
}

# smooth
# ---------------------
#' @keywords internal
.get_cp_uns_loc_get_prob_smooth <- function(
  data_mod,
  stage,
  path_project,
  chnl
) {
  stage_chnl <- file.path(stage, chnl)
  
  if (!.get_cp_uns_loc_get_prob_smooth_check_n_cell(data_mod)) {
    .int_save_nm(
      "not_enough_cells_to_smooth",
      NULL,
      .get_ind(data_mod),
      stage_chnl,
      path_project
    )
    data_mod_out <- .get_cp_uns_loc_get_prob_smooth_check_n_cell_out(data_mod)
    .int_save_nm(
      "prob_smooth_out",
      data_mod_out,
      .get_ind(data_mod),
      stage_chnl,
      path_project
    )
    return(data_mod_out)
  }

  pred_vec <- .get_cp_uns_loc_get_prob_smooth_actual(data_mod, stage)
  data_mod_out <- data_mod |> dplyr::mutate(pred = pred_vec)
  .int_save_nm(
    "prob_smooth_out",
    data_mod_out,
    .get_ind(data_mod),
    stage_chnl,
    path_project
  )
  data_mod_out
}

#' @keywords internal
.get_cp_uns_loc_get_prob_smooth_check_n_cell <- function(data_mod) {
  is.data.frame(data_mod) && nrow(data_mod) > 10
}

#' @keywords internal
.get_cp_uns_loc_get_prob_smooth_check_n_cell_out <- function(data_mod) {
  if (is.data.frame(data_mod)) {
    data_mod |> dplyr::mutate(pred = prob_smooth - 1e-4) # nolint
  } else {
    data_mod
  }
}

#' @keywords internal
.get_cp_uns_loc_get_prob_smooth_actual <- function(data_mod, stage) {
  fit_1 <- .get_cp_uns_loc_get_prob_smooth_actual_first(data_mod, stage)
  .get_cp_uns_loc_get_prob_smooth_actual_first_response(
    fit = fit_1,
    data_mod = data_mod,
    stage = stage
  )
}

#' @keywords internal
.get_cp_uns_loc_get_prob_smooth_actual_first <- function(data_mod, stage) {
  .debug("Smoothing I") # nolint
  try(
    {
      fml <- as.formula(paste0(
        "prob_smooth ~ s(x_stim), bs = 'mpi')"
      ))
      scam::scam(
        fml, # nolint
        family = "binomial",
        .data = data_mod |>
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
      )
    },
    silent = TRUE
  )
}

#' @keywords internal
.get_cp_uns_loc_get_prob_smooth_actual_first_response <- function(
  fit,
  data_mod,
  stage
) {
  if (.get_cp_uns_loc_get_prob_smooth_actual_check(fit, data_mod)) {
    .debug("Smoothed") # nolint
    return(
      .get_cp_uns_loc_get_prob_smooth_actual_response_success(
        fit = fit,
        data_mod = data_mod
      )$pred
    )
  }
  .get_cp_uns_loc_get_prob_smooth_actual_first_response_failure(
    stage = stage,
    data_mod = data_mod
  )
}

#' @keywords internal
.get_cp_uns_loc_get_prob_smooth_actual_check <- function(fit, data_mod) {
  if (inherits(fit, "try-error")) {
    return(FALSE)
  }
  out_list <- .get_cp_uns_loc_get_prob_smooth_actual_response_success(
    fit = fit,
    data_mod = data_mod
  )
  !(all(out_list$pred > 0.99) || out_list$mean_abs_error > 0.3) # nolint
}

#' @keywords internal
.get_cp_uns_loc_get_prob_smooth_actual_response_success <- function(
  fit,
  data_mod
) {
  pred_vec <- predict(fit, type = "response")
  mean_abs_error <- mean(abs(pred_vec - data_mod$prob_smooth))
  list("pred" = pred_vec, "mean_abs_error" = mean_abs_error)
}

#' @keywords internal
.get_cp_uns_loc_get_prob_smooth_actual_first_response_failure <- function(
  stage,
  data_mod
) {
  fit_2 <- .get_cp_uns_loc_get_prob_smooth_actual_second(data_mod, stage)
  if (.get_cp_uns_loc_get_prob_smooth_actual_check(fit_2, data_mod)) {
    .debug("Smoothed") # nolint
    return(
      .get_cp_uns_loc_get_prob_smooth_actual_response_success(
        fit = fit_2,
        data_mod = data_mod
      )$pred
    )
  }
  .get_cp_uns_loc_get_prob_smooth_actual_third(data_mod, stage)
}

#' @keywords internal
.get_cp_uns_loc_get_prob_smooth_actual_second <- function(data_mod, stage) {
  .debug("Smoothing II") # nolint
  try(
    {
      fml <- as.formula(paste0(
        "prob_smooth ~ s(x_stim), bs = 'micv')"
      ))
      scam::scam(
        fml,
        family = "binomial",
        .data = data_mod,
        control = scam::scam.control(
          print.warn = FALSE,
          trace = FALSE,
          devtol.fit = 0.01
        )
      )
    },
    silent = TRUE
  )
}

#' @keywords internal
.get_cp_uns_loc_get_prob_smooth_actual_third <- function(data_mod, stage) {
  .debug("Failed to smooth") # nolint
  data_mod$prob_smooth - 0.0001
}


# get cp
#' @keywords internal
.get_cp_uns_loc_get_cp <- function(
  data_mod,
  ex_tbl_stim_orig,
  ex_tbl_stim_no_min,
  ex_tbl_uns_orig,
  ex_tbl_uns_bias,
  bias,
  cp_min,
  stage,
  path_project
) {
  ind <- .get_ind(ex_tbl_stim_no_min)
  chnl <- .get_cp_uns_loc_get_chnl(ex_tbl_stim_no_min)
  stage_chnl <- file.path(stage, chnl)
  if (!is.data.frame(data_mod)) {
    .int_save_nm("no_data_mod_df", NULL, ind, stage_chnl, path_project)
    .int_save_nm("cp_ind", data_mod, ind, stage_chnl, path_project)
    return(data_mod)
  }

  data_threshold <- .get_cp_uns_loc_get_cp_data_threshold(
    data_mod = data_mod,
    ex_tbl_stim_orig = ex_tbl_stim_orig,
    ex_tbl_stim_no_min = ex_tbl_stim_no_min,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    ex_tbl_uns_orig = ex_tbl_uns_orig,
    bias = bias
  )
  .int_save(ind, stage_chnl, path_project, data_threshold)
  cp_ind <- .get_cp_uns_loc_get_cp_actual(
    data_threshold = data_threshold,
    ex_tbl_stim_no_min = ex_tbl_stim_no_min,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min,
    stage = stage
  )
  .int_save(ind, stage_chnl, path_project, cp_ind)
  .debug("Completed loc gate for single sample") # nolint
  list("cp" = cp_ind, "p_list" = .get_cp_uns_loc_p_list_empty())
}

#' @keywords internal
.get_cp_uns_loc_get_cp_data_threshold <- function(
  data_mod,
  ex_tbl_stim_orig,
  ex_tbl_stim_no_min,
  ex_tbl_uns_bias,
  ex_tbl_uns_orig,
  bias
) {
  data_count <- .get_cp_uns_loc_get_cp_data_threshold_count(data_mod)
  prob_bs_est <- .get_cp_uns_loc_get_cp_data_threshold_prop_bs_est(
    data_count = data_count,
    ex_tbl_stim_orig = ex_tbl_stim_orig
  )
  .get_cp_uns_loc_get_cp_data_threshold_actual(
    data_count = data_count,
    prop_bs_est = prob_bs_est,
    ex_tbl_stim_orig = ex_tbl_stim_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    ex_tbl_uns_orig = ex_tbl_uns_orig,
    bias = bias
  )
}

#' @keywords internal
.get_cp_uns_loc_get_cp_data_threshold_count <- function(data_mod) {
  if (nrow(data_mod) == 1L) {
    min_val <- min(.get_cut(data_mod)) - 1
  } else {
    min_val <- min(.get_cut(data_mod))
  }
  data_mod <- data_mod[.get_cut(data_mod) > min_val, ]
  data_mod <- data_mod[order(.get_cut(data_mod)), ]
  data_mod |>
    dplyr::mutate(n_row = seq_len(dplyr::n())) |>
    dplyr::filter(cumsum(pred > prob_smooth) != n_row) |> # nolint
    dplyr::select(-n_row)
}

#' @keywords internal
.get_cp_uns_loc_get_cp_data_threshold_prop_bs_est <- function(
  data_count,
  ex_tbl_stim_orig
) {
  sum(data_count$pred) / nrow(ex_tbl_stim_orig)
}

#' @keywords internal
.get_cp_uns_loc_get_cp_data_threshold_actual <- function(
  data_count,
  prop_bs_est,
  ex_tbl_stim_orig,
  ex_tbl_uns_bias,
  ex_tbl_uns_orig,
  bias
) {
  data_count <- data_count[order(.get_cut(data_count)), ]
  prop_stim_vec <- purrr::map_dbl(.get_cut(data_count), function(x) {
    sum(.get_cut(ex_tbl_stim_orig) >= x) / nrow(ex_tbl_stim_orig)
  })
  prop_uns_vec <- purrr::map_dbl(.get_cut(data_count), function(x) {
    sum(.get_cut(ex_tbl_uns_orig) >= x) / nrow(ex_tbl_uns_orig)
  })
  data_count |>
    dplyr::mutate(
      prop_stim = prop_stim_vec,
      prop_uns = prop_uns_vec,
    ) |>
    dplyr::mutate(
      prop_bs = prop_stim - prop_uns, # nolint
      prop_bs_diff = prop_bs - prop_bs_est # nolint
    )
}

#' @keywords internal
.get_cp_uns_loc_get_cp_actual <- function(
  data_threshold,
  ex_tbl_stim_no_min,
  ex_tbl_uns_bias,
  cp_min,
  stage
) {
  if (nrow(data_threshold) == 0L) {
    return(.get_cp_uns_loc_ind_check_out(
      cp_min = cp_min,
      ex_tbl_stim_no_min = ex_tbl_stim_no_min,
      ex_tbl_uns_bias = ex_tbl_uns_bias,
      stage = stage,
      msg = "Too few responding cells"
    )[["cp"]])
  }
  data_threshold <- data_threshold |>
    dplyr::filter(abs(prop_bs_diff) == min(abs(prop_bs_diff))) |> # nolint
    dplyr::slice(1) |>
    .get_cut()
  data_threshold
}

#' @keywords internal
.create_combined_identifier <- function(ind_stim) {
  if (is.null(ind_stim) || length(ind_stim) == 0) {
    "empty_batch"
  } else {
    paste(ind_stim, collapse = "_")
  }
}

#' @keywords internal
.get_cp_uns_loc_output <- function(
  cp_uns_loc_obj_list,
  ind_uns,
  ind_stim,
  stage,
  path_project,
  chnl
) {
  stage_chnl <- file.path(stage, chnl)
  cp_vec <- .get_cp_uns_loc_sample_cp_rep(
    stage = stage,
    cp_uns_loc_obj_list = cp_uns_loc_obj_list,
    ind_uns = ind_uns,
    ind_stim = ind_stim,
    path_project = path_project,
    chnl = chnl
  )
  ind_combined <- .create_combined_identifier(ind_stim)
  .int_save(ind_combined, stage_chnl, path_project, cp_vec)
  .debug("done getting loc gate at sample level") # nolint
  list(
    "loc" = cp_vec,
    "p_list" = list()
  )
}

#' @keywords internal
.get_cp_uns_loc_sample_cp_rep <- function(
  stage,
  cp_uns_loc_obj_list,
  ind_uns,
  ind_stim,
  path_project,
  chnl
) {
  .debug("Possibly re-using calculated cutpoints") # nolint
  ind_combined <- .create_combined_identifier(ind_stim)
  stage_chnl <- file.path(stage, chnl)

  cp_vec <- purrr::map_dbl(cp_uns_loc_obj_list, ~ .x[["cp"]])
  .int_save_nm(
    "cp_vec_before_rep",
    cp_vec,
    ind_combined,
    stage_chnl,
    path_project
  )

  if (length(cp_vec) != length(ind_stim)) {
    .int_save_nm(
      "prejoined_cp_used",
      NULL,
      ind_combined,
      stage_chnl,
      path_project
    )
    cp_vec <- stats::setNames(
      rep(cp_vec, length(ind_stim)),
      ind_stim
    )
  } else {
    .int_save_nm(
      "individual_cp_used",
      NULL,
      ind_combined,
      stage_chnl,
      path_project
    )
    cp_vec <- stats::setNames(cp_vec, ind_stim)
  }
  .int_save_nm(
    "cp_vec_after_rep",
    cp_vec,
    ind_combined,
    stage_chnl,
    path_project
  )

  if (!all(purrr::map_lgl(cp_vec, is.na))) {
    .int_save_nm(
      "adding_uns_cp",
      NULL,
      ind_combined,
      stage_chnl,
      path_project
    )
    cp_vec <- c(cp_vec, stats::setNames(mean(cp_vec, na.rm = TRUE), ind_uns))
  } else {
    .int_save_nm(
      "add_na_for_uns_cp",
      NULL,
      ind_combined,
      stage_chnl,
      path_project
    )
    cp_vec <- c(cp_vec, NA)
  }
  .int_save_nm(
    "cp_vec_after_uns_added",
    cp_vec,
    ind_combined,
    stage_chnl,
    path_project
  )

  cp_vec
}
