#' @keywords internal
.gate_cyt_pos <- function(
  chnl_settings,
  ind_batch_list,
  .data,
  gate_name = NULL,
  calc_cyt_pos = TRUE,
  stage,
  path_project
) {
  .debug("-------------") # nolint
  .debug("getting cytokine-positive gates") # nolint
  .debug("-------------") # nolint

  # prep
  # -------------------------------

  # vector of chanls
  chnl_vec <- .get_cyt_pos_gates_chnl_vec_from_chnl_list(
    chnl_settings = chnl_settings
  )

  # chnl_lab
  chnl_lab_vec <- .get_labs(.data = .data[[1]], chnl_cut = chnl_vec)

  # get max bw_min for densities from chnl_settings elements
  bw_min <- .gate_cyt_pos_max_bw_min(chnl_settings)

  # get original gates
  gate_tbl <- .get_cyt_pos_gates_gate_tbl_get(
    chnl_vec = chnl_vec,
    pop = chnl_settings[[1]]$pop_gate,
    path_project = path_project,
    chnl_lab = chnl_lab_vec
  )

  # keep cytokine-positive gates (gate_cyt)
  # as the original gates, if not actually
  # gating on cytokine-positive-only cells
  if (!calc_cyt_pos) {
    .debug("Returning original gates as cyt+ gates") # nolint
    return(gate_tbl |> dplyr::mutate(gate_cyt = gate)) # nolint
  }

  # get cyt+ gates for each of the different gate types
  gn_vec <- unique(gate_tbl$gate_name)
  if (any(grepl("_clust$", gn_vec))) {
    gn_vec <- gn_vec[grepl("_clust$", gn_vec)]
  } else if (any(grepl("_adj$", gn_vec))) {
    gn_vec <- gn_vec[grepl("_adj$", gn_vec)]
  }

  purrr::map_df(gn_vec, function(gn) {
    gate_tbl_gn <- gate_tbl |> dplyr::filter(gate_name == gn)
    force(gate_tbl_gn)
    .get_cyt_pos_gates_gate_name(
      gate_tbl_gn = gate_tbl_gn,
      .data = .data,
      ind_batch_list = ind_batch_list,
      chnl_vec = chnl_vec,
      chnl_lab_vec = chnl_lab_vec,
      pop_gate = chnl_settings[[1]]$pop_gate,
      bw_min = bw_min,
      calc_cyt_pos = calc_cyt_pos,
      stage = stage,
      path_project = path_project
    )
  })
}

#' @keywords internal
.get_cyt_pos_gates_gate_name <- function(
  gate_tbl_gn,
  .data,
  ind_batch_list,
  chnl_vec,
  chnl_lab_vec,
  pop_gate,
  bw_min,
  calc_cyt_pos,
  stage,
  path_project
) {
  .debug(
    "Getting cyt+ gates for gate_name: ",
    gate_tbl_gn$gate_name[[1]]
  ) # nolint
  ind_vec <- unlist(ind_batch_list)

  cp_tbl_cyt <- purrr::map_df(ind_vec, function(ind) {
    ind_uns <- .get_ind_uns(ind, ind_batch_list)
    batch <- .get_batch(ind, ind_batch_list)
    .get_cyt_pos_gates_ind(
      ind = ind,
      .data = .data,
      ind_uns = ind_uns,
      gate_tbl_gn = gate_tbl_gn,
      chnl_vec = chnl_vec,
      chnl_lab_vec = chnl_lab_vec,
      pop_gate = pop_gate,
      bw_min = bw_min,
      calc_cyt_pos = calc_cyt_pos,
      stage = stage,
      path_project = path_project,
      batch = batch
    )
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  # join gate_cyt onto gate_tbl
  gate_tbl_gn |>
    dplyr::left_join(
      cp_tbl_cyt |>
        dplyr::select(batch, ind, chnl, marker, gate_cyt), # nolint
      by = c("batch", "ind", "chnl", "marker")
    )
}

#' @keywords internal
.get_cyt_pos_gates_ind <- function(
  ind,
  .data,
  ind_uns,
  gate_tbl_gn,
  chnl_vec,
  chnl_lab_vec,
  pop_gate,
  bw_min,
  calc_cyt_pos,
  stage,
  batch,
  path_project
) {
  .debug("Getting cyt+ gates for ind: ", ind) # nolint

  # return if ind in batch is the last one, as that is the unstim ind
  if (ind == ind_uns) {
    return(NULL)
  }

  # get expression dataframe
  ex <- .get_ex(
    .data = .data[[ind]],
    pop = pop_gate,
    chnl_cut = chnl_vec,
    ind = ind,
    ind_uns = ind_uns,
    batch = batch,
    path_project = path_project
  )

  ex_uns <- .get_ex(
    .data = .data[[ind_uns]],
    pop = pop_gate, # nolint
    chnl_cut = chnl_vec,
    ind = ind_uns, # BUG FIX: Cache under unstim ind, not stim ind
    ind_uns = ind_uns,
    batch = batch,
    path_project = path_project
  )

  # gates
  # -----------------
  gate_tbl_ind <- gate_tbl_gn |>
    dplyr::filter(.data$ind == .env$ind) # nolint

  # ==============
  # Calculate cyt+ cutpoints
  # ==============
  cp_vec_cyt_pos <- purrr::map_dbl(
    seq_along(chnl_vec),
    function(i) {
      .get_cp_pos_gates_chnl(
        chnl_curr = chnl_vec[[i]],
        ex = ex,
        gate_tbl_ind = gate_tbl_ind,
        bw_min = bw_min,
        ind = ind,
        stage = stage,
        path_project = path_project
      )
    }
  )

  tibble::tibble(
    batch = batch,
    ind = attr(ex, "ind"),
    chnl = chnl_vec,
    marker = chnl_lab_vec[chnl_vec],
    gate_cyt = cp_vec_cyt_pos
  )
}

#' @keywords internal
.get_cp_pos_gates_chnl <- function(
  chnl_curr,
  ex,
  gate_tbl_ind,
  bw_min,
  ind,
  stage,
  path_project
) {
  .debug("chnl_curr: ", chnl_curr) # nolint
  if (is.na(gate_tbl_ind$gate[gate_tbl_ind$chnl == chnl_curr])) {
    return(NA)
  }

  # subset only cells pos for at least one other cyt
  # --------------
  non_na_chnl_vec <- gate_tbl_ind |>
    dplyr::filter(!is.na(gate)) |> # nolint
    dplyr::pull("chnl")
  inc_vec <- .get_pos_ind(
    ex = ex,
    gate_tbl = gate_tbl_ind,
    chnl = setdiff(non_na_chnl_vec, chnl_curr),
    chnl_alt = setdiff(non_na_chnl_vec, chnl_curr),
    gate_type_cyt_pos = "base",
    gate_type_single_pos = "base"
  )
  .int_save_nm(
    paste0(chnl_curr, "_inc_vec"),
    inc_vec,
    ind,
    stage,
    path_project
  )

  # get original cutpoint
  cp_orig <- gate_tbl_ind |>
    dplyr::filter(chnl == chnl_curr) |> # nolint
    dplyr::pull("gate")

  .int_save_nm(
    paste0(chnl_curr, "_cp_orig"),
    cp_orig,
    ind,
    stage,
    path_project
  )

  # =====================
  # Cutpoint - based on cyt+ cells
  # =====================
  cp_pos <- .get_cp_pos(
    ex = ex,
    inc = inc_vec,
    chnl = chnl_curr,
    bw_min = bw_min,
    trust_no_or_high_am = FALSE,
    min_cell = 10,
    cp_orig = cp_orig,
    n_loop = 5
  )
  .int_save_nm(
    paste0(chnl_curr, "_cp_pos"),
    cp_pos,
    ind,
    stage,
    path_project
  )

  # BUG FIX: Actually calculate cp_neg here instead of checking for its existence
  cp_neg <- .get_cp_neg(
    ex = ex,
    inc = inc_vec,
    chnl = chnl_curr,
    bw_min = bw_min,
    min_cell = 10 
  )

  # =====================
  # Final cutpoint
  # =====================
  cp_cyt_pos <- min(cp_pos, cp_orig, na.rm = TRUE)
  if (!is.na(cp_neg)) {
    cp_cyt_pos <- min(cp_cyt_pos, cp_neg, na.rm = TRUE)
  }
  
  .int_save_nm(
    paste0(chnl_curr, "_cp_cyt_pos_final"),
    cp_cyt_pos,
    ind,
    stage,
    path_project
  )

  cp_cyt_pos
}

#' @keywords internal
.gate_cyt_pos_max_bw_min <- function(chnl, quant = 0.8) {
  chnl |>
    purrr::map_dbl(~ .x$bw_min) |>
    quantile(quant)
}
