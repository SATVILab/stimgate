.gate_cyt_pos <- function(marker_list,
                          pop_gate,
                          data,
                          gate_name = NULL,
                          bw_min,
                          data_name,
                          ind_in_batch_lab_vec,
                          ind_in_batch_gate,
                          debug = FALSE,
                          calc_cyt_pos = TRUE,
                          path_project) {
  .debug(debug, "-------------") # nolint
  .debug(debug, "getting cytokine-positive gates") # nolint
  .debug(debug, "-------------") # nolint

  # prep
  # -------------------------------

  # vector of chanls
  chnl_vec <- .get_cyt_pos_gates_chnl_vec_from_marker_list( # nolint
    marker_list
  )

  # chnl_lab
  chnl_lab_vec <- .get_labs( # nolint
    data = data[[1]],
    cut = chnl_vec
  )

  # params object
  params <- list(
    pop_gate = pop_gate,
    chnl_lab = chnl_lab_vec,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    ind_in_batch_gate = ind_in_batch_gate,
    data_name = data_name
  )

  # get max bw_min for densities
  bw_min <- .gate_cyt_pos_max_bw_min(marker_list)

  # get original gates
  gate_tbl <- .get_cyt_pos_gates_gate_tbl_get( # nolint
    chnl_vec = chnl_vec,
    path_project = path_project,
    params = params,
    gate_name = gate_name,
    debug = debug,
    chnl_lab_vec = chnl_lab_vec
  )

  # keep cytokine-positive gates (gate_cyt)
  # as the original gates, if not actually
  # gating on cytokine-positive-only cells
  if (!calc_cyt_pos) {
    .debug(debug, "Returning original gates as cyt+ gates") # nolint
    return(gate_tbl |> dplyr::mutate(gate_cyt = gate)) # nolint
  }

  # get cyt+ gates for each of the different gate types
  purrr::map_df(unique(gate_tbl$gate_name), function(gn) {
    .get_cyt_pos_gates_gate_name( # nolint
      gate_tbl_gn = gate_tbl |> dplyr::filter(gate_name == gn),
      gn = gn,
      debug = debug,
      data = data,
      chnl_vec = chnl_vec,
      chnl_lab_vec = chnl_lab_vec,
      pop_gate = pop_gate,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec,
      ind_in_batch_gate = ind_in_batch_gate,
      data_name = data_name,
      bw_min = bw_min,
      calc_cyt_pos = calc_cyt_pos,
      path_project = path_project
    )
  })
}

.get_cyt_pos_gates_gate_name <- function(debug,
                                         gn,
                                         gate_tbl_gn,
                                         data,
                                         chnl_vec,
                                         chnl_lab_vec,
                                         pop_gate,
                                         ind_in_batch_lab_vec,
                                         ind_in_batch_gate,
                                         data_name,
                                         bw_min,
                                         calc_cyt_pos,
                                         path_project) {
  .debug(debug, "Getting cyt+ gates for gate_name: ", gn) # nolint
  cp_tbl_cyt <- purrr::map_df(seq_along(data), function(ind) {
    .get_cyt_pos_gates_ind( # nolint
      ind = ind,
      data = data,
      gate_tbl_gn = gate_tbl_gn,
      chnl_vec = chnl_vec,
      chnl_lab_vec = chnl_lab_vec,
      pop_gate = pop_gate,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec,
      ind_in_batch_gate = ind_in_batch_gate,
      data_name = data_name,
      bw_min = bw_min,
      debug = debug,
      calc_cyt_pos = calc_cyt_pos,
      path_project = path_project
    )
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  # join gate_cyt onto gate_tbl
  gate_tbl_gn |>
    dplyr::left_join(
      cp_tbl_cyt |>
        dplyr::select(ind, chnl, marker, gate_cyt), # nolint
      by = c("ind", "chnl", "marker")
    )
}

.get_cyt_pos_gates_ind <- function(ind,
                                   data,
                                   gate_tbl_gn,
                                   chnl_vec,
                                   chnl_lab_vec,
                                   pop_gate,
                                   ind_in_batch_lab_vec,
                                   ind_in_batch_gate,
                                   data_name,
                                   bw_min,
                                   debug,
                                   calc_cyt_pos,
                                   path_project) {
  .debug(debug, "Getting cyt+ gates for ind: ", ind) # nolint

  # prep
  # -------------------------------

  # get which ind in batch
  ind_in_batch <- ind %% length(ind_in_batch_lab_vec)

  # return if ind in batch is the last one, as that is the unstim ind
  if (ind_in_batch == 0) {
    return(NULL)
  }

  # get stim
  stim <- ind_in_batch_lab_vec[[ind_in_batch]]

  # get expression dataframe
  ex <- .get_ex( # nolint
    data = data[[ind]], pop = pop_gate,
    cut = chnl_vec, high = NULL, ind = ind,
    is_uns = FALSE, stim = stim,
    ind_in_batch = ind_in_batch, data_name = data_name
  )

  ind_uns <- (ind %/% 5 * 5) + 5
  ex_uns <- .get_ex( # nolint
    data = data[[ind_uns]], pop = pop_gate, # nolint
    cut = chnl_vec, high = NULL, ind = ind,
    is_uns = FALSE, stim = stim,
    ind_in_batch = ind_in_batch, data_name = data_name
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
        debug = debug
      )
    }
  )

  tibble::tibble(
    ind = ex$ind[1], chnl = chnl_vec, marker = chnl_lab_vec[chnl_vec],
    gate_cyt = cp_vec_cyt_pos
  )
}

.get_cp_pos_gates_chnl <- function(chnl_curr,
                                   ex,
                                   gate_tbl_ind,
                                   bw_min,
                                   debug) {
  .debug(debug, "chnl_curr: ", chnl_curr) # nolint
  if (is.na(gate_tbl_ind$gate[gate_tbl_ind$chnl == chnl_curr])) {
    return(NA)
  }


  # subset only cells pos for at least one other cyt
  # --------------
  non_na_chnl_vec <- gate_tbl_ind |>
    dplyr::filter(!is.na(gate)) |> # nolint
    dplyr::pull("chnl")
  inc_vec <- .get_pos_ind( # nolint
    ex = ex,
    gate_tbl = gate_tbl_ind,
    chnl = setdiff(non_na_chnl_vec, chnl_curr),
    chnl_alt = setdiff(non_na_chnl_vec, chnl_curr),
    gate_type_cyt_pos = "base",
    gate_type_single_pos = "base"
  )

  # subset expression matrix to only cells positive for
  # at least one other marker


  # get original cutpoint
  cp_orig <- gate_tbl_ind |>
    dplyr::filter(chnl == chnl_curr) |> # nolint
    dplyr::pull("gate")

  # =====================
  # Cutpoint - based on cyt- cells
  # =====================

  if (FALSE) {
    cp_neg <- .get_cp_neg(
      ex = ex, inc = inc_vec, chnl = chnl_curr,
      bw_min = bw_min
    )
  }

  # =====================
  # Cutpoint - based on cyt+ cells
  # =====================

  cp_pos <- .get_cp_pos( # nolint
    ex = ex,
    inc = inc_vec,
    chnl = chnl_curr,
    bw_min = bw_min,
    trust_no_or_high_am = FALSE,
    min_cell = 10,
    cp_orig = cp_orig, n_loop = 5
  )

  # =====================
  # Final cutpoint
  # =====================

  cp_cyt_pos <- min(cp_pos, cp_orig, na.rm = TRUE)
  if ("cp_neg" %in% ls()) {
    cp_cyt_pos <- min(cp_cyt_pos, cp_neg, na.rm = TRUE)
  }

  cp_cyt_pos
}

.gate_cyt_pos_max_bw_min <- function(marker, quant = 0.8) {
  marker |>
    purrr::map_dbl(~ .x$bw_min) |>
    quantile(quant)
}
