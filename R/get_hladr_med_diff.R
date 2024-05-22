get_hladr_med_diff <- function(data,
                               data_name,
                               pop_gate,
                               chnl,
                               mult,
                               gate_tbl = NULL,
                               ind_in_batch_lab_vec,
                               gate_name,
                               ind_in_batch_gate,
                               ind_in_batch_uns,
                               gate_type_cyt_pos = "cyt",
                               gate_type_single_pos = "single",
                               path_project) {
  # ==============================
  # Preparation
  # ==============================

  # Create params list
  # -----------------------------

  # dataset_name
  # data_name <- deparse(substitute(data))

  # chnl_lab
  chnl_lab_vec <- .get_labs(
    data = data[[1]],
    cut = chnl
  )

  # params object
  params <- list(
    pop_gate = pop_gate,
    chnl_lab = chnl_lab_vec,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    ind_in_batch_gate = ind_in_batch_gate,
    data_name = data_name
  )

  # ====================================
  # Get gates
  # ====================================

  if (is.null(gate_tbl)) {
    gate_tbl <- purrr::map_df(names(params$chnl_lab), function(chnl_curr) {
      params[["cut"]] <- chnl_curr
      # get base directory
      dir_base <- stimgate_dir_base_create(
        dir_base_init = path_project,
        params = params
      )
      # get stats tbl
      gate_tbl <- readRDS(file.path(dir_base, "gate_tbl.rds"))

      if (!is.null(gate_name)) gate_tbl <- gate_tbl |> dplyr::filter(gate_name == .env$gate_name)

      gate_tbl |>
        # dplyr::filter(.data$gate_name == .env$gate_name) |>
        dplyr::mutate(chnl = chnl_curr, marker = params$chnl_lab[chnl_curr]) |>
        dplyr::select(chnl, marker, gate_name, batch, ind, gate, gate_cyt, gate_single)
    })
  }

  # ====================================
  # Apply gates
  # ====================================

  ind <- 2
  stat_tbl <- purrr::map(seq_along(data), function(ind) {
    # get which ind in batch

    ind_in_batch <- ind %% length(ind_in_batch_lab_vec)

    # return if ind in batch is the last one, as that is the unstim ind
    if (ind_in_batch == 0) {
      return(NULL)
    }

    # get stim
    stim <- ind_in_batch_lab_vec[[ind_in_batch]]

    # get expression dataframe
    ex <- .get_ex(
      data = data[[ind]], pop = pop_gate,
      cut = chnl, high = NULL, ind = ind,
      is_uns = FALSE, stim = stim,
      ind_in_batch = ind_in_batch, data_name = data_name
    )


    n_cell <- nrow(ex)

    gate_tbl_ind <- gate_tbl |> dplyr::filter(.data$ind == .env$ind)

    if (!mult) {
      inc_vec <- .get_pos_ind(
        ex = ex, gate_tbl = gate_tbl_ind, chnl = chnl, chnl_alt = NULL,
        gate_type_cyt_pos = gate_type_cyt_pos,
        gate_type_single_pos = gate_type_single_pos
      )
    } else {
      inc_vec <- .get_pos_ind_mult(
        ex = ex, gate_tbl = gate_tbl_ind, chnl = chnl, chnl_alt = NULL,
        gate_type_cyt_pos = gate_type_cyt_pos
      )
    }

    if (sum(inc_vec) %in% c(0, n_cell)) {
      hladr_med_pos <- NA
      hladr_med_neg <- NA
      hladr_med_diff <- 0
    } else {
      ex_pos <- ex[inc_vec, ]
      ex_neg <- ex[!inc_vec, ]
      hladr_med_pos <- median(ex_pos$Eu151Di)
      hladr_med_neg <- median(ex_neg$Eu151Di)
      hladr_med_diff <- hladr_med_pos - hladr_med_neg
    }

    ex |>
      dplyr::select(batch:stim, ind, is_uns, ind_in_batch) |>
      dplyr::slice(1) |>
      dplyr::mutate(
        hladr_med_pos = hladr_med_pos,
        hladr_med_neg = hladr_med_neg,
        hladr_med_diff = hladr_med_diff
      )
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()
}
