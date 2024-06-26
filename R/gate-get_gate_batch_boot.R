#' @title Get cutpoints for a single batch
#'
#' @inheritParams gate # all except those listed below
#' @inheritParams .get_gate_list # pop_gate, cut, high, fdr, tol, gate_combn
.get_gate_batch_boot <- function(data,
                                 ind_batch,
                                 ind_in_batch_gate,
                                 ind_in_batch_uns,
                                 ind_in_batch_lab_vec,
                                 pop_gate,
                                 cut,
                                 high,
                                 gate_combn,
                                 pop_man,
                                 pop_man_match_exact,
                                 tol,
                                 data_name,
                                 fdr,
                                 noise_sd,
                                 bias_uns,
                                 bw_min,
                                 cp_min,
                                 boot_sd,
                                 boot_n,
                                 min_cell,
                                 params,
                                 plot,
                                 debug) {
  # ==========================================
  # Get manual gates
  # ==========================================

  # get list of dataframes for manual gating
  gate_tbl_man <- .get_gate_batch_boot_man(
    pop_man = pop_man, data = data, ind_batch = ind_batch,
    ind_in_batch_gate = ind_in_batch_gate,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    pop_gate = pop_gate, cut = cut, high = high,
    gate_combn = gate_combn, data_name = data_name,
    pop_man_match_exact = pop_man_match_exact,
    boot_n = boot_n
  )

  # ==========================================
  # Data
  # ==========================================

  # get list of dataframes
  ex_list <- .get_ex_list( # nolint
    data = data, ind_batch = ind_batch,
    ind_in_batch_gate, ind_in_batch_uns,
    ind_in_batch_lab_vec,
    pop = pop_gate,
    cut = cut, high = high,
    data_name = data_name
  )

  # ==========================================
  # Get automated cutpoints
  # ==========================================

  # set number of bootstrap samples to 1


  par_list_boot <- .get_get_batch_boot_par(boot_n, boot, boot_sd)

  gate_tbl_auto <- purrr::map_df(seq_len(par_list_boot$boot_n), function(i) {
    ex_list_boot <- .get_gate_batch_boot_ex_list(
      ex_list = ex_list, boot = par_list_boot[["boot"]],
      boot_n = par_list_boot$boot_n, boot_sd = par_list_boot$boot_sd, i = i
    )

    .get_gate_batch_ind( # nolint
      ex_list = ex_list_boot,
      ind_batch = ind_batch,
      ind_in_batch_gate = ind_in_batch_gate,
      ind_in_batch_uns = ind_in_batch_uns,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec,
      pop_gate = pop_gate,
      cut = cut,
      high = high,
      gate_combn = gate_combn,
      pop_man = pop_man,
      pop_man_match_exact = pop_man_match_exact,
      tol = tol,
      data_name = data_name,
      fdr = fdr,
      noise_sd = noise_sd,
      bias_uns = bias_uns,
      bw_min = bw_min,
      cp_min = cp_min,
      min_cell = min_cell,
      params = params,
      plot = plot,
      debug = debug
    ) |>
      dplyr::mutate(boot = par_list_boot[["boot"]], boot_ind = i) |>
      dplyr::select(boot, boot_ind, everything()) # nolint
  })

  if (is.null(pop_man)) {
    return(gate_tbl_auto)
  }

  gate_tbl_man |>
    dplyr::bind_rows(gate_tbl_auto)
}

.get_gate_batch_boot_man <- function(pop_man,
                                     data,
                                     ind_batch,
                                     ind_in_batch_gate,
                                     ind_in_batch_uns,
                                     ind_in_batch_lab_vec,
                                     pop_gate,
                                     cut,
                                     high,
                                     gate_combn,
                                     data_name,
                                     pop_man_match_exact,
                                     boot_n) {
  if (is.null(pop_man)) {
    return(NULL)
  }

  ex_list_man <- .get_ex_list( # nolint
    data = data, ind_batch = ind_batch,
    ind_in_batch_gate, ind_in_batch_uns,
    ind_in_batch_lab_vec,
    pop = pop_gate,
    cut = cut, high = high,
    data_name = data_name
  )

  # calculate manual cutpoint for each sample
  # ---
  gate_list_man <- .get_cp_man( # nolint
    data = data, ind = ind_batch[ind_in_batch_gate],
    ex_list = ex_list_man,
    pop_man = pop_man,
    pop_man_match_exact = pop_man_match_exact,
    gate_combn = gate_combn[["man"]]
  )

  tibble::tibble(
    boot = !is.null(boot_n), boot_ind = 1,
    gate_name = "man_no",
    gate_type = "man", gate_combn = "no",
    ind = names(gate_list_man[[1]]) |> as.integer(),
    gate = gate_list_man[[1]],
    gate_use = "gate"
  )
}

.get_gate_batch_boot_ex_list <- function(ex_list,
                                         boot,
                                         boot_n,
                                         boot_sd,
                                         i) {
  if (!boot) {
    return(ex_list)
  }

  if (boot) {
    if (i %% 20 == 0) print(paste0("boot # ", i, " of ", boot_n))
    # get bootstrap sample
    ex_list_boot <- ex_list |> purrr::map(
      function(ex) ex |> dplyr::sample_frac(size = 1, replace = TRUE)
    )
    # add noise to marker to gate on
    if (!boot_sd == 0) {
      ex_list_boot <- ex_list_boot |>
        purrr::map(function(ex) {
          ex |> dplyr::mutate(cut = cut + rnorm(nrow(ex), sd = boot_sd))
        })
    }
  }
  ex_list_boot
}

  .get_get_batch_boot_par <- function(boot_n, boot, boot_sd) {
    if (is.null(boot_n)) {
      boot <- FALSE
      boot_n <- 1
    } else {
      boot <- TRUE
    }
    if (is.null(boot_sd)) boot_sd <- 0
    list(boot_n = boot_n, boot = boot, boot_sd = boot_sd)
  }