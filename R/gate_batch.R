#' @title Get cutpoints for a single batch
#'
#' @inheritParams gate # all except those listed below
#' @inheritParams .get_gate_list # pop_gate, cut, high, fdr, tol, gate_combn
.gate_batch <- function(.data,
                        ind_batch,
                        pop_gate,
                        cut,
                        gate_combn,
                        tol,
                        noise_sd,
                        bias_uns,
                        bw_min,
                        cp_min,
                        min_cell,
                        params,
                        debug) {
  # ==========================================
  # .data
  # ==========================================

  # get list of dataframes
  ex_list <- .get_ex_list( # nolint
    .data = .data, ind_batch = ind_batch,
    pop = pop_gate,
    cut = cut, high = high,
    data_name = data_name
  )

  if (is.null(params$gate_tbl)) {
    .gate_batch_all(
      debug, ind_batch,
      ex_list, gate_combn, .data, noise_sd, bias_uns, cp_min,
      min_cell, cut, tol, bw_min, params, path_project
    )
  } else {
    .gate_batch_single(
      debug, ind_batch,
      ex_list, .data, noise_sd,
      bias_uns, cp_min, min_cell, cut, tol, bw_min, params,
      path_project
    )
  }
}


.gate_batch_ex_list <- function(ex_list,
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
