.get_gate_batch_ind <- function(ex_list, # nolint
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
                                min_cell,
                                params,
                                plot,
                                path_project,
                                debug) {
  .debug(debug, "running .get_gate_batch_ind") # nolint
  if (is.null(params$gate_tbl)) {
    .debug(debug, "params$gate_tbl is NULL") # nolint
    .debug( # nolint
      debug, "gating ind_batch_list", paste0(ind_batch, collapse = "-") # nolint
    )

    # create bare list
    gate_list <- list()

    # Unstim-based cutpoint - FDR
    # ---
    if ("uns" %in% names(gate_combn)) {
      .debug(debug, "getting uns-based cutpoint") # nolint
      cp_uns_list <- .get_cp_uns( # nolint
        ex_list = ex_list,
        ind_gate = ind_batch[ind_in_batch_gate],
        ind_uns = ind_batch[ind_in_batch_uns],
        fdr = fdr, gate_combn = gate_combn[["uns"]],
        data = data, noise_sd = noise_sd, bias_uns = bias_uns,
        cp_min = cp_min, min_cell = min_cell
      )
      gate_list <- gate_list |> append(cp_uns_list)
    }

    # Unstim-based cutpoint - local FDR
    # ---

    if ("loc" %in% names(gate_combn)) {
      cp_uns_loc_list <- .get_cp_uns_loc( # nolint
        ex_list = ex_list,
        ind_gate = ind_batch[ind_in_batch_gate],
        ind_uns = ind_batch[ind_in_batch_uns],
        gate_combn = gate_combn[["loc"]],
        data = data,
        bias_uns = bias_uns,
        noise_sd = NULL,
        bw_min = bw_min,
        cp_min = cp_min, min_cell = min_cell,
        params = params,
        plot = plot,
        path_project = path_project,
        debug = debug
      )
      gate_list <- gate_list |> append(cp_uns_loc_list)
    }

    if ("tg" %in% names(gate_combn)) {
      .debug(debug, "getting straight tg-based cutpoint") # nolint
      gate_list[["tg"]] <- .get_cp_tg( # nolint
        ex_list = ex_list,
        gate_combn = gate_combn[["tg"]],
        cut = cut, tol = tol,
        ind_gate = ind_batch[ind_in_batch_gate],
        exc_min = TRUE,
        min_cell = min_cell,
        cp_min = cp_min,
        bw = bw_min,
        debug = debug
      )
    }

    if (!is.null(params$tol_ctrl)) {
      for (tol in params$tol_ctrl) {
        .debug(debug, "getting tg-based cutpoint as a control") # nolint
        gate_list[[paste0("tg_ctrl_", tol)]] <- .get_cp_tg( # nolint
          ex_list = ex_list,
          gate_combn = "no",
          cut = cut, tol = params$tol_ctrl,
          ind_gate = ind_batch[ind_in_batch_gate],
          exc_min = TRUE,
          min_cell = 0,
          cp_min = 0,
          bw = bw_min
        )
      }
    }

    if (!is.null(params$tol_gate)) {
      .debug(debug, "getting tolerance gate") # nolint
      gate_list[["tg_clust"]] <- .get_cp_tg( # nolint
        ex_list = ex_list,
        gate_combn = "no",
        cut = cut, tol = params$tol_gate,
        ind_gate = ind_batch[ind_in_batch_gate],
        exc_min = TRUE,
        min_cell = 0,
        cp_min = 0,
        bw = bw_min
      )
    }

    gate_tbl <- .get_gate_batch_ind_tbl(gate_list, debug) # nolint
    return(gate_tbl)
  } else {
    .debug(debug, "params$gate_tbl is not NULL") # nolint
    # =================================
    # get pre-adj and -clust gates for each gate type
    # =================================

    gate_tbl <- params$gate_tbl

    gate_tbl <- gate_tbl |>
      dplyr::mutate(
        gate_method = ifelse(
          stringr::str_detect(gate_name, "loc"), "loc", NA # nolint
        ),
        gate_method = ifelse(
          stringr::str_detect(gate_name, "tg"), "tg", gate_method # nolint
        ),
        gate_method = ifelse(
          stringr::str_detect(gate_name, "uns"), "uns", gate_method
        )
      ) |>
      dplyr::mutate(bias = purrr::map_chr(
        gate_name,
        function(x) {
          split_vec <- stringr::str_split(x, "_")[[1]]
          rem_vec <- split_vec |>
            stringr::str_remove("locb|unsb|tg")
          rem_vec[1]
        }
      ) |>
        as.numeric()) |>
      dplyr::mutate(gate_combn = purrr::map_chr(gate_name, function(x) {
        if (stringr::str_detect(x, "_no")) {
          return("no")
        }
        if (stringr::str_detect(x, "_min")) {
          return("min")
        }
        if (stringr::str_detect(x, "_prejoin")) {
          return("prejoin")
        }
      })) |>
      dplyr::mutate(clust = purrr::map_chr(
        gate_name,
        function(x) ifelse(stringr::str_detect(x, "clust"), "clust", "")
      )) |>
      dplyr::mutate(adj = purrr::map_chr(
        gate_name,
        function(x) ifelse(stringr::str_detect(x, "adj"), "adj", "")
      ))

    # create bare list
    gate_list <- list()

    gate_name_vec <- unique(gate_tbl$gate_name)

    # get single-pos gates for each of the gate types already done
    gate_tbl_single <- purrr::map_df(gate_name_vec, function(gate_name_curr) {
      .debug(debug, "getting single-pos gate", gate_name_curr) # nolint
      gate_tbl_gn_marker <- gate_tbl |>
        dplyr::filter(
          gate_name == gate_name_curr, # nolint
          chnl == params$cut # nolint
        )

      gate_name_tbl_row <- gate_tbl_gn_marker[1, , drop = FALSE]

      gate_method <- gate_name_tbl_row$gate_method

      adj_ind <- stringr::str_detect(gate_name_curr, "_adj")

      clust_ind <- stringr::str_detect(gate_name_curr, "_clust")

      # filter to yield cells negative for all cytokine combinations
      # except possible this cytokine single-positive
      ex_list_neg_but_single_pos_curr <- purrr::map(
        seq_along(ex_list),
        function(i) {
          if (i == params$ind_in_batch_uns) {
            return(ex_list[[i]])
          }

          gate_tbl_gn_ind <- gate_tbl |>
            dplyr::filter(
              ind == ex_list[[i]]$ind[1], # nolint
              gate_name == gate_name_curr # nolint
            )

          pos_ind_vec_but_single_pos_curr <-
            .get_pos_ind_but_single_pos_for_one_cyt( # nolint
              ex = ex_list[[i]],
              gate_tbl = gate_tbl_gn_ind,
              chnl_single_exc = params$cut,
              chnl = NULL,
              gate_type_cyt_pos = ifelse(params$calc_cyt_pos_gates,
                "cyt", "base"
              ),
              gate_type_single_pos = "base"
            )
          ex_list[[i]][!pos_ind_vec_but_single_pos_curr, , drop = FALSE]
        }
      ) |>
        stats::setNames(names(ex_list))

      if (gate_method == "tg") {
        # apply method
        gate_list <- .get_cp_tg( # nolint
          ex_list = ex_list_neg_but_single_pos_curr,
          gate_combn = gate_name_tbl_row$gate_combn,
          cut = cut, tol = tol,
          ind_gate = ind_batch[ind_in_batch_gate],
          exc_min = TRUE,
          min_cell = min_cell,
          cp_min = cp_min,
          bw = bw_min
        )

        gate_tbl <- tibble::tibble(
          gate_name = gate_name_curr,
          ind = names(gate_list[[1]]),
          gate = gate_list[[1]],
          gate_use = "gate"
        )
      }

      if (gate_method == "loc") {
        gate_list <- .get_cp_uns_loc( # nolint
          ex_list = ex_list_neg_but_single_pos_curr,
          ind_gate = ind_batch[ind_in_batch_gate],
          ind_uns = ind_batch[ind_in_batch_uns],
          gate_combn = gate_name_tbl_row$gate_combn,
          pop_root = NULL,
          data = data,
          bias_uns = gate_name_tbl_row$bias,
          noise_sd = NULL,
          bw_min = bw_min,
          cp_min = cp_min,
          min_cell = min_cell,
          params = params |> append(list(
            gate_name_curr = gate_name_curr,
            ex_uns = ex_list[[params$ind_in_batch_uns]]
          )),
          plot = plot,
          path_project = path_project
        )
        gate_list <- gate_list[[1]]

        gate_tbl <- tibble::tibble(
          gate_name = gate_name_curr,
          ind = names(gate_list[[1]]),
          gate = gate_list[[1]],
          gate_use = "gate"
        )
      }

      if (gate_method == "uns") {
        gate_list <- .get_cp_uns( # nolint
          ex_list = ex_list_neg_but_single_pos_curr, # nolint
          ind_gate = ind_batch[ind_in_batch_gate],
          ind_uns = ind_batch[ind_in_batch_uns],
          fdr = fdr, gate_combn = gate_name_tbl_row$gate_combn,
          data = data, noise_sd = noise_sd, bias_uns = gate_name_tbl_row$bias,
          cp_min = cp_min, min_cell = min_cell
        )
        gate_list <- gate_list[[1]]

        gate_tbl <- tibble::tibble(
          gate_name = gate_name_curr,
          ind = names(gate_list[[1]]),
          gate = gate_list[[1]]
        )
      }

      if (adj_ind) {
        # apply method
        gate_list <- .get_cp_tg( # nolint
          ex_list = ex_list_neg_but_single_pos_curr,
          gate_combn = "no",
          cut = cut, tol = params$tol_ctrl,
          ind_gate = ind_batch[ind_in_batch_gate],
          exc_min = TRUE,
          min_cell = min_cell,
          cp_min = cp_min,
          bw = bw_min
        )

        gate_tbl_adj <- tibble::tibble(
          gate_name = gate_name_curr,
          ind = names(gate_list[[1]]),
          gate = gate_list[[1]],
          gate_use = "ctrl"
        )
        gate_tbl <- gate_tbl |> dplyr::bind_rows(gate_tbl_adj)
      }

      if (clust_ind) {
        # apply method
        gate_list <- .get_cp_tg( # nolint
          ex_list = ex_list_neg_but_single_pos_curr,
          gate_combn = "no",
          cut = cut, tol = params$tol_gate_single,
          ind_gate = ind_batch[ind_in_batch_gate],
          exc_min = TRUE,
          min_cell = min_cell,
          cp_min = cp_min,
          bw = bw_min
        )

        gate_tbl_clust <- tibble::tibble(
          gate_name = gate_name_curr,
          ind = names(gate_list[[1]]),
          gate = gate_list[[1]],
          gate_use = "tg_clust"
        )
        gate_tbl <- gate_tbl |> dplyr::bind_rows(gate_tbl_clust)
      }

      gate_tbl
    })
  }

  gate_tbl_single |>
    dplyr::mutate(
      gate_type = purrr::map_chr(
        gate_name, # nolint
        function(gn) stringr::str_split(gn, "_")[[1]][1]
      ),
      gate_combn = gate_name |>
        stringr::str_remove("_adj") |>
        stringr::str_remove("_clust") |>
        stringr::str_remove(gate_type) |> # nolint
        stringr::str_remove("_")
    ) |>
    dplyr::mutate(ind = as.numeric(ind)) # nolint
}

.get_gate_batch_ind_tbl <- function(gate_list, debug) {
  purrr::map_df(seq_along(gate_list), function(i) {
    .get_gate_batch_ind_tbl_along_type(gate_list, i, debug)
  })
}

.get_gate_batch_ind_tbl_along_type <- function(gate_list, i, debug) {
  .debug(debug, "gate list index", i) # nolint
  cp_list <- .get_gate_batch_ind_tbl_cp(gate_list[[i]])
  gate_type <- .get_gate_batch_ind_tbl_type(gate_list, i)
  purrr::map_df(seq_along(cp_list), function(j) {
    .get_gate_batch_ind_tbl_along_combn(cp_list, gate_type, j, debug)
  })
}

.get_gate_batch_ind_tbl_along_combn <- function(cp_list,
                                                gate_type,
                                                j,
                                                debug) {
  .debug(debug, "gate list sub-index", j) # nolint
  gate_combn <- .get_gate_batch_ind_tbl_combn(cp_list, j)
  tibble::tibble(
    gate_name = .get_gate_batch_ind_tbl_name(gate_type, gate_combn),
    gate_type = gate_type,
    gate_combn = gate_combn,
    ind = .get_gate_batch_ind_tbl_ind(cp_list, j),
    gate = .get_gate_batch_ind_tbl_gate(cp_list, j),
    gate_use = .get_gate_batch_ind_tbl_use(gate_type)
  )
}

.get_gate_batch_ind_tbl_name <- function(gate_type, gate_combn) {
  paste0(gate_type, "_", gate_combn)
}

.get_gate_batch_ind_tbl_cp <- function(gate_list_elem) {
  if (!"cp" %in% names(gate_list_elem)) {
    return(gate_list_elem)
  }
  gate_list_elem[["cp"]]
}

.get_gate_batch_ind_tbl_type <- function(gate_list, i) {
  names(gate_list)[i]
}
.get_gate_batch_ind_tbl_combn <- function(cp_list, j) {
  names(cp_list)[[j]]
}

.get_gate_batch_ind_tbl_ind <- function(cp_list, j) {
  as.integer(names(cp_list[[j]]))
}
.get_gate_batch_ind_tbl_gate <- function(cp_list, j) {
  cp_list[[j]]
}

.get_gate_batch_ind_tbl_use <- function(gate_type) {
  if (grepl("tg_ctrl_", gate_type)) {
    return("ctrl")
  }
  if (grepl("tg_clust", gate_type)) "tg_clust" else "gate"
}
