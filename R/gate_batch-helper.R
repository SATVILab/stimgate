#' @keywords internal
.gate_batch_all <- function(
  ind_batch,
  batch,
  ex_list,
  .data,
  chnl_settings,
  stage,
  path_project
) {
  .debug("chnl_settings$gate_tbl is NULL") # nolint
  .debug(
    "gating ",
    paste0(ind_batch, collapse = "-") # nolint
  )

  # Extract clear scoping variables from chnl_settings
  chnl_cut <- chnl_settings$chnl_cut
  exc_min <- chnl_settings$exc_min
  bw_min <- chnl_settings$bw_min
  tol_ctrl <- chnl_settings$tol_ctrl
  tol_clust <- chnl_settings$tol_clust

  # create bare list
  gate_list <- .get_cp_uns_loc(
    ex_list = ex_list,
    .data = .data,
    chnl_settings = chnl_settings,
    stage = stage,
    path_project = path_project
  )

  if (!is.null(tol_ctrl)) {
    for (tol in tol_ctrl) {
      .debug("getting tg-based cutpoint as a control") # nolint
      gate_list[[paste0("tg_ctrl_", tol)]] <- .get_cp_tg(
        ex_list = ex_list,
        chnl_settings = chnl_settings,
        tg_type = "tol_ctrl",
        stage = stage,
        path_project = path_project
      )
    }
  }

  if (!is.null(tol_clust)) {
    .debug("getting tolerance gate") # nolint
    gate_list[["tg_clust"]] <- .get_cp_tg(
      ex_list = ex_list,
      chnl_settings = chnl_settings,
      tg_type = "tol_clust",
      stage = stage,
      path_project = path_project
    )
  }

  .gate_batch_tbl(gate_list, attr(ex_list[[1]], "batch")) # nolint
}

#' @keywords internal
.gate_batch_tbl <- function(gate_list, batch) {
  purrr::map_df(seq_along(gate_list), function(i) {
    .gate_batch_tbl_along_type(gate_list, batch, i)
  })
}

#' @keywords internal
.gate_batch_tbl_along_type <- function(gate_list, batch, i) {
  .debug("gate list index", i) # nolint
  cp_list <- .gate_batch_tbl_cp(gate_list[[i]])
  gate_type <- .gate_batch_tbl_type(gate_list, i)
  purrr::map_df(seq_along(cp_list), function(j) {
    .gate_batch_tbl_along_combn(cp_list, gate_type, batch, j)
  })
}

#' @keywords internal
.gate_batch_tbl_along_combn <- function(cp_list, gate_type, batch, j) {
  .debug("gate list sub-index", j) # nolint
  gate_combn <- .gate_batch_tbl_combn(cp_list, j)
  tibble::tibble(
    gate_name = .gate_batch_tbl_name(gate_type, gate_combn),
    gate_type = gate_type,
    gate_combn = gate_combn,
    batch = batch,
    ind = .gate_batch_tbl_ind(cp_list, j),
    gate = .gate_batch_tbl_gate(cp_list, j),
    gate_use = .gate_batch_tbl_use(.env$gate_type)
  )
}

#' @keywords internal
.gate_batch_tbl_name <- function(gate_type, gate_combn) {
  paste0(gate_type, "_", gate_combn)
}

#' @keywords internal
.gate_batch_tbl_cp <- function(gate_list_elem) {
  if (!"cp" %in% names(gate_list_elem)) {
    return(gate_list_elem)
  }
  gate_list_elem[["cp"]]
}

#' @keywords internal
.gate_batch_tbl_type <- function(gate_list, i) {
  names(gate_list)[i]
}
#' @keywords internal
.gate_batch_tbl_combn <- function(cp_list, j) {
  names(cp_list)[[j]]
}

#' @keywords internal
.gate_batch_tbl_ind <- function(cp_list, j) {
  as.character(names(cp_list[[j]]))
}
#' @keywords internal
.gate_batch_tbl_gate <- function(cp_list, j) {
  cp_list[[j]]
}

#' @keywords internal
.gate_batch_tbl_use <- function(gate_type) {
  if (grepl("tg_ctrl_", gate_type)) {
    return("ctrl")
  }
  if (grepl("tg_clust", gate_type)) "tg_clust" else "gate"
}

#' @keywords internal
.gate_batch_single <- function(
  ind_batch,
  batch,
  ex_list,
  .data,
  chnl_cut,
  chnl_settings,
  stage,
  path_project
) {
  .debug("chnl_settings$gate_tbl is not NULL") # nolint
  .debug(paste0("Gating batch ", batch))

  # =================================
  # get pre-adj and -clust gates for each gate type
  # =================================
  gate_tbl <- .gate_batch_single_tbl_format(chnl_settings$gate_tbl)
  gate_name_vec <- unique(gate_tbl$gate_name)
  
  # Scoped extraction to replace old signature variables
  ind_in_batch_uns <- chnl_settings$ind_in_batch_uns
  calc_cyt_pos_gates <- chnl_settings$calc_cyt_pos_gates
  tol_ctrl <- chnl_settings$tol_ctrl
  tol_clust_single <- chnl_settings$tol_clust_single

  # get single-pos gates for each of the gate types already done
  purrr::map_df(gate_name_vec, function(gate_name_curr) {
    .debug("getting single-pos gate", gate_name_curr) # nolint
    gate_tbl_gn_marker <- gate_tbl |>
      dplyr::filter(gate_name == gate_name_curr, chnl == chnl_cut) # nolint

    gate_name_tbl_row <- gate_tbl_gn_marker[1, , drop = FALSE]
    gate_method <- gate_name_tbl_row$gate_method
    adj_ind <- stringr::str_detect(gate_name_curr, "_adj")
    clust_ind <- stringr::str_detect(gate_name_curr, "_clust")

    # filter to yield cells negative for all cytokine combinations
    # except possibly this cytokine single-positive
    ex_list_neg_but_single_pos_curr <- purrr::map(
      seq_along(ex_list),
      function(i) {
        if (i == ind_in_batch_uns) {
          return(ex_list[[i]])
        }

        gate_tbl_gn_ind <- gate_tbl |>
          dplyr::filter(
            ind == attr(ex_list[[i]], "ind"), # nolint
            gate_name == gate_name_curr # nolint
          )

        pos_ind_vec_but_single_pos_curr <-
          .get_pos_ind_but_single_pos_for_one_cyt(
            ex = ex_list[[i]],
            gate_tbl = gate_tbl_gn_ind,
            chnl_single_exc = chnl_cut,
            chnl = NULL,
            gate_type_cyt_pos = ifelse(
              calc_cyt_pos_gates,
              "cyt",
              "base"
            ),
            gate_type_single_pos = "base"
          )
        ex_list[[i]][!pos_ind_vec_but_single_pos_curr, , drop = FALSE]
      }
    ) |>
      stats::setNames(names(ex_list))

    # Temporarily overlay single execution scope overrides into local settings clone
    chnl_settings_exec <- chnl_settings
    chnl_settings_exec$gate_combn <- gate_name_tbl_row$gate_combn

    gate_list <- switch(
      gate_method,
      "tg" = .get_cp_tg(
        ex_list = ex_list_neg_but_single_pos_curr,
        chnl_settings = chnl_settings_exec,
        tg_type = "tg",
        stage = stage,
        path_project = path_project
      ),
      "loc" = {
        chnl_settings_exec$gate_name_curr <- gate_name_curr
        chnl_settings_exec$ex_uns <- ex_list[[length(ex_list)]]
        .get_cp_uns_loc(
          ex_list = ex_list_neg_but_single_pos_curr,
          .data = .data,
          chnl_settings = chnl_settings_exec,
          stage = stage,
          path_project = path_project
        )[[1]][[1]]
      }
    )

    if (names(gate_list)[[1]] == "cp") {
      gate_list <- gate_list[["cp"]]
    }

    gate_tbl_out <- switch(
      gate_method,
      "tg" = ,
      "loc" = ,
      "uns" = tibble::tibble(
        gate_name = gate_name_curr,
        ind = names(gate_list[[1]]),
        gate = gate_list[[1]],
        gate_use = "gate"
      )
    )

    if (adj_ind) {
      chnl_settings_adj <- chnl_settings
      chnl_settings_adj$gate_combn <- "no"
      chnl_settings_adj$tol <- tol_ctrl

      gate_list <- .get_cp_tg(
        ex_list = ex_list_neg_but_single_pos_curr,
        chnl_settings = chnl_settings_adj,
        tg_type = "adj",
        stage = stage,
        path_project = path_project
      )

      if (names(gate_list)[[1]] == "cp") {
        gate_list <- gate_list[["cp"]]
      }

      gate_tbl_adj <- tibble::tibble(
        gate_name = gate_name_curr,
        ind = names(gate_list[[1]]),
        gate = gate_list[[1]],
        gate_use = "ctrl"
      )
      gate_tbl_out <- gate_tbl_out |> dplyr::bind_rows(gate_tbl_adj)
    }

    if (clust_ind) {
      chnl_settings_clust <- chnl_settings
      chnl_settings_clust$gate_combn <- "no"
      chnl_settings_clust$tol <- tol_clust_single

      gate_list <- .get_cp_tg(
        ex_list = ex_list_neg_but_single_pos_curr,
        chnl_settings = chnl_settings_clust,
        tg_type = "clust",
        stage = stage,
        path_project = path_project
      )

      if (names(gate_list)[[1]] == "cp") {
        gate_list <- gate_list[["cp"]]
      }

      gate_tbl_clust <- tibble::tibble(
        gate_name = gate_name_curr,
        ind = names(gate_list[[1]]),
        gate = gate_list[[1]],
        gate_use = "tg_clust"
      )
      gate_tbl_out_final <- gate_tbl_out |>
        dplyr::bind_rows(gate_tbl_clust)
    } else {
      gate_tbl_out_final <- gate_tbl_out
    }

    gate_tbl_out_final
  }) |>
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

#' @keywords internal
.gate_batch_single_tbl_format <- function(gate_tbl) {
  gate_tbl |>
    dplyr::mutate(
      gate_method = .gate_batch_single_tbl_format_method(gate_name) # nolint
    ) |>
    dplyr::mutate(
      gate_combn = .gate_batch_single_tbl_format_combn(gate_name) # nolint
    ) |>
    dplyr::mutate(
      clust = purrr::map_chr(
        gate_name,
        function(x) ifelse(stringr::str_detect(x, "clust"), "clust", "")
      )
    ) |>
    dplyr::mutate(
      adj = purrr::map_chr(
        gate_name,
        function(x) ifelse(stringr::str_detect(x, "adj"), "adj", "")
      )
    )
}

#' @keywords internal
.gate_batch_single_tbl_format_method <- function(gate_name) {
  gate_method <- ifelse(
    stringr::str_detect(gate_name, "loc"),
    "loc",
    NA # nolint
  )
  gate_method <- ifelse(
    stringr::str_detect(gate_name, "tg"),
    "tg",
    gate_method # nolint
  )
  ifelse(
    stringr::str_detect(gate_name, "uns"),
    "uns",
    gate_method
  )
}

#' @keywords internal
.gate_batch_single_tbl_format_combn <- function(gate_name) {
  purrr::map_chr(gate_name, function(x) {
    if (stringr::str_detect(x, "_no")) {
      return("no")
    }
    if (stringr::str_detect(x, "_min")) {
      return("min")
    }
    if (stringr::str_detect(x, "_prejoin")) {
      return("prejoin")
    }
  })
}