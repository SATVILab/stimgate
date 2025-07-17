#' @title Write FCS files of marker-positive FCS files
#'
#' @description
#' Uses the gates to write FCS files of marker-positive FCS files.
#'
#' @param path_project character. Path to project directory.
#' @param .data GatingSet. GatingSet object containing the flow cytometry data.
#' @param ind_batch_list list. List of indices grouped by batch.
#' @param path_dir_save character. Directory path to save the FCS files to.
#' @param chnl character vector. Specific channels to gate on.
#' @param gate_tbl data.frame. Pre-computed gate table, if available.
#' @param trans_fn function. Transformation function to apply.
#' @param trans_chnl character vector. Columns to transform.
#' @param combn_exc list. Combinations of channels to exclude.
#' @param gate_type_cyt_pos character. Gate type to use for cytokine-positive cells.
#' @param gate_type_single_pos character. Gate type to use for single-positive cells.
#' @param mult logical. Whether cells must be multi-positive.
#' @param gate_uns_method character. Method to calculate unstimulated thresholds.
#'
#' @examples
#' \dontrun{
#'   # Write FCS files of cytokine-positive cells
#'   stimgate_fcs_write(
#'     path_project = "/path/to/project",
#'     .data = gs,
#'     ind_batch_list = list(batch1 = 1:10, batch2 = 11:20),
#'     path_dir_save = "/path/to/output"
#'   )
#' }
#' @export
stimgate_fcs_write <- function(path_project, # project directory
                               .data, # gatingset
                               ind_batch_list, # indices by batch
                               path_dir_save, # directory to save to
                               chnl = NULL, # specific channels to gate on
                               gate_tbl = NULL, # whether gate_tbl is pre-available
                               trans_fn = NULL, # transformation to apply
                               trans_chnl = NULL, # columns to transform
                               combn_exc = NULL, # combinations of chnl to exclude
                               gate_type_cyt_pos = "cyt", # gate type to use for cyt-pos cells # nolint
                               gate_type_single_pos = "single", # gate type to use for single-pos cells # nolint
                               mult = FALSE, # whether cells must be multi-positive
                               gate_uns_method = "min") { # how to calculate unstim thresholds # nolint

  # clear and create directory to save to
  if (dir.exists(path_dir_save)) {
    unlink(path_dir_save, force = TRUE, recursive = TRUE)
  }
  dir.create(path_dir_save, recursive = TRUE)

  # get gates
  gate_tbl <- .fcs_write_get_gate_tbl(
    gate_tbl, chnl, .data, ind_batch_list, gate_uns_method,
    gate_type_cyt_pos, gate_type_single_pos
  )

  n_fn <- length(.data)

  purrr::walk(seq_along(.data), function(ind) {
    message(paste0("Writing ", ind, " of ", n_fn, " files"))
    .fcs_write_impl(
      .data, ind, gate_tbl, path_dir_save, chnl, mult,
      gate_type_cyt_pos, gate_type_single_pos, combn_exc,
      trans_fn, trans_chnl
    )
  })
  invisible(path_dir_save)
}




# ================
# Get Gates
# ================

.fcs_write_get_gate_tbl <- function(gate_tbl,
                                    chnl,
                                    .data,
                                    ind_batch_list,
                                    gate_uns_method,
                                    gate_type_cyt_pos,
                                    gate_type_single_pos) {
  gate_tbl |>
    .fcs_write_gate_gate_tbl_gated(chnl) |>
    .fcs_write_get_gate_tbl_add_uns(
      gate_uns_method, gate_type_cyt_pos, gate_type_single_pos, ind_batch_list
    ) |>
    .fcs_write_get_gate_tbl_filter_chnl(chnl) |>
    .fcs_write_get_gate_tbl_add_marker(chnl, .data)
}

.fcs_write_gate_gate_tbl_gated <- function(gate_tbl,
                                           chnl) {
  if (!is.null_gate_tbl) {
    return(gate_tbl)
  }
  purrr::map_df(chnl, function(chnl_curr) {
    readRDS(file.path(path_project, chnl_curr, "gate_tbl.rds"))
  })
}

.fcs_write_get_gate_tbl_add_uns <- function(gate_tbl,
                                            gate_uns_method,
                                            gate_type_cyt_pos,
                                            gate_type_single_pos,
                                            ind_batch_list) {
  gate_tbl_uns <- .fcs_write_get_gate_tbl_add_uns_get_uns(
    gate_uns_method, ind_batch_list
  )
  
  if ("gate_cyt" %in% colnames(gate_tbl)) {
    gate_tbl_uns <- gate_tbl_uns |>
      dplyr::mutate(gate_cyt = pmin(gate, gate_cyt)) # nolint
  } else {
    gate_tbl_uns <- gate_tbl_uns |> dplyr::select(-gate_cyt) # nolint
  }
  if ("gate_single" %in% colnames(gate_tbl)) {
    gate_tbl_uns <- gate_tbl_uns |>
      dplyr::mutate(gate_single = pmax(gate, gate_single)) # nolint
  } else {
    gate_tbl_uns <- gate_tbl_uns |> dplyr::select(-gate_single) # nolint
  }

  gate_tbl |>
    dplyr::bind_rows(gate_tbl_uns)
}

.fcs_write_get_gate_tbl_add_uns_get_uns <- function(gate_uns_method,
                                                    ind_batch_list) {
  calc_uns_gate <- .fcs_write_get_gate_tbl_add_uns_get_uns_calc(
    gate_uns_method
  )

  .fcs_write_get_gate_tbl_add_uns_get_uns_impl(
    gate_tbl, calc_uns_gate, ind_batch_list
  )
}

.fcs_write_get_gate_tbl_add_uns_get_uns_calc <- function(gate_uns_method) {
  switch(gate_uns_method,
    "min" = min,
    "max" = max,
    "mean" = mean,
    "tmean" = function(x) mean(x, trim = 0.2),
    "med" = median,
    stop("gate_uns_method not recognised")
  )
}

.fcs_write_get_gate_tbl_add_uns_get_uns_impl <- function(gate_tbl,
                                                         calc,
                                                         ind_batch_list) {
  gate_tbl |>
    dplyr::group_by(chnl, marker, batch) |> # nolint
    dplyr::summarise(
      ind_stim = paste0(ind |> sort(), collapse = "_"),
      gate = calc(gate), # nolint
      gate_cyt = ifelse("gate_cyt" %in% colnames(gate_tbl),
        calc(gate_cyt), # nolint
        NA
      ),
      gate_single = ifelse("gate_single" %in% colnames(gate_tbl),
        calc(gate_single), # nolint
        NA
      )
    ) |>
    dplyr::ungroup() |>
    .fcs_write_get_gate_tbl_add_uns_get_uns_ind(ind_batch_list)
}

.fcs_write_get_gate_tbl_add_uns_get_uns_ind <- function(gate_tbl,
                                                        ind_batch_list) {
  ind_batch_vec <- lapply(ind_batch_list, function(x) {
    (x[-length(x)]) |> sort() |> paste0(collapse = "_")
  }) |>
    unlist()
  ind_uns_vec <- lapply(ind_batch_list, function(x) x[length(x)]) |>
    unlist()
  ind_vec <- lapply(seq_len(nrow(gate_tbl)), function(x) {
    ind_match <- which(ind_batch_vec == gate_tbl$ind_stim[[x]])
    stopifnot(length(ind_match) == 1L)
    ind_uns_vec[ind_match]
  })
  gate_tbl |>
    dplyr::mutate(ind = ind_vec) |>
    dplyr::select(chnl, marker, batch, ind, everything()) |> # nolint
    dplyr::arrange(chnl, marker, batch, ind)
}



.fcs_write_get_gate_tbl_filter_chnl <- function(gate_tbl, chnl) {
  if (is.null(chnl)) {
    return(gate_tbl)
  }
  gate_tbl |> dplyr::filter(chnl %in% .env$chnl)
}

.fcs_write_get_gate_tbl_add_marker <- function(gate_tbl, chnl, .data) {
  chnl_lab_vec <- .get_labs(.data = .data[[1]], cut = chnl) # nolint
  gate_tbl |>
    dplyr::mutate(marker = chnl_lab_vec[.data$chnl]) |> # nolint
    dplyr::select(
      chnl, marker,
      batch, ind, gate, gate_cyt, gate_single # nolint
    ) |>
    dplyr::arrange(chnl, marker, gate_name, batch, ind)
}

# ===============
# Implementation
# ================

.fcs_write_impl <- function(.data,
                            ind,
                            gate_tbl,
                            path_dir_save,
                            chnl,
                            mult,
                            gate_type_cyt_pos,
                            gate_type_single_pos,
                            combn_exc,
                            trans_fn,
                            trans_chnl) {

  fr <- .fcs_write_impl_load(.data, ind)
  ex <- flowCore::exprs(fr) |> tibble::as_tibble()

  if (is.na(ex[1, chnl[1]]) && nrow(ex) == 1) {
    return(invisible(FALSE))
  }

  gate_tbl_ind <- gate_tbl |>
    dplyr::filter(.data$ind == .env$ind) # nolint

  ex <- .fcs_write_impl_filter(
    ex, gate_tbl_ind, mult, chnl,
    gate_type_cyt_pos, gate_type_single_pos
  )

  if (nrow(ex) == 0) {
    message(paste0("No stimulation-positive cells. No FCS file written."))
    return(invisible(FALSE))
  }

  ex <- .fcs_write_impl_filter_exc(
    ex, combn_exc, gate_tbl_ind, chnl_pos, chnl,
    gate_type_cyt_pos, gate_type_single_pos
  )

  if (nrow(ex) == 0) {
    message(paste0(
      "No cells after excluding particular combinations. No FCS file written."
    ))
    return(invisible(FALSE))
  }

  ex <- .fcs_write_impl_trans(ex, trans_fn, trans_chnl)

  .fcs_write_impl_write(ex, fr, path_dir_save)
  invisible(TRUE)
}

.fcs_write_impl_load <- function(.data, ind) {
  fr <- flowWorkspace::gh_pop_get_data(.data[[ind]])
  if (inherits(fr, "cytoframe")) {
    fr <- flowWorkspace::cytoframe_to_flowFrame(fr)
  }
  fr
}

.fcs_write_impl_filter_inc <- function(ex,
                                       gate_tbl_ind,
                                       mult,
                                       chnl,
                                       gate_type_cyt_pos,
                                       gate_type_single_pos) {
  inc_vec <- rep(FALSE, nrow(ex))

  if (!mult) {
    inc_vec <- .get_pos_ind( # nolint
      ex = ex, gate_tbl = gate_tbl_ind, chnl = chnl, chnl_alt = NULL,
      gate_type_cyt_pos = gate_type_cyt_pos,
      gate_type_single_pos = gate_type_single_pos
    )
  } else {
    inc_vec <- .get_pos_ind_mult( # nolint
      ex = ex, gate_tbl = gate_tbl_ind, chnl = chnl, chnl_alt = NULL,
      gate_type_cyt_pos = gate_type_cyt_pos
    )
  }
  ex[inc_vec, , drop = FALSE]
}

.fcs_write_impl_filter_exc <- function(ex,
                                       combn_exc,
                                       gate_tbl_ind,
                                       chnl_pos,
                                       chnl,
                                       gate_type_cyt_pos,
                                       gate_type_single_pos) {
  if (is.null(combn_exc)) {
    return(ex)
  }
  for (chnl_pos in combn_exc) {
    if (nrow(ex) == 0) break
    exc_vec <- .get_pos_ind_cyt_combn( # nolint
      ex = ex, gate_tbl = gate_tbl_ind,
      chnl_pos = chnl_pos, chnl_neg = setdiff(chnl, chnl_pos),
      chnl_alt = NULL, gate_type_cyt_pos = gate_type_cyt_pos,
      gate_type_single_pos = gate_type_single_pos
    )
    ex <- ex[!exc_vec, , drop = FALSE]
  }
  ex
}

.fcs_write_impl_trans <- function(ex,
                                  trans_fn,
                                  trans_chnl) {
      # transform
  if (is.null(trans_fn)) {
    return(ex)
  }
  if (is.null(trans_chnl)) {
    ex <- trans_fn(ex)
  } else {
    for (nm in trans_chnl) {
      ex[, nm] <- trans_fn(ex[, nm])
    }
  }
  ex
}
.fcs_write_impl_write <- function(ex, 
                                  fr,
                                  path_dir_save) {
  flowCore::exprs(fr) <- as.matrix(ex)
  fn <- flowCore::keyword(fr)[["GUID"]] |> basename()
  fn_out <- file.path(dir_save, fn)
  if (file.exists(fn_out)) {
    invisible(file.remove(fn_out))
  }
  flowCore::write.FCS(x = fr, filename = fn_out)
  message(paste0("Wrote ", fn))
  invisible(TRUE)
}
