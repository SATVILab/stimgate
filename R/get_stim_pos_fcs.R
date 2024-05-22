#' @title Write FCS files of marker-positive FCS files
#'
#' @description
#' Uses the gates to write FCS files of marker-positive FCS files.
#'
#' @export
stimgate_fcs_write <- function(data,
                               data_name,
                               pop_gate,
                               chnl = NULL,
                               ind_in_batch_lab_vec,
                               gate_name = NULL,
                               ind_in_batch_gate,
                               dir_base = NULL,
                               ind_in_batch_uns,
                               gate_tbl = NULL,
                               trans_fn = NULL,
                               trans_chnl = NULL,
                               combn_exc = NULL,
                               gate_type_cyt_pos = "cyt",
                               gate_type_single_pos = "single",
                               mult = FALSE,
                               gate_uns = FALSE,
                               gate_uns_method = "min",
                               path_project) {
  # ==============================
  # Preparation
  # ==============================

  # Create params list
  # -----------------------------

  # chnl_lab
  chnl_lab_vec <- .get_labs( # nolint
    data = data[[1]],
    cut = chnl
  )

  # params object
  params <- list(
    pop_gate = pop_gate,
    chnl_lab = chnl_lab_vec,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    ind_in_batch_gate = ind_in_batch_gate,
    data_name = data_name,
    cut = chnl[1]
  )

  # get directory to save to
  # ------------------------------

  # get base directory
  if (is.null(dir_base)) {
    dir_base <- stimgate_dir_base_create( # nolint
      params = params, dir_base_init = path_project
    ) |>
      dirname()
  }


  # get directory to save to
  chnls <- paste0(chnl, collapse = "_")
  if (is.null(combn_exc)) {
    exc <- "exc-none"
  } else {
    if (!is.list(combn_exc)) combn_exc <- list(combn_exc)
    exc <- "exc"
    for (i in seq_along(combn_exc)) {
      combn_exc_curr <- combn_exc[[i]]
      if (i == 1) {
        exc <- paste0(exc, "-", paste0(combn_exc_curr, collapse = "_"))
      } else {
        exc <- paste0(exc, "_&_", paste0(combn_exc_curr, collapse = "_"))
      }
    }
  }

  dir_save <- file.path(dir_base, chnls, exc, gate_name, "fcs")

  if (dir.exists(dir_save)) unlink(dir_save, force = TRUE)
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)

  # ====================================
  # Get gates
  # ====================================

  if (is.null(gate_tbl)) {
    gate_tbl <- purrr::map_df(names(params$chnl_lab), function(chnl_curr) {
      # get base directory
      dir_base <- stimgate_dir_base_create( # nolint
        dir_base_init = path_project,
        params = params |> append(list(cut = chnl_curr))
      )
      # get stats tbl
      gate_tbl <- readRDS(file.path(dir_base, "gate_tbl.rds"))

      if (!is.null(gate_name)) {
        gate_tbl <- gate_tbl |>
          dplyr::filter(gate_name == .env$gate_name) # nolint
      }
    })
  }

  if (is.null(chnl)) chnl <- unique(gate_tbl$chnl)
  if (is.null(gate_name)) gate_name <- unique(gate_tbl$gate_name)

  gate_tbl <- gate_tbl |>
    # dplyr::filter(.data$gate_name == .env$gate_name) |>
    dplyr::filter(chnl %in% .env$chnl) # nolint

  gate_tbl <- gate_tbl |>
    dplyr::mutate(marker = params$chnl_lab[.data$chnl]) |> # nolint
    dplyr::select(
      chnl, marker, gate_name, # nolint
      batch, ind, gate, gate_cyt, gate_single # nolint
    )

  if (gate_uns) {
    calc_uns_gate <- switch(gate_uns_method,
      "min" = min,
      "max" = max,
      "mean" = mean,
      "tmean" = function(x) mean(x, trim = 0.2),
      "med" = median,
      stop("gate_uns_method not recognised")
    )

    gate_tbl_uns <- gate_tbl |>
      dplyr::group_by(chnl, marker, gate_name, batch) |> # nolint
      dplyr::summarise(
        gate = calc_uns_gate(gate), # nolint
        gate_cyt = ifelse("gate_cyt" %in% colnames(gate_tbl),
          calc_uns_gate(gate_cyt), # nolint
          NA
        ),
        gate_single = ifelse("gate_single" %in% colnames(gate_tbl),
          calc_uns_gate(gate_single), # nolint
          NA
        )
      ) |>
      dplyr::ungroup()

    gate_tbl_uns <- gate_tbl_uns |>
      dplyr::mutate(ind = (batch - 1) * ind_in_batch_uns + ind_in_batch_uns) |> # nolint
      dplyr::select(chnl, marker, gate_name, batch, ind, everything()) |> # nolint
      dplyr::arrange(chnl, marker, gate_name, batch, ind) # nolint

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

    gate_tbl <- gate_tbl |>
      dplyr::bind_rows(gate_tbl_uns)
  }

  gate_tbl <- gate_tbl |>
    dplyr::arrange(chnl, marker, gate_name, batch, ind) # nolint

  # ====================================
  # Apply gates
  # ====================================

  ind <- 1
  purrr::walk(seq_along(data), function(ind) {
    ind_in_batch <- ind %% length(ind_in_batch_lab_vec)

    # return if ind in batch is the last one, as that is the unstim ind
    if (ind_in_batch == 0 && !gate_uns) {
      return(NULL)
    }

    ind_in_batch <- ifelse(ind_in_batch == 0, ind_in_batch_uns, ind_in_batch)

    fr <- flowWorkspace::gh_pop_get_data(data[[ind]])
    ex <- flowCore::exprs(fr) |> tibble::as_tibble()

    if (is.na(ex[1, chnl[1]]) && nrow(ex) == 1) {
      return(invisible(TRUE))
    }

    inc_vec <- rep(FALSE, nrow(ex))

    gate_tbl_ind <- gate_tbl |>
      dplyr::filter(.data$ind == .env$ind) # nolint


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

    if (sum(inc_vec) == 0) {
      return(invisible(TRUE))
    }

    # get filtered expression matrix
    ex <- ex[inc_vec, , drop = FALSE]


    if (!is.null(combn_exc)) {
      for (chnl_pos in combn_exc) {
        if (nrow(ex) == 0) next
        exc_vec <- .get_pos_ind_cyt_combn( # nolint
          ex = ex, gate_tbl = gate_tbl_ind,
          chnl_pos = chnl_pos, chnl_neg = setdiff(chnl, chnl_pos),
          chnl_alt = NULL, gate_type_cyt_pos = gate_type_cyt_pos,
          gate_type_single_pos = gate_type_single_pos
        )
        ex <- ex[!exc_vec, , drop = FALSE]
      }
    }

    if (nrow(ex) == 0) {
      return(invisible(TRUE))
    }

    # transform
    if (!is.null(trans_fn)) {
      if (is.null(trans_chnl)) {
        ex <- trans_fn(ex)
      } else {
        for (nm in trans_chnl) {
          ex[, nm] <- trans_fn(ex[, nm])
        }
      }
    }
    flowCore::exprs(fr) <- as.matrix(ex)


    # save
    fn <- flowCore::keyword(fr)[["GUID"]] |> basename()
    fn_out <- file.path(dir_save, fn)
    flowCore::write.FCS(x = fr, filename = fn_out)

    invisible(TRUE)
  })
  dir_save
}
