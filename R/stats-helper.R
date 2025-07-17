.get_stats_chnl_lab_get <- function(chnl_lab,
                                         .data,
                                         chnl) {
  if (is.null(chnl_lab)) {
    chnl_lab <- .get_labs( # nolint
      .data = .data[[1]],
      chnl_cut = chnl
    )
  }
  chnl_lab
}

.get_stats_params_get <- function(params = NULL,
                                  chnl = NULL,
                                  pop_gate = NULL,
                                  chnl_lab = NULL,
                                  ind_batch_list = NULL,
                                  .data = NULL) {
  if (!is.null(params)) {
    return(params)
  }
  if (is.null(chnl)) {
    stop("chnl must be provided if params is NULL")
  }
  list(
    pop_gate = pop_gate,
    chnl_lab = chnl_lab,
    ind_batch_list = ind_batch_list,
    .data = .data
  )
}

.get_stats_gate_tbl_get <- function(gate_tbl,
                                    chnl_lab,
                                    path_project,
                                    params,
                                    gate_name = NULL,
                                    tol_clust = NULL) {
  if (!is.null(gate_tbl)) {
    return(gate_tbl)
  }
  purrr::map_df(
    names(chnl_lab),
    function(chnl_curr) {
      # get stats tbl
      gate_tbl <- readRDS(file.path(path_project, chnl_curr, "gate_tbl.rds"))

      if (!is.null(gate_name)) {
        gate_tbl <- gate_tbl |>
          dplyr::filter(gate_name == .env$gate_name) # nolint
      }
      if (!is.null(tol_clust)) {
        if (tol_clust) {
          gate_tbl <- gate_tbl |>
            dplyr::filter(grepl("_clust$", gate_name))
        }
      }

      gate_tbl |>
        # dplyr::filter(.data$gate_name == .env$gate_name) |>
        dplyr::mutate(
          chnl = chnl_curr,
          marker = chnl_lab[chnl_curr]
        ) |>
        dplyr::select(
          chnl, marker, gate_name, # nolint
          batch, ind, gate, gate_cyt, gate_single # nolint
        )
    }
  )
}

.get_stats_chnl_get <- function(chnl, gate_tbl) {
  if (!is.null(chnl)) {
    return(chnl)
  }
  unique(gate_tbl$chnl)
}

.get_stats_gate_name_get <- function(gate_name, gate_tbl) {
  if (!is.null(gate_name)) {
    return(gate_name)
  }
  unique(gate_tbl$gate_name)
}

.get_stats_combn_mat_list_get <- function(n_chnl, n_pos) {
  purrr::map(
    seq_len(n_chnl),
    function(n_pos) gtools::combinations(n = n_chnl, r = n_pos)
  ) |>
    stats::setNames(1:n_chnl)
}

.get_stats_cyt_combn_vec_list_get <- function(combn_mat_list,
                                              chnl) {
  purrr::map(
    names(combn_mat_list),
    function(n_pos_nm) {
      combn_mat <- combn_mat_list[[n_pos_nm]]
      purrr::map_chr(seq_len(nrow(combn_mat)), function(i) {
        chnl_pos <- chnl[combn_mat[i, , drop = TRUE]]
        purrr::map_chr(chnl, function(chnl_curr) {
          if (chnl_curr %in% chnl_pos) {
            return(paste0(chnl_curr, "~+~"))
          }
          paste0(chnl_curr, "~-~")
        }) |>
          paste0(collapse = "")
      })
    }
  ) |>
    stats::setNames(names(combn_mat_list))
}

.get_stats_gate_tbl_save <- function(gate_tbl,
                                     path_project,
                                     params,
                                     chnl,
                                     save) {
  if (!save) {
    return(invisible(FALSE))
  }
  if (!"chnl" %in% colnames(gate_tbl)) {
    gate_tbl <- gate_tbl |> dplyr::mutate(
      chnl = params$chnl_cut,
      marker = params$chnl_lab[params$chnl_cut]
    )
  }

  gate_tbl <- gate_tbl |>
    dplyr::select(
      gate_name, chnl, marker, ind, everything() # nolint
    )
  gate_tbl[, "ind"] <- as.character(gate_tbl[["ind"]])

  gate_tbl <- gate_tbl |>
    dplyr::arrange(gate_name, chnl, marker, ind) # nolint
  params[["chnl_cut"]] <- chnl[1]
  dir_save <- file.path(path_project, chnl[1])
  fn_rds <- "gate_tbl.rds"
  fn_csv <- "gate_tbl.csv"
  dir_save_fn_rds <- file.path(dir_save, fn_rds)
  dir_save_fn_csv <- file.path(dir_save, fn_csv)
  if (file.exists(dir_save_fn_rds)) file.remove(dir_save_fn_rds)
  if (file.exists(dir_save_fn_csv)) file.remove(dir_save_fn_csv)
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  write.csv(gate_tbl, dir_save_fn_csv)
  saveRDS(gate_tbl, dir_save_fn_rds)
}



.stats_save <- function(save,
                        stat_tbl,
                        path_project,
                        params,
                        chnl) {
  if (!save) {
    return(invisible(stat_tbl))
  }
  if (!dir.exists(path_project)) {
    dir.create(path_project, recursive = TRUE)
  }
  if ("ind" %in% colnames(stat_tbl)) {
    stat_tbl[, "ind"] <- as.character(stat_tbl[["ind"]])
  }
  if ("batch" %in% colnames(stat_tbl)) {
    stat_tbl[, "batch"] <- as.character(stat_tbl[["batch"]])
  }
  fn_rds <- paste0("gate_stats.rds")
  fn_csv <- paste0("gate_stats.csv")
  path_save_fn_rds <- file.path(path_project, fn_rds)
  path_save_fn_csv <- file.path(path_project, fn_csv)
  if (file.exists(path_save_fn_rds)) file.remove(path_save_fn_rds)
  if (file.exists(path_save_fn_csv)) file.remove(path_save_fn_csv)
  write.csv(stat_tbl, path_save_fn_csv)
  saveRDS(stat_tbl, path_save_fn_rds)
  invisible(path_project)
}
