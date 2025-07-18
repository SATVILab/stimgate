setup_project_postmortem <- function(min_cell = TRUE) {
  filter_method <- if (min_cell) "min_cell" else "min_uns"
  path_project <- file.path(tempdir(), "stimgate_gate_run", filter_method)
  suppressWarnings(unlink(path_project, recursive = TRUE))
  dir.create(path_project, showWarnings = FALSE, recursive = TRUE)
  path_project
}

get_fn_tbl_info_postmortem <- function(gs) {
  
  fn_vec_orig <- NULL
  n_cell_vec <- NULL
  for (i in seq_along(gs)) {
    fr <- flowWorkspace::gh_pop_get_data(gs[[i]])
    fn_vec_orig[[i]] <- basename(flowCore::keyword(fr, "FILENAME")[["FILENAME"]])
    n_cell_vec[[i]] <- nrow(flowCore::exprs(fr))
  }
  fn_vec_orig <- fn_vec_orig |> unlist()
  fn_vec <- tolower(fn_vec_orig)
  fn_vec <- gsub("_processed.*|-processed.*", "", fn_vec)
  fn_vec <- gsub(" ", "-", fn_vec)
  pid_vec <- gsub("^pm_", "", fn_vec) |>
    gsub("-.*|_.*", "", x = _)
  pid_vec <- paste0("pid_", pid_vec)
  # remove pid
  fn_vec <- gsub(
    "^pm_\\d+-+|^pm-\\d+-+|^pm-\\d+_+|^pm_\\d+_+", "", fn_vec
  )
  fn_vec <- gsub("-1$|_1$", "", fn_vec)
  fn_vec <- gsub("_", "-", fn_vec)
  location_vec <- gsub("-.*", "", fn_vec)
  stim_vec <- gsub(".*-", "", fn_vec)
  fn_tbl_info <- tibble::tibble(
    fn = fn_vec_orig,
    pid = pid_vec,
    location = location_vec,
    stim = stim_vec,
    n_cell_pop = n_cell_vec |> unlist()
  )
  fn_tbl_info
}

get_batch_list_postmortem <- function(fn_tbl_info,
                                      col_grp,
                                      col_stim,
                                      uns_chr,
                                      col_n_cell,
                                      min_cell_uns,
                                      min_cell_stim,
                                      filter_method) {
  fn_tbl_info[["row_number"]] <- seq_len(nrow(fn_tbl_info))
  grp_vec <- fn_tbl_info[[col_grp[[1]]]]
  for (i in seq_along(col_grp)[-1]) {
    grp_vec <- paste0(grp_vec, "_", fn_tbl_info[[col_grp[[i]]]])
  }
  grp_vec_unique <- unique(grp_vec)
  out_list <- lapply(grp_vec_unique, function(grp) {
    sel_vec_ind <- which(grp_vec == grp)
    sel_vec_ind_uns <- sel_vec_ind[
      which(fn_tbl_info[[col_stim]][sel_vec_ind] == uns_chr)
      ]
    # check only one unstim sample
    if (length(sel_vec_ind_uns) > 1L) {
      stop(paste0("Multiple unstim samples for group ", grp))
    }
    # no unstim, ignore
    if (length(sel_vec_ind_uns) == 0L) {
      return(NULL)
    }
    # check sufficient unstim cells
    if (fn_tbl_info[[col_n_cell]][sel_vec_ind_uns] < min_cell_uns) {
      return(NULL)
    }
    sel_vec_ind_stim <- setdiff(sel_vec_ind, sel_vec_ind_uns)
    # remove stim samples with insufficient cells
    sel_vec_ind_stim <- sel_vec_ind_stim[
      fn_tbl_info[[col_n_cell]][sel_vec_ind_stim] >= min_cell_stim
    ]
    # return NULL if there is no stim sample
    # with sufficient cells
    if (length(sel_vec_ind_stim) == 0L) {
      return(NULL)
    }
    # return with stim indices first
    c(sel_vec_ind_stim, sel_vec_ind_uns)
  })
  # remove stim-and-unstim sets with no stim
  # with sufficient cells and/or no unstim with sufficient
  # cells
  out_list <- out_list[!vapply(out_list, is.null, logical(1))]
  out_list <- out_list |>
    stats::setNames(paste0("batch_", seq_along(out_list)))

  choose_most_min_cell <- filter_method == "min_cell"
  if (choose_most_min_cell) {
    batch_vec_min_cell <- lapply(out_list, function(x) {
      fn_tbl_info[["n_cell_pop"]][x] |>
        min()
    }) |>
      stats::setNames(seq_along(out_list)) |>
      unlist() |>
      sort() |>
      rev() |>
      names() |>
      as.numeric()
    out_list <- out_list[batch_vec_min_cell]
  } else {
    batch_vec_filter_uns <- lapply(out_list, function(x) {
      n_cell_vec <- fn_tbl_info[["n_cell_pop"]][x]
      n_cell_vec[length(n_cell_vec)] >= 1e2
    }) |>
      stats::setNames(seq_along(out_list)) |>
      unlist()
    batch_vec_filter_uns <- batch_vec_filter_uns[batch_vec_filter_uns] |>
      names() |>
      as.numeric()
    out_list <- out_list[batch_vec_filter_uns]
  }
  out_list
}

.get_ind_uns <- function(ind, ind_batch_list) {
  has_ind <- vapply(
    ind_batch_list, function(x) ind %in% x, logical(1)
  )
  if (sum(has_ind) > 1L) {
    # this is an unstim, as it appears
    # in more than one batch
    return(ind)
  }
  has_ind <- which(has_ind)
  ind_batch <- ind_batch_list[has_ind] |>
    unlist()
  ind_batch[[length(ind_batch)]]
}

plot_raw_data_postmortem <- function(filter_method = "min_cell") {
  path_gs <- testthat::test_path(
    "testdata", "postmortem", "gs", "gs_cytof_acs_cd4"
  )
  gs <- flowWorkspace::load_gs(path_gs)
  fn_tbl_info <- get_fn_tbl_info_postmortem(gs)
  # .debugonce(get_batch_list_postmortem)
  batch_list <- get_batch_list_postmortem(
    fn_tbl_info,
    col_grp = c("pid", "location"),
    col_stim = "stim",
    uns_chr = "uns",
    min_cell_uns = 100,
    min_cell_stim = 100,
    col_n_cell = "n_cell_pop",
    filter_method = filter_method
  )
  batch_vec_sort <- lapply(batch_list, function(x) {
    fn_tbl_info[["n_cell_pop"]][x] |>
      min()
  }) |>
    stats::setNames(seq_along(batch_list)) |>
    unlist() |>
    sort()
  batch_vec_sel <- batch_vec_sort[
    batch_vec_sort >= 1000
  ]
  batch_vec_sel <- batch_vec_sel |>
    rev() |>
    names() |>
    as.numeric()
  batch_list <- batch_list[batch_vec_sel] |>
    stats::setNames(paste0("batch_", batch_vec_sel))
  path_project <- file.path(tempdir(), "stimgate_gate_run")
  suppressWarnings(unlink(path_project, recursive = TRUE))
  dir.create(path_project, showWarnings = FALSE, recursive = TRUE)
  marker_vec <- c("Er168Di", "Lu175Di")
  path_dir_plot <- here::here("_tmp", "fig", "raw", "tnf_vs_ifng", filter_method)
  if (dir.exists(path_dir_plot)) {
    unlink(path_dir_plot, recursive = TRUE)
  }
  dir.create(
    path_dir_plot, recursive = TRUE, showWarnings = FALSE
  )
  for (i in seq_along(batch_list)) {
    batch_vec <- batch_list[[i]]
    p_list <- lapply(batch_vec, function(ind) {
      fr <- flowWorkspace::gh_pop_get_data(gs[[ind]])
      stim <- fn_tbl_info[["stim"]][ind]
      
      ex_tbl <- flowCore::exprs(fr) |>
        tibble::as_tibble()
      plot_cyto(
        data = ex_tbl,
        marker = marker_vec,
        exc_min = TRUE,
        limits_expand = list(3),
        limits_equal = TRUE
      ) +
        theme(
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")
        ) +
        labs(x = "IFNg (Er168Di)", y = "TNF (Lu175Di)") +
        ggtitle(stim)
    })
    plot_tbl_list_uv_ifng <- lapply(batch_vec, function(ind) {
      stim <- fn_tbl_info[["stim"]][ind] 
      fr <- flowWorkspace::gh_pop_get_data(gs[[ind]])
      ex_tbl <- flowCore::exprs(fr) |>
        tibble::as_tibble()
      ex_tbl <- ex_tbl |>
        dplyr::mutate(
          stim = stim
        )
      n_row_init <- nrow(ex_tbl)
      ex_tbl_ifng <- ex_tbl |>
        dplyr::filter(Er168Di > min(Er168Di))
      if (nrow(ex_tbl_ifng) == 0L) {
        return(NULL)
      }
      n_row_fin_ifng <- nrow(ex_tbl_ifng)
      prob_g_0_ifng <- n_row_fin_ifng / n_row_init
      dens_obj_ifng <- density(ex_tbl_ifng$Er168Di)
      dens_obj_ifng_adj <- dens_obj_ifng
      dens_obj_ifng_adj$y <- dens_obj_ifng_adj$y * prob_g_0_ifng
      plot_tbl <- tibble::tibble(
        x = dens_obj_ifng$x, y = dens_obj_ifng$y, stim = stim, type = "raw"
      )
      plot_tbl_adj <- tibble::tibble(
        x = dens_obj_ifng$x, y = dens_obj_ifng_adj$y, stim = stim, type = "adj"
      )
      plot_tbl |>
        dplyr::bind_rows(plot_tbl_adj)
    })
    # remove null
    plot_tbl_list_uv_ifng <- plot_tbl_list_uv_ifng[
      vapply(plot_tbl_list_uv_ifng, Negate(is.null), logical(1))
    ]
    plot_tbl_uv_ifng <- Reduce(rbind, plot_tbl_list_uv_ifng)
    p_uv_ifn <- ggplot(
      plot_tbl_uv_ifng,
      aes(x = x, y = y, col = stim, linetype = type)
    ) +
      geom_line() +
      labs(x = "Ifng", y = "Density") +
      cowplot::theme_cowplot() +
      cowplot::background_grid(major = "x") +
      theme(
        plot.background = element_rect(fill = "white")
      )
    plot_tbl_list_uv_tnf <- lapply(batch_vec, function(ind) {
      stim <- fn_tbl_info[["stim"]][ind] 
      fr <- flowWorkspace::gh_pop_get_data(gs[[ind]])
      ex_tbl <- flowCore::exprs(fr) |>
        tibble::as_tibble()
      ex_tbl <- ex_tbl |>
        dplyr::mutate(
          stim = stim
        )
      n_row_init <- nrow(ex_tbl)
      ex_tbl_tnf <- ex_tbl |>
        dplyr::filter(Lu175Di > min(Lu175Di))
      if (nrow(ex_tbl_tnf) == 0L) {
        return(NULL)
      }
      n_row_fin_tnf <- nrow(ex_tbl_tnf)
      prob_g_0_tnf <- n_row_fin_tnf / n_row_init
      dens_obj_tnf <- density(ex_tbl_tnf$Lu175Di)
      dens_obj_tnf_adj <- dens_obj_tnf
      dens_obj_tnf_adj$y <- dens_obj_tnf_adj$y * prob_g_0_tnf
      plot_tbl <- tibble::tibble(
        x = dens_obj_tnf$x, y = dens_obj_tnf$y, stim = stim, type = "raw"
      )
      plot_tbl_adj <- tibble::tibble(
        x = dens_obj_tnf$x, y = dens_obj_tnf_adj$y, stim = stim, type = "adj"
      )
      plot_tbl |>
        dplyr::bind_rows(plot_tbl_adj)
    })
    plot_tbl_list_uv_tnf <- plot_tbl_list_uv_tnf[
      vapply(plot_tbl_list_uv_tnf, Negate(is.null), logical(1))
    ]
    plot_tbl_uv_tnf <- Reduce(rbind, plot_tbl_list_uv_tnf)
    p_uv_tnf <- ggplot(
      plot_tbl_uv_tnf,
      aes(x = x, y = y, col = stim, linetype = type)
    ) +
      geom_line() +
      labs(x = "Tnf", y = "Density") +
      cowplot::theme_cowplot() +
      cowplot::background_grid(major = "x") +
      theme(
        plot.background = element_rect(fill = "white")
      )
    p_list <- p_list |> append(list(p_uv_ifn, p_uv_tnf))
    location <- fn_tbl_info[["location"]][batch_vec[1]]
    path_plot <- file.path(path_dir_plot, paste0(location, "-", names(batch_list)[i] , ".png")) 
    p_grid <- cowplot::plot_grid(plotlist = p_list, ncol = 2) +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")
      ) +
      ggtitle(location)
    ggplot2::ggsave(
      path_plot, width = 15, height = 18, units = "cm"
    )
  }
}

.get_chnl_lab <- function(gs) {
  chnl_lab <- chnl_lab(
    flowWorkspace::gh_pop_get_data(gs[[1]])
  )
  for (i in seq_along(chnl_lab)) {
    chnl_lab[[i]] <- gsub("^.*_", "", chnl_lab[[i]])
  }
  chnl_lab
}

.get_ind_lab <- function(info_tbl) {
  paste0(
    info_tbl$location, "_",
    info_tbl$pid |> gsub("pid_", "", x = _), "_",
    info_tbl$stim
  ) |>
    stats::setNames(seq_len(nrow(info_tbl)))
}
  