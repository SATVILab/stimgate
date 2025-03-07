.get_cp_cluster_plot_thresholds <- function(cp_tbl,
                                            params,
                                            marker_list,
                                            chnl,
                                            gs,
                                            ind_batch_list,
                                            debug) {
  # diff by cp
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = cp_orig_quant_min - cp_join_lse_orig_mean) |> # nolint
    dplyr::filter(diff == max(diff)) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)

  # diff by prop
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = prop_bs_cp_diff - prop_bs_cp_diff) |> # nolint
    dplyr::filter(diff == max(abs(diff))) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)

  cp_vec <- stats::setNames(sel_tbl$cp_join_lse_orig, sel_tbl$grp)

  # tg
  # -------------------------

  # diff by cp
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = cp_orig_quant_min - cp_join_tg_orig_mean) |> # nolint
    dplyr::filter(diff == max(diff)) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)

  # diff by prop
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = prop_bs_cp_diff - prop_bs_cp_diff) |> # nolint
    dplyr::filter(diff == max(abs(diff))) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)

  # diff by cp
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = cp_orig_quant_min - cp_join_lse_orig_mean_tg) |> # nolint
    dplyr::filter(diff == max(diff)) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)

  # diff by prop
  sel_tbl <- cp_tbl |>
    dplyr::group_by(grp) |> # nolint
    dplyr::mutate(diff = prop_bs_cp_diff - prop_bs_cp_diff) |> # nolint
    dplyr::filter(diff == max(abs(diff))) |>
    dplyr::arrange(cp_join) |> # nolint
    dplyr::slice(1)


  cp_vec <- stats::setNames(sel_tbl$cp_join_lse_orig_mean_tg, sel_tbl$grp)

  plot_list <- purrr::map(sel_tbl$ind, function(i) {
    ind_batch <- which(
      purrr::map_lgl(
        params$ind_batch_list,
        function(ind_batch) i %in% ind_batch
      )
    )
    ind_batch <- params$ind_batch_list[[ind_batch]]
    ex_tbl_all <- .get_ex_list( # nolint
      .data = gs, ind_batch = ind_batch,
      pop = "root",
      cut = purrr::map_chr(marker_list, function(x) x$cut), high = NULL,
      data_name = purrr::map
    ) |>
      dplyr::bind_rows()
    ex_uns <- ex_tbl_all |> dplyr::filter(stim == "uns") # nolint
    ex_stim <- ex_tbl_all |> dplyr::filter(.data$ind == .env$i) # nolint
    cut_alt <- ifelse(chnl == "Ho165Di", "Nd146Di", "Ho165Di")
    ex_uns <- ex_uns[ex_uns[[chnl]] > min(ex_uns[[chnl]]) |
      ex_uns[[cut_alt]] > min(ex_uns[[cut_alt]]), ]
    ex_stim <- ex_stim[ex_stim[[chnl]] > min(ex_stim[[chnl]]) |
      ex_stim[[cut_alt]] > min(ex_stim[[cut_alt]]), ]
    cn_vec_uns <- colnames(ex_uns)
    chnl_ind <- which(cn_vec_uns == chnl)
    cn_vec_uns[chnl_ind] <- "x"
    cut_alt_ind <- which(cn_vec_uns == cut_alt)
    cn_vec_uns[cut_alt_ind] <- "y"
    colnames(ex_uns) <- cn_vec_uns

    cn_vec_stim <- colnames(ex_stim)
    chnl_ind <- which(cn_vec_stim == chnl)
    cn_vec_stim[chnl_ind] <- "x"
    cut_alt_ind <- which(cn_vec_stim == cut_alt)
    cn_vec_stim[cut_alt_ind] <- "y"
    colnames(ex_stim) <- cn_vec_stim

    ex_plot <- ex_uns |>
      dplyr::bind_rows(ex_stim)
    ggplot(ex_plot, aes(x = x, y = y)) + # nolint
      geom_hex() + # nolint
      scale_fill_viridis_c(trans = "log10") + # nolint
      facet_wrap(~stim, ncol = 1) + # nolint
      geom_vline( # nolint
        xintercept = cp_vec[sel_tbl$grp[sel_tbl$ind == i]],
        col = "red", size = 2
      ) +
      coord_cartesian( # nolint
        xlim = c(0, 7),
        ylim = c(0, 7)
      ) +
      labs(x = chnl, y = cut_alt) # nolint
  })

  cowplot::plot_grid(
    plotlist = plot_list, ncol = 2,
    labels = names(cp_vec)
  )
}
