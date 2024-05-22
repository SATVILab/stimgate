#' @title Plot gate, statistic and performance and create tables thereof for gates derived from a single population
plot_cp <- function(gate_tbl,
                    params,
                    pop_sub = NULL,
                    plot = TRUE,
                    path_project) {
  # plot histogram of gates
  if (plot || TRUE) {
    # ========================
    # Plot base gates
    # ========================
    gate_lab_vec <- .get_gate_lab_vec(gate_tbl = gate_tbl) # nolint

    gate_name_vec <- unique(gate_tbl$gate_name) |> sort()
    gate_name_level_vec <- gate_name_vec[
      !stringr::str_detect(gate_name_vec, "man")
    ]
    gate_name_level_vec <- c(
      gate_name_level_vec,
      gate_name_vec[stringr::str_detect(gate_name_vec, "man")]
    )
    plot_tbl <- gate_tbl |>
      dplyr::mutate(gate_name = factor(gate_tbl$gate_name,
        levels = gate_name_level_vec
      ))

    len <- length(params$data)
    len_uns <- round(len / params$ind_in_batch_uns)
    len_stim <- len - len_uns
    max_cells_per_uns_sample <- round(1e6 / len_uns)
    max_cells_per_stim_sample <- round(1e6 / len_stim)
    ex_tbl_large <- purrr::map_df(seq_along(params$data), function(i) {
      is_uns <- (i %% params$ind_in_batch_uns) == 0
      fr <- flowWorkspace::gh_pop_get_data(params$data[[i]])
      ex <- flowCore::exprs(fr) |> tibble::as_tibble()
      filter_vec <- ex[[params$cut]] > min(ex[[params$cut]])
      ex <- ex[filter_vec, ]
      max_cells_per_sample <- ifelse(is_uns, max_cells_per_uns_sample,
        max_cells_per_stim_sample
      )
      n_row_sample <- min(nrow(ex), max_cells_per_sample)
      if (n_row_sample == nrow(ex)) {
        return(ex[, params$cut] |>
          dplyr::mutate(is_uns = ifelse(is_uns, "uns", "stim")))
      }
      ex[sample.int(nrow(ex), n_row_sample), params$cut] |>
        dplyr::mutate(stim = ifelse(is_uns, "uns", "stim"))
    })

    dens_large_uns <- density(ex_tbl_large |>
      dplyr::filter(is_uns == "uns") |> # nolint
      dplyr::pull(params$cut))
    dens_tbl_large_uns <- tibble::tibble(
      x = dens_large_uns$x,
      y = dens_large_uns$y
    ) |>
      dplyr::mutate(y = y / max(y)) |>
      dplyr::mutate(stim = "uns")

    dens_large_stim <- density(ex_tbl_large |>
      dplyr::filter(is_uns == "stim") |>
      dplyr::pull(params$cut))
    dens_tbl_large_stim <- tibble::tibble(
      x = dens_large_stim$x,
      y = dens_large_stim$y
    ) |>
      dplyr::mutate(y = y / max(y)) |>
      dplyr::mutate(stim = "stim")

    dens_tbl_large <- dplyr::bind_rows(
      dens_tbl_large_uns,
      dens_tbl_large_stim
    )

    plot_list <- purrr::map(gate_name_level_vec, function(gate_name_curr) {
      plot_tbl_curr <- plot_tbl |>
        dplyr::filter(.data$gate_name == .env$gate_name_curr)
      dens_gates <- density(plot_tbl_curr$gate,
        kernel = "rectangular",
        adjust = 0.75
      )
      dens_gates_tbl <- tibble::tibble(
        x = dens_gates$x,
        y = dens_gates$y
      ) |>
        dplyr::mutate(y = y / max(y)) |>
        dplyr::mutate(stim = "gates")
      plot_tbl_lrg <- dens_tbl_large |>
        dplyr::bind_rows(dens_gates_tbl)
      ggplot(plot_tbl_lrg, aes(
        x = x, y = y,
        col = stim
      )) +
        cowplot::background_grid(major = "x", minor = "x") +
        geom_line(size = 1) +
        geom_vline(
          xintercept = min(plot_tbl_curr$gate),
          size = 1
        ) +
        geom_vline(
          xintercept = max(plot_tbl_curr$gate),
          size = 1
        ) +
        labs(title = gate_lab_vec[gate_name_curr]) +
        theme(legend.title = element_blank()) +
        scale_color_manual(values = c(
          "stim" = "orange",
          "uns" = "dodgerblue",
          "gates" = "springgreen"
        )) +
        labs(
          x = params$chnl_lab[params$cut],
          y = "Density"
        ) +
        coord_cartesian(xlim = c(min(plot_tbl_lrg$x), max(plot_tbl$gate)))
    })
    p_gates <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)

    dir_base <- stimgate_dir_base_create(
      dir_base_init = path_project,
      params = params
    )
    dir_save <- file.path(
      dir_base, "gating_plots",
      "gate", "plot_hist_gates"
    )
    if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
    fn <- "gates-base.png"
    cowplot::ggsave2(
      filename = file.path(file.path(dir_save, fn)),
      height = 10, width = 10
    )

    # ========================
    # Plot cyt-pos gates
    # ========================

    # gates
    # -----------------

    purrr::walk(unique(gate_tbl$gate_name), function(gn) {
      p <- ggplot(
        gate_tbl |>
          dplyr::filter(gate_name == gn) |>
          dplyr::mutate(diff = gate - gate_cyt),
        aes(x = gate, y = diff)
      ) +
        cowplot::theme_cowplot(font_size = 20) +
        cowplot::background_grid(major = "xy") +
        geom_hline(yintercept = 0, size = 1) +
        geom_point() +
        geom_smooth() +
        facet_wrap(gate_name ~ marker, ncol = 2)

      dir_base <- stimgate_dir_base_create(
        dir_base_init = path_project,
        params = params
      )
      dir_save <- file.path(
        dir_base, "gating_plots",
        "gate", "plot_hist_gates"
      )
      if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
      fn <- paste0("gates-cyt_pos-", gn, ".png")
      cowplot::ggsave2(
        filename = file.path(file.path(dir_save, fn)),
        height = 10, width = 10
      )
    })

    # ========================
    # Plot cyt_single gates
    # ========================

    # gates
    # -----------------

    purrr::walk(unique(gate_tbl$gate_name), function(gn) {
      p <- ggplot(
        gate_tbl |>
          dplyr::filter(gate_name == gn) |>
          dplyr::mutate(diff = gate - gate_single),
        aes(x = gate, y = diff)
      ) +
        cowplot::theme_cowplot(font_size = 20) +
        cowplot::background_grid(major = "xy") +
        geom_hline(yintercept = 0, size = 1) +
        geom_point() +
        geom_smooth() +
        facet_wrap(gate_name ~ marker, ncol = 2)

      dir_base <- stimgate_dir_base_create(
        dir_base_init = path_project,
        params = params
      )
      dir_save <- file.path(
        dir_base, "gating_plots",
        "gate", "plot_hist_gates"
      )
      if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
      fn <- paste0("gates-single-", gn, ".png")
      cowplot::ggsave2(
        filename = file.path(file.path(dir_save, fn)),
        height = 10, width = 10
      )
    })
  }

  # =========================
  # Preparation
  # =========================

  # get labels for gates
  # debugonce(.get_gate_lab_vec)
  # gate_lab_vec <- .get_gate_lab_vec(gate_tbl = gate_tbl_2)
  gate_lab_vec <- .get_gate_lab_vec(gate_tbl = gate_tbl)


  # colour vector
  col <- RColorBrewer::brewer.pal(length(gate_lab_vec), name = "Set1") |>
    stats::setNames(names(gate_lab_vec))

  # =========================
  # Gates
  # =========================

  gate_tbl <- gate_tbl |>
    dplyr::mutate(
      gate_type = purrr::map_chr(
        gate_name, function(gn) stringr::str_split(gn, "_")[[1]][1]
      ),
      gate_combn = gate_name |>
        stringr::str_remove("_adj") |>
        stringr::str_remove("_clust") |>
        stringr::str_remove(gate_type) |>
        stringr::str_remove("_")
    )

  # plots
  if (plot) {
    print("getting gate plots")
    plot_gate(
      gate_tbl = gate_tbl, params = params,
      pop_sub = pop_sub,
      gate_lab = gate_lab_vec, col = col
    )
  }


  # =========================
  # Statistics
  # =========================

  # plots
  print("getting stats plots")
  if (plot) {
    # need to add changes to get_gate_stats here
    #
    # gate_stats_tbl <- purrr::map_df(.get_gate_stats(params = params,
    #                                                         gate_tbl_calc = gate_tbl_gn |> dplyr::filter(gate_use == 'gate'),
    #                                                  gate_tbl_filter = gate_tbl |> dplyr::filter(gate_use == 'gate',
    #                                                                                        gate_name == gate_name_curr),
    #                                                  pop_sub = NULL,
    #
    #                                                  single_pos = TRUE,
    #                                                  gate_type_cyt_pos = ifelse(params$calc_cyt_pos_gates,
    #                                                                             'cyt', 'base'),
    #                                                  gate_type_single_pos = 'base',
    #                                                  gate_type_single_pos_calc = 'base',
    #                                                  path_project = path_project
    # )

    # .plot_gate_stats(gate_stats = gate_stats_tbl,
    #                                 params = params,
    #                                 col = col,
    #                                 gate_lab = gate_lab_vec)
  }

  # =========================
  # Performance
  # =========================

  # table
  if (!is.null(params$pop_man_sub)) {
    ## need to add changes to get_gate_stats here
    #
    # print("getting perf")
    # gate_perf_tbl <- .get_gate_perf(gate_stats = gate_stats_tbl,
    #                                                params = params,
    #                                 path_project = path_project)

    # plots
    # print("getting perf plots")
    # .plot_gate_perf(gate_perf = gate_perf_tbl,
    #                                params = params,
    #                                gate_lab = gate_lab_vec)
  }
  invisible(TRUE)
}


.get_cp_adj_tbl <- function(gate_stats_tbl, gate_quant, gate_tbl_ctrl) {
  gate_name_vec <- setdiff(
    unique(gate_stats_tbl$gate_name),
    "man_no"
  )

  gate_tbl <- gate_stats_tbl

  if (!is.null(gate_tbl_ctrl)) {
    gate_tbl_adj_ctrl <- purrr::map_df(gate_name_vec, function(gate_name_curr) {
      if (!is.null(names(gate_quant))) {
        gate_combn_curr <- gate_stats_tbl |>
          dplyr::filter(gate_name == gate_name_curr) |>
          dplyr::slice(1) |>
          dplyr::pull("gate_combn")
        gate_quant_curr <- gate_quant[[gate_combn_curr]]
      } else {
        gate_quant_curr <- gate_quant
      }
      purrr::map_df(unique(gate_stats_tbl$boot_ind), function(boot_ind_curr) {
        gate_tbl_curr <- gate_stats_tbl |>
          dplyr::filter(
            .data$gate_name == .env$gate_name_curr,
            .data$boot_ind == .env$boot_ind_curr
          )
        gate_tbl_curr_pos <- gate_tbl_curr |>
          dplyr::filter(freq_bs > 0)

        if (nrow(gate_tbl_curr_pos) == 0) {
          gate_tbl_curr <- gate_tbl_curr |>
            dplyr::mutate(
              gate_name = paste0(gate_name, "_adj_ctrl"),
              gate_combn = paste0(gate_combn, "_adj_ctrl")
            )
          return(gate_tbl_curr)
        }
        freq_bs_bounds <- quantile(gate_tbl_curr_pos$freq_bs, c(0.6, 0.9))
        count_bs_bounds <- quantile(gate_tbl_curr_pos$count_bs, c(0.6, 0.9))
        if (length(freq_bs_bounds) < 2 | length(count_bs_bounds) < 2) {
          freq_bs_bounds <- quantile(gate_tbl_curr_pos$freq_bs, c(0.3, 0.9))
          count_bs_bounds <- quantile(gate_tbl_curr_pos$count_bs, c(0.3, 0.9))
        }
        gate_tbl_curr_pos <- gate_tbl_curr_pos |>
          dplyr::filter(
            freq_bs > freq_bs_bounds[1],
            freq_bs <= freq_bs_bounds[2],
            count_bs > count_bs_bounds[1],
            count_bs <= count_bs_bounds[2]
          )
        if (nrow(gate_tbl_curr_pos) == 0) {
          gate_tbl_curr <- gate_tbl_curr |>
            dplyr::mutate(
              gate_name = paste0(gate_name, "_adj_ctrl"),
              gate_combn = paste0(gate_combn, "_adj_ctrl")
            )
          return(gate_tbl_curr)
        }
        gate_bounds <- quantile(gate_tbl_curr_pos$gate, gate_quant_curr)

        # remove those that might have had a batch effect
        if (!is.null(gate_tbl_ctrl)) {
          gate_tbl_ctrl <- gate_tbl_ctrl |>
            dplyr::rename(gate_ctrl = gate) |>
            dplyr::select(-c(
              batch, gate_name, gate_type, gate_combn,
              boot, gate_use
            ))

          gate_tbl_curr <- gate_tbl_curr |>
            dplyr::left_join(gate_tbl_ctrl,
              by = c("ind", "boot_ind")
            )

          gate_tbl_curr_shift <- gate_tbl_curr |>
            dplyr::filter(gate_ctrl > max(gate_bounds))

          # gate_tbl_curr_shift <- gate_tbl_curr_shift |>
          #  dplyr::mutate(gate = pmax(gate, gate_ctrl))

          gate_tbl_curr <- gate_tbl_curr |>
            dplyr::filter(gate_ctrl <= max(gate_bounds)) |>
            dplyr::mutate(
              gate = pmax(gate, gate_bounds[1]),
              gate = pmin(gate, gate_bounds[2])
            ) |>
            dplyr::bind_rows(gate_tbl_curr_shift) |>
            dplyr::select(-gate_ctrl)
        } else {
          gate_tbl_curr <- gate_tbl_curr |>
            dplyr::mutate(
              gate = pmax(gate, gate_bounds[1]),
              gate = pmin(gate, gate_bounds[2])
            )
        }

        gate_tbl_curr <- gate_tbl_curr |>
          dplyr::mutate(
            gate_name = paste0(gate_name, "_adj_ctrl"),
            gate_combn = paste0(gate_combn, "_adj_ctrl")
          )
      })
    })

    gate_tbl <- gate_tbl |>
      dplyr::bind_rows(gate_tbl_adj_ctrl)
  }

  gate_tbl_ctrl <- NULL
  gate_tbl_adj <- purrr::map_df(gate_name_vec, function(gate_name_curr) {
    if (!is.null(names(gate_quant))) {
      gate_combn_curr <- gate_stats_tbl |>
        dplyr::filter(gate_name == gate_name_curr) |>
        dplyr::slice(1) |>
        dplyr::pull("gate_combn")
      gate_quant_curr <- gate_quant[[gate_combn_curr]]
    } else {
      gate_quant_curr <- gate_quant
    }
    purrr::map_df(unique(gate_stats_tbl$boot_ind), function(boot_ind_curr) {
      gate_tbl_curr <- gate_stats_tbl |>
        dplyr::filter(
          .data$gate_name == .env$gate_name_curr,
          .data$boot_ind == .env$boot_ind_curr
        )
      gate_tbl_curr_pos <- gate_tbl_curr |>
        dplyr::filter(freq_bs > 0)

      if (nrow(gate_tbl_curr_pos) == 0) {
        gate_tbl_curr <- gate_tbl_curr |>
          dplyr::mutate(
            gate_name = paste0(gate_name, "_adj"),
            gate_combn = paste0(gate_combn, "_adj")
          )
        return(gate_tbl_curr)
      }
      freq_bs_bounds <- quantile(gate_tbl_curr_pos$freq_bs, c(0.6, 0.9))
      count_bs_bounds <- quantile(gate_tbl_curr_pos$count_bs, c(0.6, 0.9))
      if (length(freq_bs_bounds) < 2 | length(count_bs_bounds) < 2) {
        freq_bs_bounds <- quantile(gate_tbl_curr_pos$freq_bs, c(0.3, 0.9))
        count_bs_bounds <- quantile(gate_tbl_curr_pos$count_bs, c(0.3, 0.9))
      }
      gate_tbl_curr_pos <- gate_tbl_curr_pos |>
        dplyr::filter(
          freq_bs > freq_bs_bounds[1],
          freq_bs <= freq_bs_bounds[2],
          count_bs > count_bs_bounds[1],
          count_bs <= count_bs_bounds[2]
        )
      if (nrow(gate_tbl_curr_pos) == 0) {
        gate_tbl_curr <- gate_tbl_curr |>
          dplyr::mutate(
            gate_name = paste0(gate_name, "_adj"),
            gate_combn = paste0(gate_combn, "_adj")
          )
        return(gate_tbl_curr)
      }
      gate_bounds <- quantile(gate_tbl_curr_pos$gate, gate_quant_curr)

      # remove those that might have had a batch effect
      if (!is.null(gate_tbl_ctrl)) {
        gate_tbl_ctrl <- gate_tbl_ctrl |>
          dplyr::rename(gate_ctrl = gate) |>
          dplyr::select(-c(
            batch, gate_name, gate_type, gate_combn,
            boot, gate_use
          ))

        gate_tbl_curr <- gate_tbl_curr |>
          dplyr::left_join(gate_tbl_ctrl,
            by = c("ind", "boot_ind")
          )

        gate_tbl_curr_shift <- gate_tbl_curr |>
          dplyr::filter(gate_ctrl > max(gate_bounds))

        gate_tbl_curr <- gate_tbl_curr |>
          dplyr::filter(gate_ctrl <= max(gate_bounds)) |>
          dplyr::mutate(
            gate = pmax(gate, gate_bounds[1]),
            gate = pmin(gate, gate_bounds[2])
          ) |>
          dplyr::bind_rows(gate_tbl_curr_shift) |>
          dplyr::select(-gate_ctrl)
      } else {
        gate_tbl_curr <- gate_tbl_curr |>
          dplyr::mutate(
            gate = pmax(gate, gate_bounds[1]),
            gate = pmin(gate, gate_bounds[2])
          )
      }

      gate_tbl_curr <- gate_tbl_curr |>
        dplyr::mutate(
          gate_name = paste0(gate_name, "_adj"),
          gate_combn = paste0(gate_combn, "_adj")
        )
    })
  })

  gate_tbl <- gate_tbl |>
    dplyr::bind_rows(gate_tbl_adj)

  gate_tbl |>
    dplyr::select(
      gate_name, gate_type, gate_combn, batch, ind,
      gate, boot, boot_ind
    )
}
