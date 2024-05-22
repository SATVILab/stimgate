plot_cp_all <- function(data,
                        gate_name = NULL,
                        chnl_base = NULL,
                        params = NULL,
                        chnl = NULL,
                        pop_gate,
                        chnl_lab = NULL,
                        ind_in_batch_lab,
                        ind_in_batch_gate,
                        fcs, data_name,
                        ind_in_batch_uns,
                        ind_batch_list,
                        path_project) {
  if (is.null(params)) {
    if (is.null(chnl)) {
      stop("chnl must be specified if params is NULL in plot_cp_all")
    }
    if (is.null(chnl_lab)) {
      chnl_lab <- .get_labs(
        data = data[[1]],
        cut = chnl
      )
    }
    params <- list(
      pop_gate = pop_gate,
      chnl_lab = chnl_lab,
      ind_in_batch_lab_vec = ind_in_batch_lab,
      ind_in_batch_gate = ind_in_batch_gate,
      data_name = data_name,
      fcs = fcs,
      ind_in_batch_uns = ind_in_batch_uns,
      ind_batch_list = ind_batch_list,
      data = data
    )
  }

  gate_tbl <- purrr::map_df(names(params$chnl_lab), function(chnl_curr) {
    # get base directory
    params[["cut"]] <- chnl_curr
    dir_base <- stimgate_dir_base_create(
      dir_base_init = path_project,
      params = params
    )
    # get stats tbl
    gate_tbl <- readRDS(file.path(dir_base, "gate_tbl.rds"))

    if (!is.null(gate_name)) {
      gate_tbl <- gate_tbl |>
        dplyr::filter(gate_name %in% .env$gate_name)
    }

    gate_tbl |>
      # dplyr::filter(.data$gate_name == .env$gate_name) |>
      dplyr::mutate(chnl = chnl_curr, marker = params$chnl_lab[chnl_curr]) |>
      dplyr::select(
        chnl, marker, gate_name, batch,
        ind, gate, gate_cyt, gate_single
      )
  })

  if (!"gate_single" %in% colnames(gate_tbl)) {
    gate_tbl <- gate_tbl |> dplyr::mutate(gate_single = gate)
  }
  params[["cut"]] <- gate_tbl$chnl[1]
  dir_base <- stimgate_dir_base_create(
    dir_base_init = path_project,
    params = params
  )
  dir_base <- stringr::str_sub(
    dir_base,
    end = -stringr::str_length(gate_tbl$marker[1]) - 2
  )

  chnl_vec <- names(params$chnl_lab)

  if (is.null(chnl_base)) chnl_base <- chnl_vec


  # all data
  purrr::walk(params$ind_batch_list, function(ind_batch) {
    purrr::walk(ind_batch, function(ind) {
      # print(ind)
      ind_in_batch <- ind %% length(params$ind_in_batch_lab_vec)

      # return if ind in batch is the last one, as that is the unstim ind
      if (ind_in_batch == 0) {
        return(NULL)
      }

      # get stim
      stim <- params$ind_in_batch_lab_vec[[ind_in_batch]]

      # get expression dataframe
      ex <- .get_ex(
        data = data[[ind]], pop = params$pop_gate,
        cut = names(params$chnl_lab),
        high = NULL, ind = ind,
        is_uns = FALSE, stim = stim,
        ind_in_batch = ind_in_batch, data_name = params$data_name
      )

      ind_uns <- (ind %/% 5 * 5) + 5
      ex_uns <- .get_ex(
        data = data[[ind_uns]], pop = params$pop_gate,
        cut = names(params$chnl_lab),
        high = NULL, ind = ind,
        is_uns = FALSE, stim = stim,
        ind_in_batch = ind_in_batch,
        data_name = params$data_name
      )
      gate_tbl_ind <- gate_tbl |>
        dplyr::filter(.data$ind == .env$ind) #|>
      # dplyr::mutate(gate_single = pmax(gate, gate_single))

      gate_lab_vec <- .get_gate_lab_vec(gate_tbl = gate_tbl)

      dir_save <- file.path(
        dir_base, paste0(names(params$chnl_lab), collapse = "_"), "p_2d"
      )
      if (!dir.exists(dir_save)) {
        dir.create(dir_save, recursive = TRUE)
      }

      purrr::walk(chnl_base, function(chnl_curr) {
        # print(chnl_curr)
        plot_list <- purrr::map(
          setdiff(chnl_vec, chnl_curr), function(chnl_alt_curr) {
            # print(chnl_alt_curr)
            purrr::map(1:2, function(i) {
              # print(i)
              stim_lab <- ifelse(i == 1, "stim", "unstim")
              if (i == 1) {
                ex_plot <- ex
              } else {
                ex_plot <- ex_uns
              }
              range_ex <- range(ex[[chnl_curr]])
              range_ex_uns <- range(ex_uns[[chnl_curr]])
              range_x <- c(NA, NA)
              # get maximum gates here.
              max_gate_x <- gate_tbl_ind |>
                dplyr::filter(chnl == chnl_curr) |>
                tidyr::pivot_longer(
                  names_to = "type", values_to = "gate",
                  gate:gate_single
                ) |>
                dplyr::pull("gate") |>
                max()

              range_x[1] <- min(
                range_ex[1], range_ex_uns[1], max_gate_x - diff(range_ex) * 0.02
              )
              range_x[2] <- max(
                range_ex[2], range_ex_uns[2], max_gate_x + diff(range_ex) * 0.02
              )

              range_ex <- range(ex[[chnl_alt_curr]])
              range_ex_uns <- range(ex_uns[[chnl_alt_curr]])
              max_gate_y <- gate_tbl_ind |>
                dplyr::filter(chnl == chnl_alt_curr) |>
                tidyr::pivot_longer(
                  names_to = "type", values_to = "gate",
                  gate:gate_single
                ) |>
                dplyr::pull("gate") |>
                max()
              range_y <- c(NA, NA)
              range_y[1] <- min(
                range_ex[1], range_ex_uns[1], max_gate_y - diff(range_ex) * 0.02
              )
              range_y[2] <- max(
                range_ex[2], range_ex_uns[2], max_gate_y + diff(range_ex) * 0.02
              )

              line_tbl <- purrr::map_df(unique(gate_tbl_ind$gate_name), function(gn) {
                gt_tbl_chnl_curr_gn <- gate_tbl_ind |>
                  dplyr::filter(
                    gate_name == gn,
                    chnl == chnl_curr
                  )
                gt_tbl_chnl_alt_gn <- gate_tbl_ind |>
                  dplyr::filter(
                    gate_name == gn,
                    chnl == chnl_alt_curr
                  )
                # single to cyt+
                single_x_x <- rep(gt_tbl_chnl_curr_gn$gate_single, 2)
                single_x_y <- c(range_y[1], gt_tbl_chnl_alt_gn$gate_cyt)
                single_y_y <- rep(gt_tbl_chnl_alt_gn$gate_single, 2)
                single_y_x <- c(range_x[1], gt_tbl_chnl_curr_gn$gate_cyt)

                #
                sc_x_x <- c(range_x[2], gt_tbl_chnl_curr_gn$gate)
                sc_x_y <- rep(gt_tbl_chnl_alt_gn$gate_cyt, 2)
                sc_y_x <- rep(gt_tbl_chnl_curr_gn$gate_cyt, 2)
                sc_y_y <- c(range_y[2], gt_tbl_chnl_alt_gn$gate)

                mid_x_x <- rep(gt_tbl_chnl_curr_gn$gate, 2)
                mid_x_y <- c(gt_tbl_chnl_alt_gn$gate, gt_tbl_chnl_alt_gn$gate_cyt)
                mid_y_y <- rep(gt_tbl_chnl_alt_gn$gate, 2)
                mid_y_x <- c(
                  gt_tbl_chnl_curr_gn$gate, gt_tbl_chnl_curr_gn$gate_cyt
                )

                x_vec <- c(
                  single_x_x, single_y_x, sc_x_x, sc_y_x, mid_x_x, mid_y_x
                )
                y_vec <- c(
                  single_x_y, single_y_y, sc_x_y, sc_y_y, mid_x_y, mid_y_y
                )
                grp_vec <- paste0(
                  gn,
                  rep(c("single_x", "single_y", "sc_x", "sc_y", "mid_x", "mid_y"), each = 2)
                )
                tibble::tibble(x = x_vec, y = y_vec, gate_name = gn, grp = grp_vec)
              })

              vline_tbl <- purrr::map_df(unique(gate_tbl_ind$gate_name), function(gn) {
                gt_tbl_chnl_curr_gn <- gate_tbl_ind |>
                  dplyr::filter(
                    gate_name == gn,
                    chnl == chnl_curr
                  )
                tibble::tibble(x = gt_tbl_chnl_curr_gn$gate, gate_name = gn)
              })
              hline_tbl <- purrr::map_df(
                unique(gate_tbl_ind$gate_name), function(gn) {
                  gt_tbl_chnl_alt_gn <- gate_tbl_ind |>
                    dplyr::filter(
                      gate_name == gn,
                      chnl == chnl_alt_curr
                    )
                  tibble::tibble(y = gt_tbl_chnl_alt_gn$gate, gate_name = gn)
                }
              )
              p <- ggplot(ex_plot |>
                dplyr::filter(
                  !!ensym(chnl_alt_curr) > min(!!ensym(chnl_alt_curr)) |
                    !!ensym(chnl_curr) > min(!!ensym(chnl_curr))
                )) +
                cowplot::theme_cowplot(font_size = 5) +
                geom_hex(
                  aes(x = !!ensym(chnl_curr), y = !!ensym(chnl_alt_curr))
                ) +
                scale_fill_viridis_c(trans = "log10") +
                geom_vline(aes(xintercept = x, col = gate_name),
                  size = 0.6,
                  data = vline_tbl, linetype = "dotted"
                ) +
                geom_hline(aes(yintercept = y, col = gate_name),
                  size = 0.6,
                  data = hline_tbl, linetype = "dotted"
                ) +
                geom_line(aes(x = x, y = y, group = grp, col = gate_name),
                  size = 0.6,
                  data = line_tbl
                ) +
                labs(
                  x = params$chnl_lab[chnl_curr],
                  y = params$chnl_lab[chnl_alt_curr],
                  title = paste0(params$chnl_lab[chnl_alt_curr], "- ", stim_lab)
                ) +
                coord_cartesian(xlim = range_x, ylim = range_y) +
                scale_colour_brewer(palette = "Set1", labels = gate_lab_vec) +
                theme(legend.title = element_blank())
              if (!(chnl_alt_curr == setdiff(chnl_vec, chnl_curr)[1] & i == 2)) {
                p <- p + guides(colour = "none")
              }
              p
            })
          }
        ) |>
          flatten()

        p <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
        dir_save_fn <- file.path(
          dir_save, paste0(
            params$chnl_lab[chnl_curr], "-", params$fcs[ind], "-", stim, ".png"
          )
        )
        if (file.exists(dir_save_fn)) {
          file.remove(dir_save_fn)
        }
        cowplot::ggsave2(
          dir_save_fn,
          units = "cm", width = 10, height = length(chnl_vec) * 4.25
        )
      })
    })
  })

  # singles only
  invisible(TRUE)
}
