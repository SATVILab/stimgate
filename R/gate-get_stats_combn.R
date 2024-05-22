#' @title Get counts for cells for combinations of cytokines
#'
#' @inheritParams gate
#'
#' @return A dataframe.
#'
#' @examples
#' get_stats_combn(
#'   data = gs, pop_gate = "/CD3+/CD4 T Cells",
#'   cut = c("Ho165Di", "Nd146Di"),
#'   ind_in_batch_lab_vec = c("1" = "ebv", "2" = "mtb", "3" = "uns"),
#'   ind_in_batch_gate = 1:3,
#'   gate_name = "locb5_no"
#' )
get_stats_combn_cp <- function(data,
                               pop_gate,
                               cut,
                               ind_in_batch_lab_vec,
                               gate_name,
                               ind_in_batch_gate,
                               ind_in_batch_uns,
                               debug = FALSE,
                               use_cyt_pos = TRUE,
                               use_single_pos = TRUE,
                               path_project) {
  # ==============================
  # Preparation
  # ==============================

  # Create params list
  # -----------------------------

  # dataset_name
  data_name <- deparse(substitute(data))

  # chnl_lab
  chnl_lab_vec <- .get_labs(
    data = data[[1]],
    cut = cut
  )

  # params object
  params <- list(
    pop_gate = pop_gate,
    chnl_lab = chnl_lab_vec,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    ind_in_batch_gate = ind_in_batch_gate,
    data_name = data_name
  )

  # ====================================
  # Get gates
  # ====================================

  gate_tbl <- purrr::map_df(cut, function(cut_curr) {
    # get base directory
     params[["cut"]] <- cut_curr
    dir_base <- stimgate_dir_base_create(
      dir_base_init = path_project,
      params = params
    )
    # get stats tbl
    gate_stats <- readRDS(file.path(dir_base, "gate_tbl.rds"))

    gate_stats |>
      dplyr::filter(.data$gate_name == .env$gate_name) |>
      dplyr::mutate(cut = cut_curr, marker = chnl_lab_vec[cut_curr]) |>
      dplyr::select(cut, marker, everything())
  })

  if (single_only) {
    stat_tbl <- purrr::map(seq_along(data), function(ind) {
      # get which ind in batch
      if (debug) print(ind)

      ind_in_batch <- ind %% length(ind_in_batch_lab_vec)

      # return if ind in batch is the last one, as that is the unstim ind
      if (ind_in_batch == 0) {
        return(NULL)
      }

      # get stim
      stim <- ind_in_batch_lab_vec[[ind_in_batch]]

      # get expression dataframe
      ex <- .get_ex(
        data = data[[ind]], pop = pop_gate,
        cut = cut, high = NULL, ind = ind,
        is_uns = FALSE, stim = stim,
        ind_in_batch = ind_in_batch, data_name = data_name
      )

      gate_tbl_ind <- gate_tbl |> dplyr::filter(.data$ind == .env$ind)

      inc_vec <- rep(FALSE, nrow(ex))

      purrr::map_df(gate_name_vec, function(gn) {
        gate_tbl_ind_gn <- gate_tbl_ind |> dplyr::filter(.data$gate_name == gn)

        count_vec <- purrr::map_dbl(cut, function(chnl_curr) {
          pos_ind_vec_but_single_pos_curr <-
            .get_pos_ind_but_single_pos_for_one_cyt(
              ex = ex,
              gate_tbl = gate_tbl_ind_gn,
              chnl = chnl_curr
            )
          ex_neg_but_curr <- ex[
            !pos_ind_vec_but_single_pos_curr, ,
            drop = FALSE
          ]

          sum(ex_neg_but_curr[[chnl_curr]] >
            gate_tbl_ind_gn$gate[gate_tbl_ind_gn$chnl == chnl_curr])
        })

        ind_uns <- (ind %/% length(ind_in_batch_lab_vec)) *
          length(ind_in_batch_lab_vec) +
          ind_in_batch_uns
        ex_uns <- .get_ex(
          data = data[[ind_uns]], pop = pop_gate,
          cut = cut, high = NULL, ind = ind_uns,
          is_uns = TRUE, stim = "uns",
          ind_in_batch = ind_in_batch, data_name = data_name
        )

        count_vec_uns <- purrr::map_dbl(cut, function(chnl_curr) {
          pos_ind_vec_but_single_pos_curr <-
            .get_pos_ind_but_single_pos_for_one_cyt(
              ex = ex_uns,
              gate_tbl = gate_tbl_ind_gn,
              chnl = chnl_curr
            )
          ex_neg_but_curr <- ex_uns[
            !pos_ind_vec_but_single_pos_curr, ,
            drop = FALSE
          ]

          sum(
            ex_neg_but_curr[[chnl_curr]] >
              gate_tbl_ind_gn$gate[gate_tbl_ind_gn$chnl == chnl_curr]
          )
        })

        tibble::tibble(ind = ind, gate_name = gn, count = count_vec, count_vec_uns)
      })
    })
  }

  # ====================================
  # Apply gates
  # ====================================

  ind <- 1
  stat_tbl <- purrr::map(seq_along(data), function(ind) {
    .get_stat_combn <- function(ex, cut, gate_tbl, ind, use_cyt_pos, use_single_pos) {
      if (!use_cyt_pos & !use_single_pos) {
        n_cell <- nrow(ex)

        pos_tbl <- purrr::map(
          cut,
          function(cut_curr) {
            stats::setNames(
              tibble::tibble(x = rep(0, n_cell)), cut_curr
            )
          }
        ) |>
          bind_cols()

        for (cut_curr in cut) {
          inc_vec <- ex[[cut_curr]] > gate_tbl[
            gate_tbl$cut == cut_curr & gate_tbl$ind == ind, "gate",
            drop = TRUE
          ]
          pos_tbl[[cut_curr]] <- inc_vec |> as.numeric()
        }

        ind_list <- purrr::map(seq_along(cut), function(i) 0:1)
        combn_mat <- expand.grid(ind_list)
        all_neg_ind_vec <- purrr::map_lgl(1:nrow(combn_mat), function(i) {
          purrr::map_lgl(
            1:ncol(combn_mat), function(j) !combn_mat[i, j, drop = TRUE]
          ) |>
            all()
        })
        combn_mat <- combn_mat[!all_neg_ind_vec, ]
        combn_tbl <- combn_mat |>
          tibble::as_tibble() |>
          stats::setNames(cut) |>
          dplyr::mutate(count = NA_integer_)
        for (i in 1:nrow(combn_tbl)) {
          ind_vec <- rep(TRUE, nrow(pos_tbl))
          for (j in 1:ncol(combn_mat)) {
            ind_vec <- ind_vec & pos_tbl[[j]] == combn_tbl[i, j, drop = TRUE]
          }
          combn_tbl[i, "count"] <- sum(ind_vec)
        }
      } else {
        n_cell <- nrow(ex)

        pos_tbl <- purrr::map(
          cut,
          function(cut_curr) {
            stats::setNames(tibble::tibble(x = rep(0, n_cell)), cut_curr)
          }
        ) |>
          bind_cols()

        for (cut_curr in cut) {
          inc_vec <- ex[[cut_curr]] > gate_tbl[
            gate_tbl$cut == cut_curr & gate_tbl$ind == ind, "gate",
            drop = TRUE
          ]
          pos_tbl[[cut_curr]] <- inc_vec |> as.numeric()
        }

        ind_list <- purrr::map(seq_along(cut), function(i) 0:1)
        combn_mat <- expand.grid(ind_list)
        all_neg_ind_vec <- purrr::map_lgl(1:nrow(combn_mat), function(i) {
          purrr::map_lgl(
            1:ncol(combn_mat), function(j) !combn_mat[i, j, drop = TRUE]
          ) |>
            all()
        })
        combn_mat <- combn_mat[!all_neg_ind_vec, ]
        combn_tbl <- combn_mat |>
          tibble::e() |>
          stats::setNames(cut) |>
          dplyr::mutate(count = NA_integer_)
        for (i in 1:nrow(combn_tbl)) {
          ind_vec <- rep(TRUE, nrow(pos_tbl))
          for (j in 1:ncol(combn_mat)) {
            ind_vec <- ind_vec & pos_tbl[[j]] == combn_tbl[i, j, drop = TRUE]
          }
          combn_tbl[i, "count"] <- sum(ind_vec)
        }
      }

      stat_tbl <- combn_tbl |>
        dplyr::mutate(
          n_cell = n_cell,
          prop = count / n_cell,
          freq = prop * 100
        )

      stat_tbl
    }

    stat_tbl <- .get_stat_combn(
      ex = ex, cut = cut, gate_tbl = gate_tbl, ind = ind
    )

    # get expression dataframe
    ind_uns <- (ind %/% length(ind_in_batch_lab_vec)) *
      length(ind_in_batch_lab_vec) +
      ind_in_batch_uns
    ex_uns <- .get_ex(
      data = data[[ind_uns]], pop = pop_gate,
      cut = cut, high = NULL, ind = ind_uns,
      is_uns = TRUE, stim = "uns",
      ind_in_batch = ind_in_batch, data_name = data_name
    )

    stat_tbl_uns <- .get_stat_combn(
      ex = ex_uns, cut = cut, gate_tbl = gate_tbl, ind = ind
    )

    stat_tbl <- stat_tbl |>
      dplyr::mutate(
        prop_stim = prop,
        freq_stim = freq,
        count_stim = count,
        n_cell_stim = n_cell,
        prop_uns = stat_tbl_uns$prop,
        freq_uns = stat_tbl_uns$freq,
        count_uns = stat_tbl_uns$count,
        n_cell_uns = stat_tbl_uns$n_cell
      ) |>
      dplyr::mutate(
        prop_bs = prop_stim - prop_uns,
        freq_bs = freq_stim - freq_uns,
        count_bs = round(prop_bs * n_cell_stim)
      )

    stat_tbl |>
      dplyr::mutate(
        pop = pop_gate, gate_name = gate_name,
        batch = ex$batch[1], batch_sh = ex$batch_sh[1],
        fcs = ex$fcs[1], stim = ex$stim[1],
        ind = ex$ind[1], is_uns = ex$is_uns[1],
        ind_in_batch = ex$ind_in_batch[1]
      ) |>
      dplyr::select(
        pop, gate_name, batch, batch_sh, fcs, stim, ind, is_uns,
        ind_in_batch, everything()
      )
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()
}


#' @export
get_stats_combn_nb <- function(data,
                               pop_gate,
                               cut,
                               ind_in_batch_lab_vec,
                               gate_name,
                               ind_in_batch_gate,
                               ind_in_batch_uns,
                               debug = FALSE,
                               path_project) {
  # ==============================
  # Preparation
  # ==============================

  # Create params list
  # -----------------------------

  # dataset_name
  data_name <- deparse(substitute(data))

  # chnl_lab
  chnl_lab_vec <- .get_labs(
    data = data[[1]],
    cut = cut
  )

  # params object
  params <- list(
    pop_gate = pop_gate,
    chnl_lab = chnl_lab_vec,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    ind_in_batch_gate = ind_in_batch_gate,
    data_name = data_name
  )

  # get base directory
  # -----------------------------------
  # get base directory
  params[["cut"]] <- cut[1]
  dir_base <- stimgate_dir_base_create(
    dir_base_init = path_project,
    params = params
  )

  marker <- chnl_lab_vec[cut[1]]
  marker_len <- stringr::str_length(marker)

  dir_base <- stringr::str_sub(
    dir_base,
    end = -marker_len - 2
  )

  # ====================================
  # Get gates
  # ====================================

  gate_tbl <- purrr::map_df(cut, function(cut_curr) {
    # get base directory
    params[["cut"]] <- cut_curr
    dir_base <- stimgate_dir_base_create(
      dir_base_init = path_project,
      params = params
    )
    # get stats tbl
    gate_stats <- readRDS(file.path(dir_base, "stats", "gate_stats_tbl"))

    gate_stats |>
      dplyr::filter(.data$gate_name == .env$gate_name) |>
      dplyr::mutate(cut = cut_curr, marker = chnl_lab_vec[cut_curr]) |>
      dplyr::select(cut, marker, everything())
  })

  gate_tbl_man <- purrr::map_df(cut, function(cut_curr) {
    # get base directory
    params[["cut"]] <- cut_curr
    dir_base <- stimgate_dir_base_create(
      dir_base_init = path_project,
      params = params
    )
    # get stats tbl
    gate_stats <- readRDS(file.path(dir_base, "stats", "gate_stats_tbl"))

    gate_stats |>
      dplyr::filter(.data$gate_name == "man_no") |>
      dplyr::mutate(cut = cut_curr, marker = chnl_lab_vec[cut_curr]) |>
      dplyr::select(cut, marker, everything())
  })

  # ====================================
  # Get probs
  # ====================================

  # get base directory
  count_tbl <- purrr::map(seq_along(data), function(ind) {
    if (debug) print(ind)

    ind_in_batch <- ind %% length(ind_in_batch_lab_vec)

    # return if ind in batch is the last one, as that is the unstim ind
    if (ind_in_batch == 0) {
      return(NULL)
    }

    # get stim condition
    stim <- ind_in_batch_lab_vec[[ind_in_batch]]

    # get expression dataframe
    ex <- .get_ex(
      data = data[[ind]], pop = pop_gate,
      cut = cut, high = NULL, ind = ind,
      is_uns = FALSE, stim = stim,
      ind_in_batch = ind_in_batch, data_name = data_name
    )

    n_cell <- nrow(ex)

    # ==========================================
    # Get cell-specific responder probabilities
    # ==========================================

    probs_tbl <- purrr::map(cut, function(cut_curr) {
      if (debug) print(cut_curr)
      params_curr <- params
      params_curr[["cut"]] <- cut_curr
      dir_base <- stimgate_dir_base_create(
        dir_base_init = path_project,
        params = params_curr
      )
      # wd <- getwd()
      # on.exit(setwd(wd))
      # setwd(file.path(dir_base, "stats", "cell_resp_probs"))
      path_probs <- file.path(
        dir_base, "stats", "probs",
        paste0(ex$batch[1], "_", ex$stim[1], ".rds")
      )
      # fn_probs <- paste0(ex$batch[1], "_", ex$stim[1], ".rds")
      path_probs <- normalizePath(path_probs)
      if (!file.exists(path_probs)) {
        return(NULL)
      }

      ex_exc_min <- ex[ex[[cut_curr]] > (min(ex[[cut_curr]]) + 1), ]
      density_exc_min <- density(ex_exc_min[[cut_curr]])
      dens_tbl <- tibble::tibble(
        x = density_exc_min$x,
        y = density_exc_min$y
      )
      peak <- dens_tbl |>
        dplyr::filter(y == max(y)) |>
        dplyr::pull("x")

      pred_tbl <- readRDS(path_probs)
      pred_tbl <- pred_tbl |> dplyr::mutate(
        cut = cut_curr,
        marker = chnl_lab_vec[cut]
      )
      pred_tbl <- pred_tbl |> dplyr::mutate(expr = cut_stim)

      pred_tbl |>
        dplyr::select(cut, marker, ind_cell, everything()) |>
        dplyr::filter(pred > 0.05) |>
        dplyr::filter(.data$cut_stim > .env$peak)
    }) |>
      purrr::compact() |>
      dplyr::bind_rows()

    probs_tbl_combn_non_zero <- probs_tbl |>
      dplyr::group_by(ind_cell) |>
      dplyr::mutate(pred_non = 1 - pred) |>
      dplyr::summarise(prob = 1 - Reduce(`*`, pred_non)) |>
      dplyr::filter(prob > 0.1) |>
      dplyr::ungroup()

    # ex_resp <- probs_tbl |>
    #  dplyr::left_join(probs_tbl_combn_non_zero |>
    #              dplyr::mutate(prob_resp = prob) |>
    #              dplyr::select(ind_cell, prob_resp),
    #            by = 'ind_cell') |>
    #  dplyr::filter(!is.na(prob_resp))


    prob_resp_lab_vec <- stats::setNames(
      probs_tbl_combn_non_zero$prob,
      probs_tbl_combn_non_zero$ind_cell
    )

    ex_resp <- ex |>
      dplyr::mutate(ind_sample = ind) |>
      dplyr::filter(.data$ind_cell %in% probs_tbl$ind_cell) |>
      dplyr::mutate(
        prob_resp = prob_resp_lab_vec[as.character(.data$ind_cell)]
      ) |>
      dplyr::mutate(type = "resp") |>
      dplyr::select(-cut) |>
      dplyr::select(
        prob_resp, type, ind_sample, stim, ind_cell, Ho165Di:Eu153Di
      ) |>
      dplyr::arrange(type, prob_resp) |>
      tidyr::pivot_longer(
        names_to = "chnl",
        values_to = "expr",
        c(Ho165Di, Gd158Di, Nd146Di)
      ) |>
      dplyr::mutate(marker = chnl_lab_vec[chnl])

    ex_all <- ex |>
      dplyr::mutate(ind_sample = ind) |>
      # dplyr::filter(.data$ind_cell %in% probs_tbl$ind_cell) |>
      dplyr::mutate(
        prob_resp = prob_resp_lab_vec[as.character(.data$ind_cell)]
      ) |>
      dplyr::mutate(type = "all") |>
      dplyr::select(-cut) |>
      dplyr::select(
        prob_resp, type, ind_sample, stim, ind_cell, Ho165Di:Eu153Di
      ) |>
      dplyr::arrange(type, prob_resp) |>
      tidyr::pivot_longer(
        names_to = "chnl",
        values_to = "expr",
        c(Ho165Di, Gd158Di, Nd146Di)
      ) |>
      dplyr::mutate(marker = chnl_lab_vec[chnl]) |>
      dplyr::group_by(chnl) |>
      dplyr::filter(expr > min(expr)) |>
      dplyr::ungroup()

    # ==========================================
    # Plot unweighted density
    # ==========================================

    # plot unweighted density
    ex_plot <- ex_resp |>
      dplyr::bind_rows(ex_all)

    gate_tbl_plot <- gate_tbl |>
      dplyr::filter(.data$ind == .env$ind) |>
      dplyr::bind_rows(gate_tbl_man |>
        dplyr::filter(.data$ind == .env$ind))

    p_dens_resp_vs_all <- ggplot(
      ex_plot,
      aes(x = expr, fill = type)
    ) +
      geom_density(alpha = 0.3) +
      scale_fill_brewer(palette = "Set1") +
      geom_vline(
        data = gate_tbl_plot,
        aes(
          xintercept = gate,
          linetype = gate_name
        ), size = 2
      ) +
      facet_wrap(~marker, scales = "free") +
      geom_rug(data = ex_plot |>
        dplyr::filter(type == "resp"))

    # ==========================================
    # Construct weighted density
    # ==========================================

    min_chnl_vec <- purrr::map_dbl(cut, function(x) min(ex[[x]])) |>
      stats::setNames(cut) + 1

    ex_dens_calc <- ex_resp |>
      dplyr::filter(!is.na(prob_resp)) |>
      dplyr::group_by(chnl) |>
      dplyr::filter(expr > min_chnl_vec[chnl]) |>
      dplyr::ungroup()


    cp_vec <- purrr::map_dbl(cut, function(cut_curr) {
      print(cut_curr)

      if (!cut_curr %in% unique(ex_dens_calc$chnl)) {
        cut <- (ex_all |> dplyr::filter(
          chnl == cut_curr
        ) |>
          dplyr::pull("expr") |>
          median()
        )
        return(cut)
      }
      expr_vec_all <- ex_all |>
        dplyr::filter(chnl == cut_curr) |>
        dplyr::pull("expr")

      expr_vec_resp <- ex_dens_calc |>
        dplyr::filter(chnl == cut_curr) |>
        dplyr::pull("expr")
      # should actually set breaks according to ex_resp
      # rather than ex_all
      hist <- hist(expr_vec_resp,
        plot = FALSE
      )
      break_vec_hist <- hist$breaks
      break_vec_hist[1] <- min(expr_vec_all - 1)
      break_vec_hist[length(break_vec_hist)] <- max(expr_vec_all + 1)
      break_vec_hist <- floor(pmax(break_vec_hist, min(expr_vec_all)))
      break_vec_hist <- ceiling(pmin(break_vec_hist, max(expr_vec_all)))

      cut_vec <- cut(
        ex_resp |>
          dplyr::filter(chnl == cut_curr) |>
          dplyr::pull("expr"),
        breaks = break_vec_hist
      )

      cut_vec <- cut(
        seq(
          break_vec_hist[1],
          max(break_vec_hist)
        ),
        breaks = break_vec_hist
      )
      uni_sort_cut_vec <- unique(cut_vec) |> sort()

      #' @title Get the midpoint of a bin returned by base::cut
      .get_bin_mid <- function(bin) {
        purrr::map_dbl(bin, function(bin_curr) {
          comma_loc <- stringr::str_locate(bin_curr, ",")[1, "start"][[1]]
          num_1 <- stringr::str_sub(bin_curr, 2, comma_loc - 1) |> as.numeric()
          num_2 <- stringr::str_sub(
            bin_curr, comma_loc + 1,
            stringr::str_length(bin_curr) - 1
          ) |> as.numeric()
          mean(c(num_1, num_2))
        })
      }

      .get_bin_lb <- function(bin) {
        purrr::map_dbl(bin, function(bin_curr) {
          comma_loc <- stringr::str_locate(bin_curr, ",")[1, "start"][[1]]
          num_1 <- stringr::str_sub(bin_curr, 2, comma_loc - 1) |> as.numeric()
          # num_2 <- stringr::str_sub(
          #  bin_curr, comma_loc + 1,
          #  stringr::str_length(bin_curr) - 1)
          #    |> as.numeric
          # mean(c(num_1, num_2))
        })
      }

      bin_lab_vec <- stats::setNames(
        seq_along(uni_sort_cut_vec), uni_sort_cut_vec
      )
      bin_lb_lab_vec <- stats::setNames(
        1:(length(break_vec_hist) - 1),
        break_vec_hist[-length(break_vec_hist)]
      )

      bin_length_vec <- purrr::map_dbl(
        seq_along(break_vec_hist[-length(break_vec_hist)]),
        function(i) {
          len <- break_vec_hist[i + 1] - break_vec_hist[i]
        }
      )
      bin_len_lab_vec <- stats::setNames(
        bin_length_vec, seq_len(length(break_vec_hist) - 1)
      )
      ex_hist_weight <- ex_dens_calc |>
        dplyr::filter(!is.na(prob_resp)) |>
        dplyr::filter(chnl == cut_curr) |>
        dplyr::mutate(
          bin = cut(expr, breaks = break_vec_hist),
          bin_mid = purrr::map_dbl(bin, .get_bin_mid),
          bin_lb = purrr::map_dbl(bin, .get_bin_lb)
        )
      ex_hist_weight <- ex_hist_weight |>
        dplyr::mutate(bin_ind = bin_lb_lab_vec[as.character(bin_lb)])

      ex_hist_weight <- ex_hist_weight |>
        dplyr::mutate(bin_len = bin_len_lab_vec[bin_ind])

      ex_hist_ctb <- ex_hist_weight |>
        dplyr::group_by(bin, bin_mid, bin_ind, bin_len) |>
        dplyr::summarise(
          hist_ctb = sum(prob_resp) * bin_len_lab_vec[bin_ind]
        ) |>
        dplyr::ungroup()

      area_under_hist <- sum(ex_hist_ctb$hist_ctb)

      ex_hist_dens <- ex_hist_ctb |>
        dplyr::mutate(dens = hist_ctb / area_under_hist)

      dens_all <- density(ex_all |>
        dplyr::filter(chnl == cut_curr) |>
        dplyr::pull("expr"))

      dens_tbl <- tibble::tibble(
        expr = dens_all$x,
        dens = dens_all$y
      )

      ex_hist_dens <- ex_hist_dens |>
        dplyr::mutate(
          dens_all = purrr::map_dbl(ex_hist_dens$bin_mid, function(bin) {
            bin <- max(bin, min(dens_tbl$expr))
            bin <- min(bin, max(dens_tbl$expr))
            .interp(
              val = bin,
              x = dens_tbl$expr,
              y = dens_tbl$dens
            )
          })
        )


      ex_hist_plot <- dens_tbl |>
        dplyr::mutate(type = "all") |>
        dplyr::bind_rows(ex_hist_dens |>
          dplyr::rename(expr = bin_mid) |>
          dplyr::select(expr, dens) |>
          dplyr::mutate(type = "resp"))

      p_dens_weight_resp_vs_all <- ggplot(
        ex_hist_plot,
        aes(x = expr, y = dens, col = type)
      ) +
        geom_line()

      n_resp_by_cut_alone <- probs_tbl |>
        dplyr::filter(cut == cut_curr) |>
        dplyr::pull("pred") |>
        sum()

      n_resp <- probs_tbl_combn_non_zero$prob |>
        sum()

      prior_prob <- n_resp_by_cut_alone / n_resp

      ex_hist_prob <- ex_hist_dens |>
        # dplyr::group_by(bin) |>
        # dplyr::slice(1) |>
        dplyr::mutate(prob_pos = 1 - (1 - prior_prob) * dens_all / dens) |>
        dplyr::mutate(
          prob_pos = pmax(0, prob_pos),
          prob_pos = pmin(1, prob_pos)
        )

      ex_hist_prob_lab_tbl <- ex_hist_prob |>
        dplyr::group_by(bin) |>
        dplyr::slice(1) |>
        dplyr::ungroup()

      prob_pos_lab_vec <- stats::setNames(
        ex_hist_prob_lab_tbl$prob_pos,
        ex_hist_prob_lab_tbl$bin_mid
      )

      ex_hist_weight <- ex_hist_weight |>
        dplyr::mutate(prob_pos = prob_pos_lab_vec[as.character(bin_mid)])

      lb_ind <- max(which(prob_pos_lab_vec < 0.5))
      ub_ind <- which(prob_pos_lab_vec >= 0.5)
      ub_ind <- ub_ind[min(which(ub_ind > lb_ind))]

      if (lb_ind == -Inf && ub_ind == 1) {
        cp <- as.numeric(names(prob_pos_lab_vec)[1]) - 1
      } else if (lb_ind == length(prob_pos_lab_vec) && ub_ind == Inf) {
        cp <- as.numeric(names(prob_pos_lab_vec)[length(prob_pos_lab_vec)]) + 1
      } else {
        lb_prob <- prob_pos_lab_vec[lb_ind]
        ub_prob <- prob_pos_lab_vec[ub_ind]
        lb_expr <- names(prob_pos_lab_vec)[lb_ind] |> as.numeric()
        ub_expr <- names(prob_pos_lab_vec)[ub_ind] |> as.numeric()
        increase_from_lb_expr <- (0.5 - lb_prob) /
          (ub_prob - lb_prob) *
          (ub_expr - lb_expr)
        cp <- lb_expr + increase_from_lb_expr
      }
      cp_prob <- cp

      # get cp based on mid-point
      if (length(expr_vec_resp) >= 2) {
        density_exc_min <- density(expr_vec_all)
        dens_tbl <- tibble::tibble(
          x = density_exc_min$x,
          y = density_exc_min$y
        )
        peak_all <- dens_tbl |>
          dplyr::filter(y == max(y)) |>
          dplyr::pull("x")


        density_exc_min <- density(expr_vec_resp)
        dens_tbl <- tibble::tibble(
          x = density_exc_min$x,
          y = density_exc_min$y
        )
        peak_resp <- dens_tbl |>
          dplyr::filter(x > peak_all + 30) |>
          dplyr::filter(y == max(y)) |>
          dplyr::pull("x")

        min_pt_tbl <- dens_tbl |>
          dplyr::filter(x >= peak_all & x <= peak_resp) |>
          dplyr::filter(y == min(y)) |>
          dplyr::slice(1)
        if (nrow(min_pt_tbl) == 0) {
          cp_min_dens <- cp_prob
        } else {
          cp_min_dens <- min_pt_tbl[["x"]][1]
        }
      } else {
        cp_min_dens <- cp_prob
      }


      cp <- max(cp_prob, cp_min_dens)

      ex_hist_weight <- ex_hist_weight |> dplyr::arrange(desc(expr))
      n_pos <- sum(ex_hist_weight$prob_pos * ex_hist_weight$prob_resp)
      k <- 0
      n_pos_tally <- 0

      ggplot(
        ex_plot |>
          dplyr::filter(chnl == cut_curr) |>
          dplyr::filter(
            expr > min_chnl_vec[cut_curr],
            !is.na(prob_resp) | type == "all"
          ),
        aes(x = expr, fill = type)
      ) +
        geom_density(alpha = 0.3) +
        scale_fill_brewer(palette = "Set1") +
        geom_vline(
          data = gate_tbl_plot |>
            dplyr::filter(cut == cut_curr),
          aes(
            xintercept = gate,
            linetype = gate_name
          ), size = 2
        ) +
        geom_vline(xintercept = cp, col = "red", size = 2) +
        facet_wrap(~marker, scales = "free") +
        geom_rug(
          data = ex_dens_calc |>
            dplyr::filter(chnl == cut_curr),
          aes(x = expr)
        )


      p_dens_weight_resp_vs_all <- ggplot(
        ex_hist_plot |>
          dplyr::mutate(
            dens = ifelse(
              type == "all",
              dens * (1 - prior_prob),
              dens
            )
          ),
        aes(x = expr, y = dens, col = type)
      ) +
        geom_line() +
        geom_vline(
          data = gate_tbl_plot |>
            dplyr::filter(cut == cut_curr),
          aes(
            xintercept = gate,
            linetype = gate_name
          ), size = 2
        ) +
        geom_vline(xintercept = cp, col = "red")

      cp
    }) |>
      stats::setNames(cut)

    gate_tbl_plot <- gate_tbl_plot |>
      dplyr::select(cut, marker, gate_name, gate)

    gate_tbl_plot <- gate_tbl_plot |>
      dplyr::bind_rows(tibble::tibble(
        cut = names(cp_vec),
        marker = chnl_lab_vec[names(cp_vec)],
        gate = cp_vec,
        gate_name = "adjusted"
      ))

    ggplot(ex_plot, aes(x = expr, fill = type)) +
      cowplot::theme_cowplot(font_size = 36) +
      geom_density(alpha = 0.7) +
      scale_fill_manual(
        values = c(
          "all" = "dodgerblue",
          "resp" = "orange"
        ),
        labels = c(
          "all" = "All cells",
          "resp" = "Responders only"
        ),
        name = ""
      ) +
      scale_linetype_manual(
        labels = c(
          "adjusted" = "Adjusted",
          "locb5_no" = "Local FDR",
          "locb5_min" = "Local FDR",
          "man_no" = "Manual"
        ),
        values = c(
          "adjusted" = "21",
          "locb5_no" = "solid",
          "man_no" = "11",
          "locb5_min" = "solid"
        ),
        name = ""
      ) +
      geom_vline(
        data = gate_tbl_plot,
        aes(
          xintercept = gate,
          linetype = gate_name
        ), size = 3
      ) +
      facet_wrap(~marker, scales = "free") +
      geom_rug(data = ex_plot |>
        dplyr::filter(type == "resp")) +
      theme(legend.key.height = unit(3, "cm")) +
      labs(x = "Marker level", y = "Density")

    .add_pos_col <- function(data, cut, cp) {
      for (cut_curr in cut) {
        data[[paste0(cut_curr, "_pos")]] <-
          as.numeric(data[[cut_curr]] > cp[cut_curr])
      }
      data
    }

    ex_combn <- ex_resp |>
      dplyr::filter(!is.na(prob_resp)) |>
      dplyr::select(-marker) |>
      tidyr::pivot_wider(names_from = chnl, values_from = expr)

    ex_combn <- .add_pos_col(data = ex_combn, cut = cut, cp = cp_vec)

    ex_combn <- ex_combn |>
      dplyr::mutate(
        combn = paste0(
          Ho165Di_pos, "_", Gd158Di_pos, "_", Nd146Di_pos
        )
      ) |>
      dplyr::group_by(combn) |>
      dplyr::summarise(count = sum(prob_resp))

    ex_combn_2 <- ex_resp |>
      dplyr::filter(!is.na(prob_resp)) |>
      dplyr::select(-marker) |>
      tidyr::pivot_wider(names_from = chnl, values_from = expr)

    gate_tbl_man_ind <- gate_tbl_man |>
      dplyr::filter(.data$ind == .env$ind)
    cp_vec_man <- stats::setNames(gate_tbl_man_ind$cp, gate_tbl_man_ind$cut)

    gate_tbl_auto_ind <- gate_tbl |>
      dplyr::filter(.data$ind == .env$ind)
    cp_vec_auto <- stats::setNames(gate_tbl_auto_ind$cp, gate_tbl_auto_ind$cut)

    ex_combn_2 <- .add_pos_col(data = ex_combn_2, cut = cut, cp = cp_vec_auto)

    exc_combn_2 <- ex_combn_2 |>
      dplyr::mutate(
        combn = paste0(Ho165Di_pos, "_", Gd158Di_pos, "_", Nd146Di_pos)
      ) |>
      dplyr::group_by(combn) |>
      dplyr::summarise(count = sum(prob_resp))

    combn_tbl <- tibble::tibble(combn = .get_all_combn(length(cut))) |>
      dplyr::mutate(
        ind = ind, cut = paste0(cut, collapse = "_"),
        pop = pop_gate, gate_name = gate_name
      )

    .add_pos_col_2 <- function(data, cut, chnl_lab) {
      for (i in seq_along(cut)) {
        cut_curr <- cut[i]
        if (i == 1) {
          ind <- 1
        } else {
          ind <- 1 + (i - 1) * 2
        }
        data[[cut_curr]] <- stringr::str_sub(data$combn, ind, ind)
        if (i == 1) {
          data <- data |> dplyr::mutate(cyt_combn = "")
        }
        marker <- chnl_lab[cut_curr] |>
          stringr::str_remove("-beads")
        cyt_lab_vec <- c(
          "0" = paste0(marker, "-"),
          "1" = paste0(marker, "+")
        )
        data <- data |>
          dplyr::mutate(
            cyt_combn = paste0(
              cyt_combn, cyt_lab_vec[
                data[[as.character(cut_curr)]]
              ]
            )
          )
      }
      data
    }

    combn_tbl <- .add_pos_col_2(
      data = combn_tbl, cut = cut, chnl_lab = chnl_lab_vec
    )

    ex_combn_diff <- combn_tbl |>
      dplyr::left_join(ex_combn, by = "combn") |>
      dplyr::mutate(count = ifelse(is.na(count), 0, count)) |>
      dplyr::left_join(ex_combn_2 |>
        dplyr::rename(count_man = count), by = "combn") |>
      dplyr::mutate(count_man = ifelse(is.na(count_man), 0, count_man))

    all_neg_obs <- paste0(rep(0, length(cut)), collapse = "_")
    ex_combn_diff <- exc_combn_diff |>
      dplyr::mutate(diff = ifelse(combn == all_neg_obs,
        0, count - count_man
      )) |>
      dplyr::mutate(
        diff_sum = sum(diff),
        diff_sum_abs = sum(abs(count - count_man))
      ) |>
      dplyr::mutate(
        freq_bs = count / n_cell * 1e2,
        freq_bs_man = count_man / n_cell * 1e2,
        n_cell = n_cell
      )

    ex_combn_diff
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()
  params[["cut"]] <- cut_curr
  dir_base <- stimgate_dir_base_create(
    dir_base_init = path_project,
    params = params
  )
  # get stats tbl
  gate_stats <- readRDS(file.path(dir_base, "stats", "gate_stats_tbl"))


  gate_tbl <- purrr::map_df(cut, function(cut_curr) {
    # get base directory
    params[["cut"]] <- cut_curr
    dir_base <- stimgate_dir_base_create(
      dir_base_init = path_project,
      params = params
    )
    # get stats tbl
    gate_stats <- readRDS(file.path(dir_base, "stats", "gate_stats_tbl"))

    gate_stats |>
      dplyr::filter(.data$gate_name == .env$gate_name) |>
      dplyr::mutate(cut = cut_curr, marker = chnl_lab_vec[cut_curr]) |>
      dplyr::select(cut, marker, everything())
  })
}

.get_all_combn <- function(n_elem) {
  elem_list <- list(c("0", "1"))
  if (n_elem == 1) {
    return(elem_list |> unlist())
  }
  n_elem <- 3
  for (i in 2:n_elem) {
    elem_list <- purrr::map(elem_list, function(x) {
      list(paste0(x, "_0"), paste0(x, "_1")) |>
        unlist()
    })
  }
  elem_list |>
    unlist()
}
