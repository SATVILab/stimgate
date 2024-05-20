.get_inc_vec <- function(chnl_curr,
                         chnl_vec,
                         ex,
                         gate_tbl_ind) {
  inc_vec <- rep(FALSE, nrow(ex))

  for (chnl_alt in setdiff(chnl_vec, chnl_curr)) {
    cp <- gate_tbl_ind |>
      dplyr::filter(.data$chnl == chnl_alt) |> # nolint
      dplyr::pull("gate")
    inc_vec <- inc_vec | ex[[chnl_alt]] > cp
  }

  inc_vec
}

.get_cp_neg <- function(ex,
                        inc,
                        chnl,
                        bw_min,
                        min_cell = 1e3,
                        max_peak_ratio = 1e3) {
  # get cytokine-negative cells
  ex_neg <- ex[!inc & ex[[chnl]] > min(ex[[chnl]]), ][[chnl]]

  # if too few, then return
  if (length(ex_neg) <= min_cell) {
    return(NA)
  }

  # calculate density
  dens_neg <- density(ex_neg, bw = "SJ")
  if (dens_neg$bw < bw_min) dens_neg <- density(ex_neg, bw = bw_min)

  # calculate mode furthest to the right
  dens_len <- length(dens_neg$x)

  mode_ind_vec <- (2:(dens_len - 1))[
    (dens_neg$y[-1] > dens_neg$y[-dens_len])[-(dens_len - 1)] &
      (dens_neg$y[-dens_len] > dens_neg$y[-1])[-1]
  ]
  mode_vec <- dens_neg$x[mode_ind_vec]
  mode_vec_y <- dens_neg$y[mode_ind_vec]
  mode_ind_vec <- mode_ind_vec[mode_vec_y > 0.01 * max(mode_vec_y)]
  mode_vec <- dens_neg$x[mode_ind_vec]

  # get all points that are to the
  # right of the right-most mode that
  # also are much smaller than peak
  low_neg_dens_pts_vec <- dens_neg$x[
    dens_neg$y < max(dens_neg$y) / max_peak_ratio &
      dens_neg$x >= max(mode_vec)
  ]

  # return NA if no such no low
  # neg dens points found, other return left-most such point
  ifelse(length(low_neg_dens_pts_vec) == 0, NA, min(low_neg_dens_pts_vec))
}


.get_cp_pos_ind <- function(ex,
                            inc,
                            chnl,
                            bw_min,
                            adjust = 1,
                            trust_no_or_high_am = FALSE,
                            min_cell = 10,
                            cp_orig) {
  # data
  ex_pos <- ex[inc & ex[[chnl]] > min(ex[[chnl]]), ][[chnl]]
  if (length(ex_pos) < min_cell) {
    return(NA)
  }
  ex_neg <- ex[!inc & ex[[chnl]] > min(ex[[chnl]]), ][[chnl]]

  # bw
  bw_dens <- density(ex_pos, bw = "SJ")$bw
  bw_sd <- sd(ex_neg)
  bw_final <- max(bw_dens, bw_sd, bw_min)

  # calculate density
  # ------------------

  dens_pos <- density(ex_pos, bw = bw_final, adjust = adjust)
  dens_neg <- density(ex_neg,
    bw = bw_final, adjust = adjust,
    from = min(dens_pos$x), to = max(dens_pos$x)
  )

  # calculate modes and antimodes
  # -----------------------------

  dens_len <- length(dens_pos$y)

  # antimodes
  am_ind_vec <- (2:(dens_len - 1))[
    (dens_pos$y[-1] < dens_pos$y[-dens_len])[-(dens_len - 1)] &
      (dens_pos$y[-dens_len] < dens_pos$y[-1])[-1]
  ]
  am_vec <- dens_pos$x[am_ind_vec]
  am_vec_height <- dens_pos$y[am_ind_vec]

  # modes
  mode_ind_vec <- (2:(dens_len - 1))[
    (dens_pos$y[-1] > dens_pos$y[-dens_len])[-(dens_len - 1)] &
      (dens_pos$y[-dens_len] > dens_pos$y[-1])[-1]
  ]
  mode_vec <- dens_pos$x[mode_ind_vec]
  mode_vec_height <- dens_pos$y[mode_ind_vec]
  mode_ind_vec <- mode_ind_vec[mode_vec_height > 0.01 * max(mode_vec_height)]
  mode_vec <- dens_pos$x[mode_ind_vec]
  mode_vec_height <- dens_pos$y[mode_ind_vec]

  # calculate lowest mode for cyt-neg cells
  # --------------------------
  mode_ind_vec_neg <- (2:(dens_len - 1))[
    (dens_neg$y[-1] > dens_neg$y[-dens_len])[-(dens_len - 1)] &
      (dens_neg$y[-dens_len] > dens_neg$y[-1])[-1]
  ]
  mode_vec_height_neg <- dens_neg$y[mode_ind_vec_neg]
  mode_ind_vec_neg <- mode_ind_vec_neg[
    mode_vec_height_neg > 0.01 * max(mode_vec_height_neg)
  ]
  highest_mode_neg <- max(dens_neg$x[mode_ind_vec_neg])

  # calculate cp_shape
  # -------------------

  # if no antimode
  if (length(am_vec) == 0) {
    if (!trust_no_or_high_am) {
      return(NA)
    }
    if (trust_no_or_high_am) {
      cp_shape <- ifelse(cp_orig < max(mode_vec), highest_mode_neg, NA)
      return(cp_shape)
    }
  }

  # nearest mode to cp_orig
  mode_vec_above_cp_orig <- mode_vec[mode_vec > cp_orig]

  # return NA there are none
  if (length(mode_vec_above_cp_orig) == 0) {
    return(NA)
  }
  mode_above_cp_orig_min <- min(mode_vec_above_cp_orig)

  # now there is one mode above cp_orig
  am_vec_more_than_cp_orig <- am_vec[am_vec > cp_orig]
  am_right_min_ind <- ifelse(
    length(am_vec_more_than_cp_orig) > 0,
    which(am_vec == min(am_vec_more_than_cp_orig)), NA
  )
  am_right_min <- am_vec[am_right_min_ind]
  # return am between cp_orig and mode to the right of cp_orig as
  # cp_shape if there is one
  if (!is.na(am_right_min[1])) {
    if (am_right_min < mode_above_cp_orig_min) {
      return(am_right_min)
    }
  }

  # now we know that there is no antimode between nearest mode to
  # cp_orig on right and cp_orig
  # so now we want to
  # if the am is deep enough, then use it to return it
  # get peak on left and right of am

  # get left-most antimode of antimodes less than cp_orig
  am_vec_less_than_cp_orig <- am_vec[am_vec < cp_orig]
  max_left_am_ind <- ifelse(
    length(am_vec_less_than_cp_orig) > 0,
    which(am_vec == max(am_vec_less_than_cp_orig)), NA
  )
  max_left_am <- am_vec[max_left_am_ind]
  max_left_am_height <- am_vec_height[max_left_am_ind]

  # if there are none, then if
  # trust no am then return peak_neg_dens_x, else return NA
  if (length(max_left_am) == 0 || all(is.na(max_left_am))) {
    cp_shape <- ifelse(trust_no_or_high_am, highest_mode_neg, NA)
    return(cp_shape)
  }

  # if there is one, but it's too small
  # (relative to mode immediately to left of it)
  # and it's not trust_no_or_high_am, then return NA
  max_mode_less_than_cp_orig <- max(mode_vec[mode_vec < cp_orig])
  max_mode_less_than_cp_orig_height <- mode_vec_height[
    which(mode_vec == max_mode_less_than_cp_orig)
  ]
  if (max_left_am_height > 0.5 * max_mode_less_than_cp_orig_height) {
    if (trust_no_or_high_am) {
      min_mode_above_cp_orig_ind <- which(
        mode_vec == min(mode_vec[mode_vec > cp_orig])
      )
      min_mode_above_cp_orig_height <- mode_vec_height[
        min_mode_above_cp_orig_ind
      ]
      min_mode_above_cp_orig <- mode_vec[min_mode_above_cp_orig_ind]

      prop_move_to_right <- max_mode_less_than_cp_orig_height /
        (max_mode_less_than_cp_orig_height + min_mode_above_cp_orig_height)
      cp_shape <- max(
        max_left_am,
        max_mode_less_than_cp_orig +
          prop_move_to_right *
            (min_mode_above_cp_orig - max_mode_less_than_cp_orig)
      )
      return(cp_shape)
    } else {
      return(NA)
    }
  }

  # return max_left_am
  # this only happens if there is a left-most antimode and it is deep enough
  max_left_am
}

.get_cp_pos <- function(ex,
                        inc,
                        chnl,
                        bw_min,
                        trust_no_or_high_am = FALSE,
                        min_cell = 10,
                        cp_orig, n_loop = 5) {
  cp_pos <- .get_cp_pos_ind(
    ex = ex, inc = inc, chnl = chnl, bw_min = bw_min, adjust = 1,
    trust_no_or_high_am = FALSE, min_cell = 10,
    cp_orig = cp_orig
  )

  k <- 1
  while (is.na(cp_pos) && k <= n_loop) {
    if (is.na(cp_pos)) {
      cp_pos <- .get_cp_pos_ind(
        ex = ex, inc = inc, chnl = chnl, bw_min = bw_min, adjust = 0.5^k,
        trust_no_or_high_am = FALSE, min_cell = 10,
        cp_orig = cp_orig
      )
    }
    k <- k + 1
  }

  cp_pos
}


.get_cyt_pos_gates_chnl_vec_from_marker_list <- function(marker_list) {
  purrr::map_chr(marker_list, function(x) x$cut)
}

.get_cyt_pos_gates_gate_tbl_get <- function(chnl_vec,
                                            path_project,
                                            params,
                                            gate_name,
                                            debug,
                                            chnl_lab_vec) {
  .debug(debug, "Getting gate_tbl") # nolint
  purrr::map_df(chnl_vec, function(chnl_curr) {
    # get base directory
    dir_base <- stim_gate_dir_base_create( # nolint
      dir_base_init = path_project,
      params = params |> append(list(cut = chnl_curr))
    )
    # get stats tbl
    gate_tbl <- readRDS(file.path(dir_base, "gate_tbl_init.rds"))

    if (!is.null(gate_name)) {
      gate_tbl <- gate_tbl |>
        dplyr::filter(gate_name == .env$gate_name) # nolint
    }

    gate_tbl |>
      # dplyr::filter(.data$gate_name == .env$gate_name) |>
      dplyr::mutate(chnl = chnl_curr, marker = chnl_lab_vec[chnl_curr]) |>
      dplyr::select(chnl, marker, gate_name, batch, ind, gate) # nolint
  })
}
