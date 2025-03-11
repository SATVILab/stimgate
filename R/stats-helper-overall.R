.get_stats_overall <- function(ind_batch_list,
                               .data,
                               pop_gate,
                               gate_tbl,
                               gate_name,
                               chnl,
                               chnl_cut,
                               chnl_lab,
                               filter_other_cyt_pos,
                               combn,
                               gate_type_cyt_pos_filter,
                               gate_type_single_pos_filter,
                               gate_type_single_pos_calc,
                               gate_type_cyt_pos_calc,
                               combn_mat_list, # missing
                               cyt_combn_vec_list, # missing
                               debug) {
  stat_tbl <- purrr::map_df(
    seq_along(ind_batch_list),
    function(i) {
      .get_stats_overall_progress(
        ind_batch_list = ind_batch_list,
        i = i,
        debug = debug,
        combn = combn,
        filter_other_cyt_pos = filter_other_cyt_pos
      )
      .get_stats_batch(
        ind_batch = ind_batch_list[[i]],
        batch = names(ind_batch_list)[i],
        debug = debug,
        .data = .data,
        pop_gate = pop_gate,
        gate_tbl = gate_tbl,
        chnl = chnl,
        filter_other_cyt_pos = filter_other_cyt_pos,
        combn = combn,
        gate_type_cyt_pos_filter = gate_type_cyt_pos_filter,
        gate_type_single_pos_filter = gate_type_single_pos_filter,
        gate_type_single_pos_calc = gate_type_single_pos_calc,
        gate_type_cyt_pos_calc = gate_type_cyt_pos_calc,
        combn_mat_list = combn_mat_list,
        cyt_combn_vec_list = cyt_combn_vec_list,
        gate_name = gate_name
      )
    }
  )

  stat_tbl <- stat_tbl |>
    dplyr::mutate(
      prop_stim = count_stim / n_cell_stim, # nolint
      prop_uns = count_uns / n_cell_uns, # nolint
      prop_bs = prop_stim - prop_uns, # nolint
      freq_stim = prop_stim * 1e2, # nolint
      freq_uns = prop_uns * 1e2, # nolint
      freq_bs = freq_stim - freq_uns # nolint
    )

  stat_tbl <- .get_stats_update_combn_n( # nolint
    combn = combn,
    stat_tbl = stat_tbl,
    chnl_cut,
    chnl_lab = chnl_lab
  )

  if ("ind" %in% colnames(stat_tbl)) {
    stat_tbl[, "ind"] <- as.character(stat_tbl[["ind"]])
  }

  # add SampleID and update stim labels
  # if desired
  .get_stats_label( # nolint
    stat_tbl = stat_tbl
  )
}

.get_stats_overall_progress <- function(ind_batch_list,
                                        i,
                                        debug,
                                        combn,
                                        filter_other_cyt_pos) {
  ind_batch <- ind_batch_list[[i]]
  .debug( # nolint
    debug, "ind_batch: ", paste0(ind_batch, collapse = "-")
  )
  # print progress
  if (i %% 10 == 0 || i == length(ind_batch_list)) {
    if (combn && !filter_other_cyt_pos) {
      print(paste0("batch ", i, " of ", length(ind_batch_list)))
    }
  }
}


.get_stats_batch <- function(ind_batch,
                             batch,
                             debug,
                             .data,
                             pop_gate,
                             gate_tbl,
                             chnl,
                             filter_other_cyt_pos,
                             combn,
                             gate_type_cyt_pos_filter,
                             gate_type_single_pos_filter,
                             gate_type_single_pos_calc,
                             gate_type_cyt_pos_calc,
                             combn_mat_list,
                             cyt_combn_vec_list,
                             gate_name) {
  # calculate n_cell_[stim/uns] and count_[stim/uns]
  # for each gate type (gn) for either:
  # individual cytokines (combn = FALSE) or
  # combinations of cytokines (combn = TRUE), or
  # when using only cells positive for other cytokines
  # (filter_other_cyt_pos = TRUE)
  .debug(debug, "Getting gate stats for a batch") # nolint
  .debug(debug, "ind_batch: ", paste0(ind_batch, collapse = "-")) # nolint

  # read in .data
  ex_list <- .get_ex_list( # nolint
    .data = .data,
    ind_batch = ind_batch,
    batch = batch,
    pop = pop_gate,
    chnl_cut = unique(gate_tbl$chnl)
  )

  purrr::map_df(gate_name, function(gn) {
    .get_stats_batch_gn(
      gn = gn,
      debug = debug,
      ex_list = ex_list,
      gate_tbl = gate_tbl,
      chnl = chnl,
      filter_other_cyt_pos = filter_other_cyt_pos,
      gate_type_cyt_pos_filter = gate_type_cyt_pos_filter,
      gate_type_single_pos_filter = gate_type_single_pos_filter,
      gate_type_single_pos_calc = gate_type_single_pos_calc,
      gate_type_cyt_pos_calc = gate_type_cyt_pos_calc,
      combn = combn,
      combn_mat_list = combn_mat_list,
      cyt_combn_vec_list = cyt_combn_vec_list,
      ind_batch = ind_batch
    )
  })
}

.get_stats_batch_gn <- function(gn,
                                debug,
                                ex_list,
                                gate_tbl,
                                chnl,
                                filter_other_cyt_pos,
                                gate_type_cyt_pos_filter,
                                gate_type_single_pos_filter,
                                gate_type_single_pos_calc,
                                gate_type_cyt_pos_calc,
                                combn,
                                combn_mat_list,
                                cyt_combn_vec_list,
                                ind_batch) {
  .debug(debug, "gate name: ", gn) # nolint
  gate_tbl_gn <- gate_tbl |> dplyr::filter(gate_name == gn) # nolint
  if (filter_other_cyt_pos || !combn) {
    # get stats when not calculating cytokine-combinations
    # or when filtering to yield cells positive
    # for all other cytokines
    stat_tbl_gn <- .get_stats_batch_gn_filter_or_non_combn(
      debug = debug,
      ex_list = ex_list,
      ind_batch = ind_batch,
      gate_tbl_gn = gate_tbl_gn,
      gn = gn,
      chnl = chnl,
      filter_other_cyt_pos = filter_other_cyt_pos,
      gate_type_single_pos_calc = gate_type_single_pos_calc,
      gate_type_cyt_pos_filter = gate_type_cyt_pos_filter,
      gate_type_single_pos_filter = gate_type_single_pos_filter
    )
    return(stat_tbl_gn)
  }
  # get stats when calculating cytokine-combinations
  # note: it doesn't make sense to have filter_other_cyt_pos = TRUE
  # and combn = TRUE, as the cytokine combinations
  # calculation already involves the other cytokines, so also
  # filtering for other cytokines doesn't make sense.
  .get_stats_batch_gn_combn_loop_ind(
    ex_list = ex_list,
    gate_tbl_gn = gate_tbl_gn,
    gn = gn,
    chnl = chnl,
    combn_mat_list = combn_mat_list,
    cyt_combn_vec_list = cyt_combn_vec_list,
    gate_type_cyt_pos_calc = gate_type_cyt_pos_calc,
    gate_type_single_pos_calc = gate_type_single_pos_calc,
    debug = debug
  )
}

.get_stats_batch_gn_combn_loop_ind <- function(ex_list,
                                               gate_tbl_gn,
                                               gn,
                                               chnl,
                                               combn_mat_list,
                                               cyt_combn_vec_list,
                                               gate_type_cyt_pos_calc,
                                               gate_type_single_pos_calc,
                                               debug) {
  # filter to yield cells negative for all cytokine combinations
  ex_list_stim <- ex_list[-length(ex_list)]
  ex_uns <- ex_list[[length(ex_list)]]
  n_cell_uns <- nrow(ex_uns) # nolint
  purrr::map_df(seq_along(ex_list_stim), function(i) {
    .debug(debug, "i: ", i) # nolint
    ex <- ex_list_stim[[i]]
    gate_tbl_gn_ind <- gate_tbl_gn |>
      dplyr::filter(ind == attr(ex, "ind")) # nolint
    combn_tbl <- purrr::map_df(names(combn_mat_list), function(j) {
      .get_stats_batch_gn_combn(
        j = j,
        debug = debug,
        ex = ex,
        ex_uns = ex_uns,
        gate_tbl_gn_ind = gate_tbl_gn_ind,
        gn = gn,
        chnl = chnl,
        combn_mat_list = combn_mat_list,
        cyt_combn_vec_list = cyt_combn_vec_list,
        gate_type_cyt_pos_calc = gate_type_cyt_pos_calc,
        gate_type_single_pos_calc = gate_type_single_pos_calc
      )
    }) |>
      dplyr::mutate(
        n_cell_stim = nrow(ex),
        n_cell_uns = .env$n_cell_uns # nolint
      )
    combn_tbl |> .get_stats_batch_gn_combn_neg(chnl)
  })
}

.get_stats_batch_gn_combn <- function(j,
                                           debug,
                                           ex,
                                           ex_uns,
                                           gate_tbl_gn_ind,
                                           gn,
                                           chnl,
                                           combn_mat_list,
                                           cyt_combn_vec_list,
                                           gate_type_cyt_pos_calc,
                                           gate_type_single_pos_calc) {
  .debug(debug, "number of cytokines positive: ", j) # nolint
  combn_mat <- combn_mat_list[[j]]
  cyt_combn <- cyt_combn_vec_list[[j]]
  stat_tbl_gn_ind <- tibble::tibble(
    ind = attr(ex, "ind"),
    gate_name = gn,
    cyt_combn = cyt_combn,
    count_stim = NA_integer_,
    n_cell_stim = NA_integer_,
    count_uns = NA_integer_,
    n_cell_uns = NA_integer_
  )

  for (i in seq_len(nrow(stat_tbl_gn_ind))) {
    .debug(debug, "i: ", i) # nolint
    chnl_pos <- chnl[combn_mat[i, , drop = TRUE]]
    chnl_neg <- chnl[
      setdiff(seq_along(chnl), combn_mat[i, , drop = TRUE])
    ]
    stat_tbl_gn_ind[i, "count_stim"] <- sum(
      .get_pos_ind_cyt_combn( # nolint
        ex = ex, gate_tbl = gate_tbl_gn_ind,
        chnl_pos = chnl_pos, chnl_neg = chnl_neg,
        chnl_alt = NULL,
        gate_type_cyt_pos = gate_type_cyt_pos_calc,
        gate_type_single_pos = gate_type_single_pos_calc
      )
    )
    stat_tbl_gn_ind[i, "count_uns"] <- sum(
      .get_pos_ind_cyt_combn( # nolint
        ex = ex_uns, gate_tbl = gate_tbl_gn_ind |>
          dplyr::mutate(ind = attr(ex_uns, "ind")),
        chnl_pos = chnl_pos, chnl_neg = chnl_neg,
        chnl_alt = NULL,
        gate_type_cyt_pos = gate_type_cyt_pos_calc,
        gate_type_single_pos = gate_type_single_pos_calc
      )
    )
  }
  stat_tbl_gn_ind
}

.get_stats_batch_gn_combn_neg <- function(.data, chnl) {
  all_neg_row <- .data |>
    dplyr::mutate(cyt_combn = paste0(paste0(chnl, collapse = "~-~"), "~-~")) |>
    dplyr::group_by(ind, cyt_combn, gate_name) |>
    dplyr::summarise(
      count_stim = n_cell_stim[[1]] - sum(count_stim), n_cell_stim = n_cell_stim[[1]],
      count_uns = n_cell_uns[[1]] - sum(count_uns), n_cell_uns = n_cell_uns[[1]],
      .groups = "drop"
    )
  .data |> dplyr::bind_rows(all_neg_row)
}

.get_stats_batch_gn_filter_or_non_combn <- function(debug,
                                                         ex_list,
                                                         ind_batch,
                                                         gate_tbl_gn,
                                                         gn,
                                                         chnl,
                                                         filter_other_cyt_pos,
                                                         gate_type_single_pos_calc, # nolint
                                                         gate_type_cyt_pos_filter, # nolint
                                                         gate_type_single_pos_filter) {
  .debug(debug, "filtering or not working out combinations") # nolint
  purrr::map_df(chnl, function(chnl_curr) {
    .debug(debug, "chnl_curr: ", chnl_curr) # nolint

    # remove cells positive for other cytokines
    if (filter_other_cyt_pos) {
      ex_list <-
        .get_stats_batch_gn_filter_or_non_combn_filter(
          ex_list = ex_list,
          gate_tbl_gn = gate_tbl_gn,
          chnl_curr = chnl_curr,
          gate_type_cyt_pos_filter = gate_type_cyt_pos_filter,
          gate_type_single_pos_filter = gate_type_single_pos_filter,
          debug = debug
        )
    }

    # get statistics for the individuals
    stat_tbl_gn_ind <- tibble::tibble(
      ind = ind_batch[-length(ind_batch)],
      gate_name = gn,
      chnl = chnl_curr,
      count_stim = NA,
      n_cell_stim = NA,
      count_uns = NA,
      n_cell_uns = NA
    )
    for (j in seq_len(nrow(stat_tbl_gn_ind))) {
      .debug(debug, "j: ", j) # nolint
      ex <- ex_list[[j]]
      gate_tbl_gn_ind <- gate_tbl_gn |>
        dplyr::filter(ind == attr(ex, "ind")) # nolint
      nothing_to_gate <- nrow(ex) == 0 ||
        nrow(gate_tbl_gn_ind) == 0 ||
        all(is.na(ex[[chnl_curr]]))
      if (nothing_to_gate) {
        .debug(debug, "filling in NAs") # nolint
        stat_tbl_gn_ind[j, "count_stim"] <- NA_integer_
        stat_tbl_gn_ind[j, "n_cell_stim"] <- min(
          sum((!is.na(ex[[chnl_curr]])) & (!is.nan(ex[[chnl_curr]])))
        )
        stat_tbl_gn_ind[j, "count_uns"] <- NA_integer_
        stat_tbl_gn_ind[j, "n_cell_uns"] <-
          nrow(ex_list[[length(ex_list)]])
        next
      }
      cn_vec <- colnames(gate_tbl_gn_ind)
      gate_col_ind <- switch(gate_type_single_pos_calc,
        "base" = which(cn_vec == "gate"),
        "single" = which(cn_vec == "gate_single"),
        stop(paste0(
          "gate_type_single_pos_calc value of ",
          gate_type_single_pos_calc,
          " not either 'base' or 'single'."
        ))
      )
      gate_gn_ind_chnl <-
        gate_tbl_gn_ind[[gate_col_ind]][gate_tbl_gn_ind$chnl == chnl_curr]
      stat_tbl_gn_ind[j, "count_stim"] <-
        sum(ex[[chnl_curr]] > gate_gn_ind_chnl)
      stat_tbl_gn_ind[j, "n_cell_stim"] <-
        nrow(ex)
      ex_uns <- ex_list[[length(ex_list)]]
      if (filter_other_cyt_pos) {
        pos_ind_vec_but_single_pos_curr <-
          .get_pos_ind_but_single_pos_for_one_cyt( # nolint
            ex = ex_uns |>
              dplyr::mutate(is_uns = FALSE),
            gate_tbl = gate_tbl_gn_ind |>
              dplyr::mutate(ind = attr(ex, "ind")),
            chnl_single_exc = chnl_curr,
            chnl = NULL,
            gate_type_cyt_pos = gate_type_cyt_pos_filter,
            gate_type_single_pos = gate_type_single_pos_filter
          )
        ex_uns <- ex_uns[!pos_ind_vec_but_single_pos_curr, , drop = FALSE]
      }
      stat_tbl_gn_ind[j, "count_uns"] <-
        sum(ex_uns[[chnl_curr]] > gate_gn_ind_chnl)
      stat_tbl_gn_ind[j, "n_cell_uns"] <-
        nrow(ex_uns)
    }
    stat_tbl_gn_ind
  })
}


.get_stats_batch_gn_filter_or_non_combn_filter <- function(ex_list,
                                                           gate_tbl_gn,
                                                           chnl_curr,
                                                           gate_type_cyt_pos_filter, # nolint
                                                           gate_type_single_pos_filter, # nolint
                                                           debug) {
  .debug(debug, "filtering other cyt pos") # nolint

  # loop over individual samples
  purrr::map(seq_along(ex_list), function(i) {
    .debug(debug, "i: ", i) # nolint

    # return early if nothing to filter
    return_early <-
      .get_stats_batch_gn_filter_or_non_combn_filter_check_early(
        i = i,
        ex_list = ex_list,
        chnl_curr = chnl_curr,
        gate_tbl_gn = gate_tbl_gn
      )
    if (return_early) {
      return(ex_list[[i]])
    }

    # identify cells positive for other cytokine(s)
    pos_ind_vec_but_single_pos_curr <-
      .get_pos_ind_but_single_pos_for_one_cyt( # nolint
        ex = ex_list[[i]],
        gate_tbl = gate_tbl_gn |>
          dplyr::filter(ind == attr(ex_list[[i]], "ind")), # nolint
        chnl_single_exc = chnl_curr,
        chnl = NULL,
        gate_type_cyt_pos = gate_type_cyt_pos_filter,
        gate_type_single_pos = gate_type_single_pos_filter
      )
    # remove cells positive for other cytokine(s)
    ex_list[[i]][!pos_ind_vec_but_single_pos_curr, , drop = FALSE]
  }) |>
    stats::setNames(names(ex_list))
}

.get_stats_batch_gn_filter_or_non_combn_filter_check_early <- function(i,
                                                                       ex_list, # nolint
                                                                       chnl_curr, # nolint
                                                                       gate_tbl_gn) { # nolint
  if (i == length(ex_list)) {
    return(TRUE)
  }
  gate_tbl_gn_ind <- gate_tbl_gn |>
    dplyr::filter(ind == attr(ex_list[[i]], "ind")) # nolint

  # return early if nothing to filter
  all(is.na(ex_list[[i]][[chnl_curr]])) ||
    nrow(gate_tbl_gn_ind) == 0 ||
    nrow(ex_list[[i]]) == 0
}

.get_stats_update_combn_n <- function(combn,
                                      stat_tbl,
                                      chnl_cut,
                                      chnl_lab) {
  if (combn) {
    return(stat_tbl)
  }
  if ((!"chnl" %in% colnames(stat_tbl))) {
    stat_tbl <- stat_tbl |>
      dplyr::mutate(
        chnl = chnl_cut,
        marker = chnl_lab[chnl_cut]
      )
  }
  if ((!"marker" %in% colnames(stat_tbl))) {
    stat_tbl <- stat_tbl |>
      dplyr::mutate(marker = chnl_lab[.data$chnl]) # nolint
  }
  stat_tbl
}

.get_stats_label <- function(stat_tbl) {
  cn_vec_order <- c(
    "gate_name", "chnl", "marker", "ind"
  )
  cn_vec_order_curr <- cn_vec_order[cn_vec_order %in% colnames(stat_tbl)]
  stat_tbl |>
    dplyr::select(all_of(cn_vec_order_curr), everything()) # nolint
}
