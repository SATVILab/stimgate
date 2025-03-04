#' @title Get labels for the gates
#'
#' @inheritParams plot_list # gate_list
#'
#' @return A named character vector.
.get_gate_lab_vec <- function(gate_tbl) {
  gate_name_sh_vec <- unique(gate_tbl$gate_name) # .get_gate_name_vec_sh(gate_tbl = gate_tbl)

  purrr::map_chr(gate_name_sh_vec, function(gate_name_sh) {
    gate_name_sh_split <- stringr::str_split(gate_name_sh, "_")[[1]]
    if (any(stringr::str_detect(gate_name_sh_split, "ctrl"))) {
      ctrl_add <- "Ctrl"
      gate_name_sh_split <- gate_name_sh_split[-which(gate_name_sh_split == "ctrl")]
    } else {
      ctrl_add <- NULL
    }

    if (any(stringr::str_detect(gate_name_sh_split, "adj"))) {
      adj_add <- "Adj"
      gate_name_sh_split <- gate_name_sh_split[-which(gate_name_sh_split == "adj")]
    } else {
      adj_add <- NULL
    }

    if (any(stringr::str_detect(gate_name_sh_split, "clust"))) {
      clust_add <- "Clust"
      gate_name_sh_split <- gate_name_sh_split[-which(gate_name_sh_split == "clust")]
    } else {
      clust_add <- NULL
    }


    gate_type_sh <- gate_name_sh_split[1]
    gate_combn_sh <- gate_name_sh_split[2]

    gate_type_sh_lab_vec_auto <- c(
      "man" = "Manual", tg = "Finak (2014)",
      "midp" = "Mid probability", "scp" = "Single changepoint"
    )
    if (gate_type_sh %in% names(gate_type_sh_lab_vec_auto)) {
      gate_type_fh <- gate_type_sh_lab_vec_auto[[gate_type_sh]]
    } else if (stringr::str_detect(gate_type_sh, "uns")) {
      bias_loc <- stringr::str_locate(gate_type_sh, "b")[[1]]
      bias <- stringr::str_sub(gate_type_sh, start = bias_loc + 1)
      fdr <- stringr::str_sub(gate_type_sh, start = 4, end = bias_loc - 1)

      gate_type_fh <- paste0("FDR (fdr ", fdr, "; bias ", bias, ")")
    } else if (stringr::str_detect(gate_type_sh, "loc")) {
      bias_loc <- stringr::str_locate(gate_type_sh, "b")[[1]]
      bias <- stringr::str_sub(gate_type_sh, start = bias_loc + 1)
      gate_type_fh <- paste0("Local FDR (bias ", bias, ")")
    }


    gate_combn_sh_lab_vec <- c(
      "no" = "Ind", "prejoin" = "Batch",
      "min" = "Min", "max" = "Max",
      "med" = "Median", "mean" = "Mean",
      "trim20" = "Trim"
    )

    gate_combn_fh <- gate_combn_sh_lab_vec[[gate_combn_sh]]
    gate_combn_fh <- paste0(c(
      gate_combn_fh, adj_add, ctrl_add, clust_add
    ), collapse = " - ")

    paste0(gate_type_fh, " - ", gate_combn_fh)
  }) |>
    stats::setNames(gate_name_sh_vec)
}

#' @title Get names of gates as gate types and combination methods
#'
#' @inheritParams .get_gate_list # fdr
#' @inheritParams plot_cp # gate_list
#'
#' @return Character vector where each element is a gate type and combination method (shorthand, separated by an underscore) for
#' which we have calculated a gate on these data.
.get_gate_name_vec_sh <- function(gate_tbl) {
  gate_tbl <- gate_tbl |>
    dplyr::group_by(gate_type, gate_combn) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  paste0(gate_tbl$gate_method, "_", gate_tbl$gate_combn)
  # purrr::map(unique(gate_tbl$gate_name), function(cp){
  #  paste0(cp, "_", names(gate_list[[cp]]))
  # }) |>
  #  unlist()
}

#' @title Get a table of gates for each sample for each gate for a given gating population
#'
#' @param gate_list
#'
.get_gate_tbl <- function(gate_list, ind_batch_list) {
  ind_batch_lab_vec <- c()
  for (i in seq_along(ind_batch_list)) {
    ind_batch_lab_vec <- c(ind_batch_lab_vec, stats::setNames(
      rep(i, length(ind_batch_list[[i]])),
      ind_batch_list[[i]]
    ))
  }
  purrr::map_df(names(gate_list), function(gate_type) {
    purrr::map_df(names(gate_list[[gate_type]]), function(gate_combn) {
      tibble::tibble(
        gate_name = paste0(.env$gate_type, "_", .env$gate_combn),
        gate_type = .env$gate_type,
        gate_combn = .env$gate_combn,
        ind = names(gate_list[[.env$gate_type]][[.env$gate_combn]]) |>
          as.integer(),
        batch = ind_batch_lab_vec[ind],
        gate = gate_list[[.env$gate_type]][[.env$gate_combn]]
      ) |>
        dplyr::select(gate_name:gate_combn, batch, ind, gate)
    })
  })
}


#' @title Plot histogram
#'
#' @description
#' Plots a histogram of a specific channel. Allows removal of the cells with the lowest value, and display of cutpoints.
#'
#' @inheritParams plot_cp
#' @inheritParams .check_if_high
#' @param bin integer. Number of bins of histogram. Default is 300.
#' @param min numeric. Minimum value to plot for histogram. If \code{NULL}, then all values are plotted. Default is \code{NULL}.
#' @param chnl character. Name of column in \code{ex} to plot.
#' @param axis_lab named character vector. Names are names of channels and values are corresponding names of markers. If provided, then the x-axis label is channged from the channel name (in \code{chnl}) to the common name. If \code{NULL}, then the axis label is kept as the channel name. Default is \code{NULL}.
#' @param col named character vector. Names are names of cutpoints in \code{cp}, and values are corresponding colours. If provided, then the colours of the vertical lines indicating the cutpoints are set to these colours. If \code{NULL}, then the colours are set automatically using \code{RColorBrewer}. Default is \code{NULL}.
#' @param cp_lab named character vector. Names are names of cutpoints, and values are common names of cutpoints to be used in the cutpoint type legend. If provided, then the cutpoint type legend's values are set to these. If \code{NULL}, then the untransformed cutpoint names are used. Default is \code{NULL}.
#'
#' @return A \code{ggplot2} object.
.plot_hist <- function(ex, plot_var, gate_tbl, exc_min = FALSE,
                       axis_lab = NULL, col = NULL, gate_lab = NULL,
                       lb = NULL, title) {
  # plot histogram
  plot_tbl <- ex[, plot_var]
  colnames(plot_tbl) <- "x"
  if (!is.null(lb)) plot_tbl <- plot_tbl |> dplyr::filter(x > lb)
  if (exc_min) plot_tbl <- plot_tbl |> dplyr::filter(x > min(x))
  if (nrow(plot_tbl) == 0) {
    return(ggplot())
  }

  hist_obj <- hist(plot_tbl[["x"]], plot = FALSE)
  bin_vec <- max(c(hist_obj$breaks |> length() |> multiply_by(2), 2))

  ggplot(plot_tbl, aes(x)) +
    geom_histogram(
      fill = "gray65", col = "gray85",
      bins = bin_vec
    ) +
    cowplot::background_grid(major = "x", minor = "x") +
    labs(y = "Count", x = axis_lab[plot_var]) +
    geom_vline(
      data = gate_tbl,
      aes(xintercept = gate, col = gate_name),
      size = 1.5
    ) +
    scale_colour_manual(
      values = col,
      labels = gate_lab
    ) +
    theme(legend.title = element_blank())
}

#' @title Create hex plots for the cut and high channels
#'
#' @description
#' Create 2D hex plots for the cut channel with each of the high channels.
#'
#' @inheritParams plot_cp
#' @inheritParams
#'
#' @return A list of ggplot2 objects, one for each
#' of the high channels.
.plot_2d <- function(ex, plot_var_x, gate_tbl, high,
                     axis_lab = NULL, col = NULL,
                     gate_lab = NULL, exc_min = TRUE,
                     title) {
  # create a list of plots for each channel "high" is specified for
  high_ind_vec <- setdiff(seq_along(high), which(names(high) == plot_var_x))
  high_name_vec <- names(high)[high_ind_vec]
  purrr::map(high_ind_vec, function(i) {
    # get table to create hex's from
    high_chnl <- high[i]
    points_tbl <- tibble::tibble(
      cut = ex[[plot_var_x]],
      high = ex[[names(high)[i]]]
    )
    if (exc_min) {
      points_tbl <- points_tbl |> dplyr::filter(cut > min(cut) |
        high > min(high))
    }

    # create table to display line for
    high_line_tbl <- tibble::tibble(
      type = "high",
      chnl = names(high)[i],
      min = high[i]
    )
    # base plot
    ggplot() +
      geom_hex(
        aes(cut, high),
        points_tbl
      ) +
      scale_fill_continuous(
        trans = "log10",
        type = "viridis"
      ) +
      # geom_hline(aes(yintercept = min),
      #           high_line_tbl, size = 2, col = 'black') +
      geom_vline(
        data = gate_tbl,
        aes(xintercept = gate, col = gate_name),
        size = 1.5
      ) +
      scale_colour_manual(
        values = col,
        labels = gate_lab
      ) +
      theme(legend.title = element_blank()) +
      labs(
        x = axis_lab[plot_var_x],
        y = axis_lab[names(high)[i]]
      )
  }) |>
    stats::setNames(paste0("plot2d-", axis_lab[high_name_vec]))
}

#' @title Plot frequencies and counts by a range of cutpoints
#'
#' @description
#' Plot frequencies and counts for any given cutpoint in a range. Also plots
#' the standardised frequencies and counts.
#'
#' @param ex dataframe.
#' @inheritParams get_cp
.plot_count_stats_by_cp <- function(ex, plot_var, gate_tbl,
                                    cp, lb, axis_lab, col,
                                    gate_lab, remove_uns = FALSE) {
  # --------------------
  # Preparation
  # --------------------

  # return if an unstim sample only
  if (length(setdiff(unique(ex[["is_uns"]]), TRUE)) == 0) {
    out_list <- purrr::map(1:4, function(i) {
      ggplot() |>
        stats::setNames(paste0("plot_gate_stat-", "p_freq"))
    })
    return(out_list)
  }


  # get ex for stim and uns dbns
  ex_uns <- ex |> dplyr::filter(is_uns)
  ex_stim <- ex |> dplyr::filter(!is_uns)

  # get ecdfs for stim and uns dbns
  ecdf_uns <- ecdf(ex_uns[["cut"]])
  ecdf_stim <- ecdf(ex_stim[["cut"]])
  n_cell <- nrow(ex_stim)

  # get cutpoints at which to calculate counts and frequencies for
  max_stim <- max(ex_stim[["cut"]])
  if (max_stim > lb + 3) {
    plot_range_vec <- c(lb, max(ex_stim[["cut"]])) |> stats::setNames(NULL)
    stat_cp_vec <- seq(plot_range_vec[1], plot_range_vec[2] + 3, length.out = 50)
  } else {
    # plot_range_vec <- c(lb, lb + 3) |> stats::setNames(NULL)
    stat_cp_vec <- seq(lb, max_stim + max_stim / 25, length.out = 50)
  }


  # --------------------
  # Calculate counts and frequencies
  # --------------------

  # calculate counts and frequencies at each cutpoint
  # for current stim and unstim
  cp_working_stats_tbl <- purrr::map_df(stat_cp_vec, function(cp) {
    count_stim <- (1 - ecdf_stim(cp)) * n_cell
    count_uns <- (1 - ecdf_uns(cp)) * n_cell
    count_diff <- count_stim - count_uns
    prop_uns <- count_uns / count_stim
    prob_uns <- count_uns / n_cell
    prob_stim <- count_stim / n_cell
    prob_diff <- prob_stim - prob_uns
    prob_diff_std <- (prob_diff) /
      sqrt(n_cell * prob_uns * (1 - prob_uns) +
        n_cell * prob_stim * (1 - prob_uns))
    freq_stim <- count_stim / n_cell * 1e2
    freq_uns <- count_uns / n_cell * 1e2
    freq_diff <- freq_stim - freq_uns
    data.frame(
      cp = cp, count_stim = count_stim, count_uns = count_uns,
      count_diff = count_diff, prop_uns = prop_uns,
      freq_stim = freq_stim, freq_uns = freq_uns,
      freq_diff = freq_diff,
      prob_stim = prob_stim, prob_uns = prob_uns,
      prob_diff = prob_diff, prob_diff_std = prob_diff_std
    )
  }) |>
    tibble::as_tibble()

  if (remove_uns) {
    cps_working_stats_tbl <- cps_working_stats_tbl |>
      dplyr::select(-c(
        count_stim, count_diff, prop_uns, freq_uns,
        freq_diff, prob_uns, prob_diff, prob_diff_std
      ))
  }


  # --------------------
  # Plots
  # --------------------

  # Proportion unstim
  # ---

  if (!remove_uns) {
    p_prop <- ggplot(cp_working_stats_tbl, aes(x = cp, y = prop_uns)) +
      geom_line(size = 2) +
      geom_hline(yintercept = 0, size = 1) +
      scale_y_continuous(breaks = function(x) .get_axis_breaks(x, by = 0.1)) +
      cowplot::background_grid(minor = "xy") +
      labs(y = "Ratio of unstim count to stim count")
  } else {
    p_prop <- ggplot()
  }


  # Count
  # ---

  if (!remove_uns) {
    plot_tbl <- cp_working_stats_tbl |>
      dplyr::select(-prop_uns) |>
      tidyr::pivot_longer(count_stim:count_diff,
        names_to = "count_type",
        values_to = "count"
      )
    p_count <- ggplot(plot_tbl, aes(x = cp)) +
      geom_hline(yintercept = 0, size = 1) +
      geom_line(aes(y = count, linetype = count_type), size = 2) +
      cowplot::background_grid(minor = "xy") +
      labs(y = "Count") +
      theme(legend.key.width = unit(20, "mm")) +
      scale_linetype_manual(
        values = c(
          "count_uns" = "11",
          "count_stim" = "31",
          "count_diff" = "solid"
        ),
        labels = c(
          "count_stim" = "Stim",
          "count_uns" = "Unstim",
          "count_diff" = "Stim - unstim"
        )
      ) +
      theme(legend.title = element_blank())
  } else {
    p_count <- ggplot(
      cp_working_stats_tbl,
      aes(x = cp)
    ) +
      geom_hline(yintercept = 0, size = 1) +
      geom_line(aes(y = count_stim), size = 2) +
      cowplot::background_grid(minor = "xy") +
      labs(y = "Count")
  }

  # Frequency
  # ---

  if (!remove_uns) {
    p_freq <- ggplot(
      cp_working_stats_tbl |>
        dplyr::select(-prop_uns) |>
        tidyr::pivot_longer(freq_stim:freq_diff,
          names_to = "freq_type",
          values_to = "freq"
        ),
      aes(x = cp, linetype = freq_type)
    ) +
      geom_hline(yintercept = 0, size = 1) +
      geom_line(aes(y = freq), size = 2) +
      cowplot::background_grid(minor = "xy") +
      theme(
        legend.key.width = unit(20, "mm"),
        legend.title = element_blank()
      ) +
      scale_linetype_manual(
        values = c(
          "freq_uns" = "11",
          "freq_stim" = "31",
          "freq_diff" = "solid"
        ),
        labels = c(
          "freq_stim" = "Stim",
          "freq_uns" = "Unstim",
          "freq_diff" = "Stim - unstim"
        )
      ) +
      labs(y = "Frequency")
  } else {
    p_freq <- ggplot(
      cp_working_stats_tbl,
      aes(x = cp)
    ) +
      geom_hline(yintercept = 0, size = 1) +
      geom_line(aes(y = freq_stim), size = 2) +
      cowplot::background_grid(minor = "xy") +
      labs(y = "Frequency")
  }


  # Standardised proportions
  # ---

  if (!remove_uns) {
    p_prob_std <- ggplot(
      cp_working_stats_tbl,
      aes(x = cp, y = prob_diff_std)
    ) +
      geom_line(size = 2) +
      geom_hline(yintercept = 0, size = 1) +
      cowplot::background_grid(minor = "xy") +
      labs(y = "Standardised difference in proportion")
  } else {
    p_prob_std <- ggplot()
  }

  # Change axis labels and add cutpoints, if provided
  plot_list <- list( # p_prop = p_prop,
    p_count = p_count,
    p_freq = p_freq
  )
  # p_prob_std = p_prob_std,
  # p_stim_vs_uns = list())

  plot_list <- purrr::map(plot_list, function(p) {
    # return straight away if plot essentially empty
    if (length(p$layers) == 0) {
      return(p)
    }

    if (!is.null(axis_lab)) p <- p + labs(x = axis_lab[plot_var])

    # add cutpoints
    p +
      geom_vline(
        data = gate_tbl,
        aes(xintercept = gate, col = gate_name),
        size = 1.5
      ) +
      scale_colour_manual(
        values = col,
        labels = gate_lab
      ) +
      theme(legend.title = element_blank())
  }) |>
    stats::setNames(paste0("plot_gate_stat-", c("p_freq")))

  # --------------------
  # Output
  # --------------------

  plot_list
}

#' @title Plot proportin of positive cells for a given marker expression level
#'
#' @description
#' Plots the proportion of cells for the marker for which a cutpoint is sought that are high for the specified markers.
#'
#' @inheritParams .plot_hist
#' @inheritParams .check_if_high
#' @param rug numeric between 0 and 1 (inclusive). If provided, then a rug plot with this proportion of the cells is plotted. Note that this is from the high end, so for example if \code{rug = 0.1}, then a rug for the cells from the 90th quantile upwards is plotted. Default is 0.05.
#'
#' @return A ggplot2 object.
.plot_prop <- function(high_ind_tbl, cp,
                       cut,
                       axis_lab = NULL, col = NULL, cp_lab = NULL,
                       rug = 0.05) {
  # get binned marker expression
  bin_tbl <- high_ind_tbl |>
    dplyr::mutate(bin = cut(cut, breaks = 50)) |>
    dplyr::mutate(comma_loc = stringr::str_locate(bin, ",")[, 1]) |>
    dplyr::mutate(
      lb = stringr::str_sub(bin, 2, comma_loc - 1) |> as.numeric(),
      ub = stringr::str_sub(bin, comma_loc + 1, -2) |> as.numeric()
    ) |>
    dplyr::select(-comma_loc) |>
    dplyr::mutate(mid = (lb + ub) / 2)

  # get proportion positive for each bin
  sum_tbl <- bin_tbl |>
    dplyr::group_by(bin, lb, ub, mid) |>
    dplyr::summarise(
      count_high = sum(high),
      prop_high = count_high / dplyr::n(),
      count_total = dplyr::n(),
      prop_total = count_total / nrow(high_ind_tbl)
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(desc(bin)) |>
    dplyr::mutate(
      count_cum = cumsum(count_total),
      prop_cum = cumsum(prop_total)
    ) |>
    dplyr::arrange(desc(bin))

  # get smoothed proportion positive for each bin
  smooth_vec <- vapply(3:(nrow(sum_tbl) - 2), function(i) {
    ind_vec <- c(i - 2, i - 1, i, i + 1, i + 2)
    prop_high_vec <- sum_tbl[["prop_high"]][ind_vec]
    prox_weight_vec <- c(0.05, 0.25, 0.4, 0.25, 0.05)
    sum(prop_high_vec * prox_weight_vec)
  }, numeric(1))

  high_vec <- sum_tbl[["prop_high"]]
  n_bin <- nrow(sum_tbl)
  smooth_vec <- c(
    0.7 * high_vec[1] + 0.25 * high_vec[2] + 0.05 * high_vec[3],
    0.3 * high_vec[1] + 0.4 * high_vec[2] + 0.25 * high_vec[3] + 0.05 * high_vec[4],
    smooth_vec,
    0.3 * high_vec[n_bin] + 0.4 * high_vec[n_bin - 1] +
      0.25 * high_vec[n_bin - 2] + 0.05 * high_vec[n_bin - 2],
    0.7 * high_vec[n_bin] + 0.25 * high_vec[n_bin - 1] + 0.05 * high_vec[n_bin - 2]
  )

  smooth_tbl <- tibble::tibble(mid = sum_tbl[["mid"]], prop_smooth = smooth_vec)

  # base plot
  p <- ggplot() +
    geom_bar(
      data = sum_tbl,
      aes(x = mid, y = prop_high),
      stat = "identity",
      fill = "gray65",
      col = "gray80"
    ) +
    geom_line(
      data = smooth_tbl,
      aes(x = mid, y = prop_smooth),
      size = 1.5
    )

  if (rug != 0) {
    p <- p +
      geom_rug(
        data = high_ind_tbl |>
          dplyr::filter(.data$cut >= quantile(.data$cut, 1 - rug)[[1]]),
        aes(x = cut)
      ) +
      labs(x = cut, y = "Proportion high") +
      scale_y_continuous(limits = c(0, 1))
  }

  # plot axis label, if needed
  if (!is.null(axis_lab) & !is.null(cut)) p <- p + labs(x = axis_lab[cut])

  if (is.null(cp)) {
    return(p)
  }

  # add cutpoints
  cp_plot_tbl <- tibble::tibble(cp = cp, cp_type = names(cp))
  if (is.null(col)) {
    col <- RColorBrewer::brewer.pal(length(cp), "Set1") |>
      stats::setNames(names(cp))
  }
  if (is.null(cp_lab)) cp_lab <- stats::setNames(names(cp$cp), names(cp$cp))

  p +
    geom_vline(
      data = cp_plot_tbl,
      aes(xintercept = cp, col = cp_type),
      size = 1.5
    ) +
    scale_colour_manual(
      values = col,
      labels = cp_lab
    ) +
    theme(legend.title = element_blank())
}


smooth_line <- function(num) {
  smooth_vec <- vapply(3:(length(num) - 2), function(i) {
    ind_vec <- c(i - 2, i - 1, i, i + 1, i + 2)
    prop_high_vec <- num[ind_vec]
    prox_weight_vec <- c(0.05, 0.25, 0.4, 0.25, 0.05)
    sum(prop_high_vec * prox_weight_vec)
  }, numeric(1))

  n_bin <- length(num)
  smooth_vec <- c(
    0.7 * num[1] + 0.25 * num[2] + 0.05 * num[3],
    0.3 * num[1] + 0.4 * num[2] + 0.25 * num[3] + 0.05 * num[4],
    smooth_vec,
    0.3 * num[n_bin] + 0.4 * num[n_bin - 1] +
      0.25 * num[n_bin - 2] + 0.05 * num[n_bin - 2],
    0.7 * num[n_bin] + 0.25 * num[n_bin - 1] + 0.05 * num[n_bin - 2]
  )

  smooth_vec
}


#' @title Plot log-likelihood
#'
#' @description
#' Plots log-likelihood and displays cutpoints.
#'
#' @inheritParams plot_cp
#' @param ll numeric vector. Log-likelihood values.
#' @param cps numeric vector. Cutpoints at which the log-likelihood was evalutated.
#'
#' @return A ggplot2 object.
.plot_ll <- function(ll, cps, cp, cut = NULL, axis_lab = NULL, col = NULL,
                     cp_lab = NULL) {
  ll_plot_tbl <- data.frame(ll = ll, cps = cps)
  p <- ggplot() +
    geom_line(data = ll_plot_tbl, aes(x = cps, y = ll)) +
    labs(y = "Log-likelihood")

  if (!is.null(axis_lab) & !is.null(cut)) p <- p + labs(x = axis_lab[cut])

  if (is.null(cp)) {
    return(p)
  }

  # add cutpoints
  cp_plot_tbl <- tibble::tibble(cp = cp, cp_type = names(cp))
  if (is.null(col)) {
    col <- RColorBrewer::brewer.pal(length(cp), "Set1") |>
      stats::setNames(names(cp))
  }
  if (is.null(cp_lab)) cp_lab <- stats::setNames(names(cp$cp), names(cp$cp))

  p <- p +
    geom_vline(
      data = cp_plot_tbl,
      aes(xintercept = cp, col = cp_type),
      size = 1.5
    ) +
    scale_colour_manual(
      values = col,
      labels = cp_lab
    ) +
    theme(legend.title = element_blank())

  p
}


#' @title Plot histogram of all values above a given level
#'
#' @description
#' Plots a histogram of all values above a certain value. May add cutpoints.
#'
#' @inheritParams plot_cp
#' @inheritParams .check_if_high
#'
#' @param bins integer. Number of bins of the histogram.
#'
#' @return A ggplot2 object.
.plot_dcp_hist <- function(ex, cut, min, bins = 50) {
  hist_tbl <- tibble::tibble(x = ex[[cut]][ex[[cut]] > cp_scp])

  cut_line_tbl <- tibble::tibble(
    type = names(cp),
    xint = cp
  )
  col_vec <- RColorBrewer::brewer.pal(length(cp), "Set1") |> stats::setNames(names(cp))

  ggplot(hist_tbl, aes(x = x)) +
    geom_histogram(bins = bins, fill = "gray20", col = "gray85") +
    geom_vline(aes(xintercept = xint, col = type),
      cut_line_tbl,
      size = 2
    ) +
    labs(x = lab[cut]) +
    scale_colour_manual(
      values = col_vec,
      name = "Gate type"
    )
}





# not used anymore
.plot_cp_pwder <- function(high_ind_tbl, cp) {
  cp_pwder_tbl <- .get_cp_pwder_tbl(
    high_ind_tbl = high_ind_tbl,
    cp_scp = cp[["scp"]]
  )

  p_der <- ggplot(cp_pwder_tbl, aes(x = cut, y = der)) +
    geom_line(size = 2) +
    geom_vline(
      xintercept = cp[["pwder"]],
      size = 2, col = "red"
    )

  # plot probability at each point
  plot_prop <- .plot_prop(high_ind_tbl, cp = cp)

  plot_prop$layers[[2]] <- NULL

  p_prop <- plot_prop +
    geom_line(aes(x = cut, y = pred),
      cp_pwder_tbl,
      size = 2
    )
  list(
    "p_der" = p_der,
    "p_prop" = p_prop
  )
}


#' @title Plot proportion positive with fitted log reg line
#'
#' @description
#' This is the same as \code{plot_prop}, except that instead of a
#' smoothed line of proportion positive, the fitted values
#' from the logistic regression used to obtain the mid-probability cut
#' is displayed instead.
#'
#' @inheritParams .plot_prop.p
#' @param plot_prop ggplot2 object. Plot returned by \code{.plot_prop}.
#'
#' @return A ggplot2 object.
.plot_logit <- function(high_ind_tbl, cp, plot_prop) {
  # get table to model - all values above changepoint
  mod_tbl <- high_ind_tbl |>
    dplyr::filter(cut >= cp["scp"])

  # check if too little data to model with here
  if (nrow(mod_tbl) < 5 || length(unique(high_ind_tbl$high)) == 1) {
    # print("Too few obs to make logistic plot.")
    return(ggplot())
  }

  # model the probability of positivity above this point
  if (nrow(mod_tbl) < 40) {
    fit_pw <- glm(high ~ cut,
      family = binomial, data = mod_tbl
    )
  } else { # if(nrow(mod_tbl) < 35){
    fit_pw <- glm(high ~ splines::ns(cut, df = 3),
      family = binomial, data = mod_tbl
    )
  } # else if(nrow(mod_tbl) < 100){
  # fit_pw <- glm(high ~ splines::ns(cut, df = 4) ,
  #              family = binomial, data = mod_tbl)
  # } else{
  #  fit_pw <- glm(high ~ splines::ns(cut, df = 5) ,
  #                family = binomial, data = mod_tbl)
  # }

  # get predictions over a range of cut values
  pred_tbl <- tibble::tibble(cut = seq(min(mod_tbl$cut), max(mod_tbl$cut)))
  pred_vec <- predict(fit_pw, pred_tbl, type = "response")
  pred_tbl <- pred_tbl |> dplyr::mutate(pred = pred_vec)

  # find prediction in between middle and max
  min <- mean(high_ind_tbl |> dplyr::filter(
    cut < cp["scp"]
  ) |>
    dplyr::pull("high"))
  max <- max(pred_tbl$pred)
  mid_prob <- mean(c(max, min))

  # plot probability at each point
  plot_prop$layers[[2]] <- NULL

  # add logistic regression fit
  plot_prop +
    geom_line(aes(x = cut, y = pred),
      pred_tbl,
      size = 2
    )
}







#' @title Get axis breaks given axis limits
#'
#' @description
#' This function is intended to be provided to the
#' breaks parameter of the ggplot2::scale_[]_continuous function.
#' It automatically determines the break points given the limits of
#' the axis.
#'
#' @param lims numeric vector. The elements are the minimum and maximum
#' of the range, respectively (in that order, to be clear).
#' @param by numeric. Gap between breaks.
#'
#' @return A numeric vector.
.get_axis_breaks <- function(lims, by = 0.1) {
  exponent <- -log10(by)
  factor <- 10^exponent
  min <- floor(lims[1] * factor) / factor
  max <- ceiling(lims[2] * factor) / factor
  seq(min, max, by = by)
}
