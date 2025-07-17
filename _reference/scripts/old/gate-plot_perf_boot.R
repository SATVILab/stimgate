#' @title Plot gate, statistic and performance and create tables thereof for gates derived from a single population
plot_cp_boot <- function(gate_tbl,
                         params,
                         pop_sub,
                         plot = TRUE,
                         path_project) {
  # =========================
  # Get stats
  # =========================

  print("getting boot stats")

  if (is.null(gate_tbl)) {
    dir_base <- stimgate_dir_base_create(
      dir_base_init = path_project, params = params
    )
    path_tbl <- file.path(dir_base, "stats", "gate_stats_tbl_boot")
    if (!file.exists(path_tbl)) {
      print("no boot tbl stats found after skipping calculation thereof")
      return(invisible(TRUE))
    }
    gate_stats_tbl <- readRDS(path_tbl)
  } else {
    # table
    gate_stats_tbl <- .get_gate_stats(
      params = params,
      gate_tbl = gate_tbl,
      pop_sub = pop_sub,
      boot = TRUE,
      path_project = path_project
    )
  }



  # =========================
  # Performance
  # =========================

  if (!is.null(params$pop_man_sub) && plot) {
    # plots
    print("getting boot perf plots")

    # get labels for gates
    gate_lab_vec <- .get_gate_lab_vec(gate_tbl = gate_stats_tbl)

    .plot_gate_perf_boot(
      gate_stats = gate_stats_tbl,
      params = params,
      gate_lab = gate_lab_vec,
      path_project = path_project
    )
  }

  invisible(TRUE)
}

.plot_gate_perf_boot <- function(gate_stats,
                                 params,
                                 gate_lab,
                                 path_project) {
  # Calculate signal-to-noise
  # ---------------------
  .get_s2n <- function(stat, type = "mean") {
    if (type == "median") {
      stat_mad <- mad(stat)
      stat_med <- median(stat)
      if (stat_med == 0) {
        return(0)
      }
      if (stat_mad == 0) {
        return(NA)
      }
      return(stat_med / stat_mad)
    } else if (type == "mix") {
      stat_med <- median(stat)
      stat_se <- sd(stat)
      if (stat_med == 0) {
        return(0)
      }
      if (is.na(stat_se)) {
        return(NA)
      }
      return(stat_med / stat_se)
    } else if (type == "mean") {
      stat_mean <- mean(stat)
      stat_se <- sd(stat)
      if (stat_mean == 0) {
        return(0)
      }
      if (is.na(stat_se)) {
        return(NA)
      }
      return(stat_mean / stat_se)
    }
  }

  # Calculate proportion of signal left after subtracting one se of mean
  # ---------------------------------------------------
  .get_se_prop <- function(stat) {
    if (mean(stat) == 0) {
      return(0)
    }
    max(0, (mean(stat) - sd(stat)) / mean(stat))
  }

  plot_list <- purrr::map(unique(gate_stats$pop), function(pop) {
    gate_stats <- gate_stats |> dplyr::filter(.data$pop == .env$pop)

    s2n_tbl <- gate_stats |>
      dplyr::filter(
        !is.na(count_bs),
        !gate_type == "man"
      ) |>
      dplyr::mutate(
        freq_bs = pmax(0, freq_bs),
        count_bs = pmax(0, count_bs)
      ) |>
      dplyr::group_by(gate_name, gate_type, gate_combn, ind, stim, batch_sh) |>
      dplyr::summarise(freq_bs_s2n = .get_s2n(stat = freq_bs, type = "median")) |>
      dplyr::ungroup()

    x_seq <- seq(0.33, 10, length.out = 1e3)
    plot_tbl <- purrr::map_df(s2n_tbl$gate_name, function(gate_name_curr) {
      s2n_tbl_curr <- s2n_tbl |>
        dplyr::filter(gate_name == gate_name_curr)
      ecdf_curr <- ecdf(s2n_tbl_curr$freq_bs_s2n)
      tibble::tibble(gate_name = gate_name_curr, x = x_seq, surv = 1 - ecdf_curr(x_seq))
    })
    plot_tbl <- plot_tbl |>
      dplyr::filter(surv > 0 | x_seq >= 0.25)
    p_s2n_ecdf <- ggplot(plot_tbl, aes(x = x, y = surv, col = gate_name)) +
      cowplot::background_grid(major = "xy", minor = "xy") +
      geom_line(size = 2) +
      scale_x_continuous(breaks = c(0.33, 0.67, 1, 1.5, 2, 3, 4, 5, 7.5, 10)) +
      labs(x = "Signal-to-noise ratio", y = "Proportion of samples with\n signal-to-noise ratio larger than a given value") +
      lims(y = c(0, 1)) +
      scale_colour_brewer(
        palette = "Set1",
        labels = gate_lab
      ) +
      theme(legend.title = element_blank())

    se_prop_tbl <- gate_stats |>
      dplyr::filter(
        !is.na(count_bs),
        !gate_type == "man"
      ) |>
      dplyr::mutate(
        freq_bs = pmax(0, freq_bs),
        count_bs = pmax(0, count_bs)
      ) |>
      dplyr::group_by(gate_name, gate_type, gate_combn, ind, stim, batch_sh) |>
      dplyr::summarise(freq_bs_seprop = .get_se_prop(freq_bs)) |>
      dplyr::ungroup()

    # Plot
    # ------------------------------

    batch_vec <- unique(s2n_tbl$batch_sh)
    shape_vec <- stats::setNames(
      c(
        15, 16, 17, 18, 14, 10, 6,
        7, 8, 12, 23
      )[seq_along(batch_vec)],
      batch_vec
    )
    library(ggforce)
    p_seprop_freq_bs <- ggplot(
      se_prop_tbl,
      aes(x = gate_name, y = freq_bs_seprop)
    ) + # , fill = gate_name)) +
      cowplot::background_grid(major = "y", minor = "y") +
      geom_hline(yintercept = 0, linetype = "dotted", size = 0.5) +
      # geom_violin(aes(col = gate_name), size = 1, scale = FALSE) +
      geom_boxplot(outlier.size = -1, alpha = 0.2, size = 1) +
      geom_sina(aes(col = stim, shape = batch_sh),
        size = 3, scale = FALSE,
        maxwidth = 0.5
      ) +
      scale_shape_manual(
        values = shape_vec,
        name = "Batch"
      ) +
      scale_colour_brewer(
        palette = "Set1",
        name = "Stimulation"
      ) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_blank()
      ) +
      labs(x = "Gate method", y = "(Mean - standard error of mean)/mean\n of background-subtracted frequency for individual stim samples") +
      scale_x_discrete(labels = gate_lab)

    p_s2n_freq_bs <- ggplot(
      s2n_tbl |>
        dplyr::filter(!is.na(freq_bs_s2n)) |>
        dplyr::mutate(freq_bs_s2n = pmax(freq_bs_s2n, 0.25)),
      aes(x = gate_name, y = freq_bs_s2n)
    ) + # , fill = gate_name)) +
      cowplot::background_grid(major = "y", minor = "y") +
      geom_hline(yintercept = 0, linetype = "dotted", size = 0.5) +
      # geom_violin(aes(col = gate_name), size = 1, scale = FALSE) +
      geom_boxplot(outlier.size = -1, alpha = 0.2, size = 1) +
      geom_sina(aes(col = stim, shape = batch_sh),
        size = 3, scale = FALSE,
        maxwidth = 0.5
      ) +
      # geom_bar(stat = 'identity') +
      scale_shape_manual(
        values = shape_vec,
        name = "Batch"
      ) +
      scale_colour_brewer(
        palette = "Set1",
        name = "Stimulation"
      ) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_blank()
      ) +
      labs(x = "Gate method", y = "Signal-to-noise ratio (mean/sd) of\nbackground-subtracted frequency for individual stim samples") +
      scale_y_continuous(trans = "log2") +
      scale_x_discrete(labels = gate_lab)

    list(
      p_s2n_freq_bs = p_s2n_freq_bs,
      p_seprop_freq_bs = p_seprop_freq_bs,
      p_s2n_ecdf = p_s2n_ecdf
    )
  }) |>
    stats::setNames(unique(gate_stats$pop))

  # Save
  # -------------------------------------

  print("saving perf plots")
  .save_plot_gate_stats_or_perf(
    plot_list = plot_list,
    params = params,
    plot_type = "perf",
    unlink = TRUE,
    path_project = path_project
  )
}
