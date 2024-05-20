#' @title Plot auto vs man count and frequency selected statistics for multiple populations
.plot_gate_stats <- function(gate_stats = NULL,
                             params = NULL,
                             stat = NULL,
                             col = NULL,
                             gate_lab = NULL,
                             path_project) {
  # get marker name
  marker <- params$chnl_lab[params$cut][[1]]

  # get plots across pops and stats
  plot_list <- .plot_gate_stats_pop_stat(
    gate_stats = gate_stats,
    stat = stat,
    col = col,
    gate_lab = gate_lab,
    marker = marker,
    pop_gate = params$pop_gate
  )

  # save
  print("saving stats plots")
  .save_plot_gate_stats_or_perf(
    plot_list = plot_list,
    params = params,
    plot_type = "stats",
    path_project = path_project
  )

  invisible(TRUE)
}


#' @title Plot selected statistics for multiple populatoins
.plot_gate_stats_pop_stat <- function(gate_stats,
                                      stat = NULL,
                                      col,
                                      gate_lab,
                                      marker,
                                      pop_gate) {
  if (is.null(stat)) stat <- c("count", "freq", "count_bs", "freq_bs")

  plot_list <- list()

  # print("auto vs man plots")
  purrr::map(unique(gate_stats$pop), function(pop_curr) {
    pop_title <- ifelse(pop_curr == "gate", pop_gate, paste0(pop_gate, "/", pop_curr))
    # print(pop_curr)
    stat_auto_vs_man_list <- purrr::map(stat, function(stat_curr) {
      # print(stat_curr)
      .plot_gate_stats_auto_vs_man(
        gate_stats = gate_stats,
        stat = stat_curr,
        pop = pop_curr,
        col = col,
        gate_lab = gate_lab,
        marker = marker,
        pop_title = pop_title
      )
    }) |>
      stats::setNames(paste0("auto_vs_man-", stat))
    # print("stim vs uns plots")
    stat_stim_vs_uns_list_dodge_gm <- purrr::map(stat, function(stat_curr) {
      # print(stat_curr)
      .plot_gate_stats_stim_vs_uns_dodge_gm(
        gate_stats = gate_stats,
        stat = stat_curr,
        pop = pop_curr,
        col = col,
        gate_lab = gate_lab,
        marker = marker,
        pop_title = pop_title
      )
    }) |>
      stats::setNames(paste0("stim_vs_uns_dodge_gm-", stat))
    stat_auto_vs_man_list <- stat_auto_vs_man_list |>
      append(stat_stim_vs_uns_list_dodge_gm)

    stat_stim_vs_uns_list_dodge_stim <- purrr::map(stat, function(stat_curr) {
      # print(stat_curr)
      .plot_gate_stats_stim_vs_uns_dodge_stim(
        gate_stats = gate_stats,
        stat = stat_curr,
        pop = pop_curr,
        col = col,
        gate_lab = gate_lab,
        marker = marker,
        pop_title = pop_title
      )
    }) |>
      stats::setNames(paste0("stim_vs_uns_dodge_stim-", stat))

    stat_auto_vs_man_list <- stat_auto_vs_man_list |>
      append(stat_stim_vs_uns_list_dodge_stim)

    stat_stim_vs_uns_list_jc <- purrr::map(
      c("count", "freq"),
      function(stat_curr) {
        # print(stat_curr)
        .plot_gate_stats_stim_vs_uns_jc(
          gate_stats = gate_stats,
          stat = stat_curr,
          pop = pop_curr,
          col = col,
          gate_lab = gate_lab,
          marker = marker,
          pop_title = pop_title
        )
      }
    ) |>
      stats::setNames(paste0("stim_vs_uns_jc-", c("count", "freq")))
    stat_auto_vs_man_list |> append(stat_stim_vs_uns_list_jc)
  }) |>
    stats::setNames(unique(gate_stats$pop))
}


#' @title Plot a single statistic for a single population
.plot_gate_stats_auto_vs_man <- function(gate_stats,
                                         pop,
                                         stat,
                                         col,
                                         gate_lab,
                                         marker,
                                         pop_title) {
  # remember to change plot_tbl to cp_stats later
  plot_tbl <- gate_stats |> dplyr::filter(.data$pop == .env$pop)

  # adjust column name for programming sake
  col_name_vec <- colnames(plot_tbl)
  stat_ind <- which(col_name_vec == stat)
  col_name_vec[stat_ind] <- "stat"
  colnames(plot_tbl) <- col_name_vec

  # remove NA entries
  plot_tbl <- plot_tbl |> dplyr::filter(!is.na(stat))

  # get axis limits
  max_axis_lim_init <- max(plot_tbl[["stat"]])
  min_axis_lim_init <- min(plot_tbl[["stat"]])
  range_length <- max_axis_lim_init - min_axis_lim_init
  min_axis_lim <- min_axis_lim_init - 0.05 * range_length
  max_axis_lim <- max_axis_lim_init + 0.05 * range_length
  axis_lim_vec <- c(min_axis_lim, max_axis_lim)

  # get manual cut values
  man_stat_temp_tbl <- plot_tbl |>
    dplyr::filter(gate_name == "man_no")
  man_stat_vec <- stats::setNames(
    man_stat_temp_tbl[["stat"]], man_stat_temp_tbl[["fcs"]]
  )

  # get plot tbl
  plot_tbl <- plot_tbl |>
    dplyr::filter(!gate_name == "man") |>
    dplyr::mutate(man_stat = man_stat_vec[fcs])

  # get shape label vector
  batch_vec <- unique(plot_tbl$batch_sh)
  shape_vec <- stats::setNames(
    c(15, 16, 17, 18, 14, 10, 6, 7, 8, 12, 23)[seq_along(batch_vec)],
    batch_vec
  )

  # plot these
  if (stringr::str_detect(stat, "count")) {
    stat_min <- -5
  } else if (stringr::str_detect(stat, "freq")) {
    stat_min <- -0.05
  }
  ggplot(
    plot_tbl |>
      dplyr::mutate(stat = pmax(stat_min, stat)),
    aes(x = man_stat, y = stat, col = gate_name)
  ) +
    cowplot::theme_cowplot(font_size = 20) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(aes(shape = batch_sh), size = 4, alpha = 0.7) +
    ggrepel::geom_text_repel(aes(label = stim),
      data = plot_tbl |>
        dplyr::group_by(batch) |>
        dplyr::filter(gate == min(gate)) |>
        dplyr::ungroup(),
      col = "gray70"
    ) +
    scale_color_manual(
      values = col,
      labels = gate_lab
    ) +
    scale_shape_manual(values = shape_vec) +
    geom_abline(intercept = 0, slope = 1) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
    coord_equal(xlim = axis_lim_vec, ylim = axis_lim_vec) +
    labs(
      x = paste0("Manual ", stat, " of ", marker), y = paste0("Auto ", stat, " of ", marker),
      title = pop_title
    ) +
    cowplot::background_grid(major = "xy", minor = "xy") +
    theme(legend.title = element_blank())
}

.plot_gate_stats_stim_vs_uns_dodge_gm <- function(gate_stats, pop, stat, col, gate_lab, marker, pop_title) {
  force(stat)
  ggplot(
    gate_stats |>
      dplyr::filter(.data$pop == .env$pop) |>
      dplyr::mutate(
        freq_bs = pmax(0, freq_bs),
        count_bs = pmax(0, count_bs)
      ),
    aes(x = stim, y = !!ensym(stat), fill = gate_name)
  ) +
    cowplot::theme_cowplot(font_size = 20) +
    geom_hline(yintercept = 0, size = 2) +
    geom_bar(position = "dodge", stat = "identity") +
    facet_wrap(~batch_sh, scales = "free_x") +
    cowplot::background_grid(major = "y") +
    labs(
      x = "Stimulation", y = paste0(stat, " of ", marker),
      title = pop_title
    ) +
    theme(legend.title = element_blank()) +
    theme(axis.text.x = element_text(size = 12)) +
    scale_fill_manual(
      values = col,
      labels = gate_lab
    )
}

.plot_gate_stats_stim_vs_uns_dodge_stim <- function(gate_stats, pop, stat, col, gate_lab, marker, pop_title) {
  force(stat)
  col <- RColorBrewer::brewer.pal(length(col), "RdBu") |>
    stats::setNames(names(col))

  if (stat %in% c("count", "freq")) {
    plot_tbl <- gate_stats
    colnames(plot_tbl) <- colnames(plot_tbl) |>
      stringr::str_replace(stat, "plot_stat")
    # calculate unstim for each
    plot_tbl <- plot_tbl |>
      dplyr::filter(.data$pop == .env$pop) |>
      dplyr::mutate(plot_stat_uns = plot_stat - plot_stat_bs)


    plot_tbl <- purrr::map_df(split(plot_tbl, plot_tbl$batch_sh), function(plot_tbl_batch) {
      plot_tbl_batch |>
        dplyr::filter(!stim == "uns") |>
        dplyr::select(batch_sh, gate_name, stim, plot_stat, plot_stat_uns) |>
        tidyr::pivot_longer(c(plot_stat, plot_stat_uns),
          names_to = "resp_type",
          values_to = "resp"
        )
    })

    p <- ggplot(
      plot_tbl,
      aes(x = gate_name, y = resp, alpha = resp_type, fill = gate_name)
    ) +
      cowplot::theme_cowplot(font_size = 16) +
      scale_alpha_manual(
        values = c("plot_stat" = 1, "plot_stat_uns" = 0),
        labels = c(
          "plot_stat" = "Stimulation (has colour)",
          "plot_stat_uns" = "Unstim (is blank)"
        ),
        name = "Stim or unstim"
      ) +
      geom_hline(yintercept = 0, size = 2) +
      facet_grid(batch_sh ~ stim, scales = "free") +
      geom_bar(position = "dodge", stat = "identity", size = 1.25, col = "black") +
      cowplot::background_grid(major = "y") +
      labs(
        x = "Gate method", y = paste0(stat, " of ", marker),
        title = pop_title
      ) +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_x_discrete(labels = gate_lab) +
      scale_fill_manual(
        values = col,
        labels = gate_lab,
        name = "Gating method"
      ) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    return(p)
  }

  if (stat %in% c("count_bs", "freq_bs")) {
    plot_tbl <- gate_stats
    colnames(plot_tbl) <- colnames(plot_tbl) |>
      stringr::str_replace(stat, "plot_stat")
    # calculate unstim for each
    plot_tbl <- plot_tbl |>
      dplyr::filter(
        .data$pop == .env$pop,
        stim != "uns"
      )

    p <- ggplot(
      plot_tbl |>
        dplyr::mutate(plot_stat = pmax(0, plot_stat)),
      aes(x = gate_name, y = plot_stat, fill = gate_name)
    ) +
      cowplot::theme_cowplot(font_size = 16) +
      geom_hline(yintercept = 0, size = 2) +
      facet_grid(batch_sh ~ stim, scales = "free") +
      geom_bar(position = "dodge", stat = "identity", size = 1.25, col = "black") +
      cowplot::background_grid(major = "y") +
      labs(
        x = "Gate method", y = paste0(stat, " of ", marker),
        title = pop_title
      ) +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_x_discrete(labels = gate_lab) +
      scale_fill_manual(
        values = col,
        labels = gate_lab,
        name = "Gating method"
      ) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    return(p)
  }
}

.plot_gate_stats_stim_vs_uns_jc <- function(gate_stats, pop, stat, col, gate_lab, marker, pop_title) {
  if (stat %in% c("count_bs", "freq_bs")) {
    return(ggplot())
  }

  if (stat %in% c("count", "freq")) {
    # extract pop of relevance and set negative counts and freqs to zero
    gate_stats <- gate_stats |>
      dplyr::filter(.data$pop == .env$pop) |>
      dplyr::mutate(
        freq_bs = pmax(0, freq_bs),
        count_bs = pmax(0, count_bs)
      ) |>
      dplyr::mutate(stim_facet = stim)

    plot_tbl_uns <- purrr::map_df(setdiff(gate_stats$stim, "uns"), function(stim) {
      gate_stats |>
        dplyr::filter(stim == "uns") |>
        dplyr::mutate(stim_facet = .env$stim)
    })

    plot_tbl <- gate_stats |>
      dplyr::filter(stim_facet != "uns") |>
      dplyr::bind_rows(plot_tbl_uns)

    col_name_vec <- colnames(plot_tbl)
    colnames(plot_tbl)[which(col_name_vec == stat)] <- "stat"

    plot_tbl <- plot_tbl |>
      dplyr::group_by(batch, batch_sh, stim_facet, gate_name) |>
      dplyr::arrange(batch, batch_sh, stim_facet, gate_name, stim) |>
      dplyr::summarise(jc = ifelse(sum(stat) == 0, 1, sum(stat[which(stim != "uns")] / sum(stat)))) |>
      dplyr::ungroup()

    p <- ggplot(plot_tbl, aes(x = gate_name, y = jc, fill = stim_facet)) +
      cowplot::theme_cowplot(font_size = 20) +
      geom_hline(yintercept = 0, size = 2) +
      # geom_boxplot() +
      geom_bar(position = "dodge", stat = "identity") +
      scale_fill_brewer(type = "div", palette = "PuOr") +
      cowplot::background_grid(major = "y") +
      labs(
        x = "Gate method", y = paste0(stat, " of ", marker),
        title = pop_title
      ) +
      theme(axis.text.x = element_text(size = 12, angle = 90)) +
      scale_x_discrete(labels = gate_lab) +
      facet_wrap(~batch_sh, ncol = 3) +
      theme(legend.title = element_blank())
    return(p)
  }

  ggplot()
}


# save plots
.save_plot_gate_stats_or_perf <- function(plot_list,
                                          params,
                                          plot_type,
                                          unlink = FALSE,
                                          path_project) {
  # create base directory if need be
  dir_base <- stim_gate_dir_base_create(params = params, dir_base_init = path_project)
  dir_save <- file.path(dir_base, plot_type)

  # if(dir.exists(dir_save)) unlink(dir_save, recursive = TRUE)
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)

  # save plots
  purrr::walk(names(plot_list), function(pop_save) {
    dir_save_pop <- file.path(dir_save, pop_save)
    if (dir.exists(dir_save_pop) && unlink) {
      purrr::walk(
        list.files(dir_save_pop, full.names = TRUE),
        file.remove
      )
    }

    if (!dir.exists(dir_save_pop)) dir.create(dir_save_pop)

    purrr::walk(names(plot_list[[pop_save]]), function(stat) {
      # print(stat)
      fp <- file.path(dir_save_pop, paste0(stat, ".png"))
      if (stringr::str_detect(fp, "dodge_stim") && (!(stringr::str_detect(fp, "count_bs") || stringr::str_detect(fp, "freq_bs")))) {
        p_height <- 15
        p_width <- 15
      } else {
        p_height <- 10.4
        p_width <- 15
      }
      if (stringr::str_detect(params$data_name, "gs_cytof_acs") && !stringr::str_detect(params$data_name, "test")) {
        if (!str_detect_any(stat, c("stim_vs_uns_dodge_stim-freq", "stim_vs_uns_dodge_stim-freq_bs"))) {
          return(invisible(TRUE))
        }
        p_height_mult <- length(params$data) / 30
        p_height <- 250
      }
      # print(fp)
      cowplot::ggsave2(fp, plot_list[[pop_save]][[stat]],
        height = p_height, width = p_width,
        limitsize = FALSE
      )
    })
  })
}
