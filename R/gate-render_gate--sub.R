#' @title Print plots in a directory
#' @export
.print_plots_in_dir <- function(dir_plot, print_header = TRUE, pattern = NULL, header_level = 5) {
  # get plot paths
  plot_path_vec <- list.files(dir_plot, full.names = TRUE)
  plot_path_vec <- plot_path_vec[!stringr::str_detect(plot_path_vec, "count")]
  plot_path_vec <- plot_path_vec[!stringr::str_detect(plot_path_vec, "dodge_gm")]
  # get plot anmes
  plot_name_vec <- list.files(dir_plot, full.names = FALSE) |>
    stringr::str_remove(".png")
  plot_name_vec <- plot_name_vec[!stringr::str_detect(plot_name_vec, "count")]
  plot_name_vec <- plot_name_vec[!stringr::str_detect(plot_name_vec, "dodge_gm")]
  # replace names with reasonable values
  plot_name_vec <- plot_name_vec |>
    stringr::str_replace("stim_vs_uns_jc", "Proportion of of stim over (stim + unstim) ") |>
    stringr::str_replace("auto_vs_man", "Auto vs manual") |>
    stringr::str_replace(
      "stim_vs_uns_dodge_stim-freq_bs",
      "Frequency (background subtracted)"
    ) |>
    stringr::str_replace(
      "stim_vs_uns_dodge_stim-freq",
      "Stim vs unstim - Frequency"
    ) |>
    stringr::str_replace("dodge_gm", " - gate method as dodge") |>
    stringr::str_replace("dodge_stim", " - stim as dodge") |>
    stringr::str_replace("-count_bs", "- count (background subtracted)") |>
    stringr::str_replace("-count", "- count") |>
    stringr::str_replace("-freq_bs", " frequency (background subtracted)") |>
    stringr::str_replace("-freq", "- frequency") |>
    stringr::str_replace("ccc", "Concordance correlation coefficient") |>
    stringr::str_replace("pcc", "Pearson's correlation coefficient") |>
    stringr::str_replace("std_diff_ind", "Standardised difference from manual") |>
    stringr::str_replace("std_diff_ind", "Absolute standardised difference from manual") |>
    stringr::str_replace("p_s2n_freq_bs", "Signal-to-noise for background-subtracted frequencies") |>
    stringr::str_replace("p_seprop_freq_bs", "Mean signal less one standard error of mean for background-subtracted frequencies") |>
    stringr::str_replace("p_s2n_ecdf", "Survival function for signal-to-noise ratio")
  # print plots for each batch
  for (i in seq_along(plot_name_vec)) {
    if (!is.null(pattern)) {
      if (!stringr::str_detect(plot_name_vec[i], pattern)) next
    }
    # print fcs name
    if (print_header) pander::pandoc.header(plot_name_vec[i], level = header_level)
    # print plot
    cat(paste0("![](", plot_path_vec[i], ")"), "\n")
  }
  invisible(TRUE)
}


#' @title Get performance stats
.get_results <- function(results,
                         boot,
                         .dir_base = dir_base,
                         stim = params$stim,
                         marker_lab = marker_lab_vec,
                         pop_lab = pop_lab_vec) {
  if (results == "stats") {
    if (boot) {
      tbl_name <- "gate_stats_tbl_boot"
    } else {
      tbl_name <- "gate_stats_tbl"
    }
  } else if (results == "perf") {
    tbl_name <- "gate_perf_tbl"
  } else {
    stop("results must be stats or perf")
  }
  dir_vec_pop <- list.dirs(.dir_base, recursive = FALSE)
  pop_vec <- list.dirs(.dir_base, recursive = FALSE, full.names = FALSE)
  purrr::map_df(seq_along(dir_vec_pop), function(i) {
    dir_vec_marker <- list.dirs(
      file.path(
        dir_vec_pop[i],
        stim
      ),
      recursive = FALSE
    )

    marker_vec <- list.dirs(
      file.path(
        dir_vec_pop[i],
        stim
      ),
      recursive = FALSE,
      full.names = FALSE
    )
    # print(pop_lab[pop_vec[i]])
    perf_list <- purrr::map(seq_along(dir_vec_marker), function(j) {
      # marker <- marker_vec[j]
      # print(marker_lab[marker])
      path_tbl <- file.path(
        dir_vec_marker[j], results,
        tbl_name
      )
      if (!file.exists(path_tbl)) {
        return(NULL)
      }
      readRDS(path_tbl) |>
        dplyr::mutate(marker = marker_lab[marker_vec[j]])
    }) |>
      purrr::compact() |>
      dplyr::bind_rows() |>
      tibble::as_tibble() |>
      dplyr::mutate(pop = pop_lab[pop_vec[i]])
  })
}

#' @title Plot a correlation plot for prepared data
#'
#' @param plot_tbl dataframe. A dataframe with columns x_factor, gate_name and score that are the discrete x- and y-axis variables and the fill variable, respectively.

.plot_heat_map_ind <- function(data, measure = "ccc") {
  if (measure %in% c("ccc", "pcc") || stringr::str_detect(measure, "prop")) {
    if (measure %in% c("ccc", "pcc")) col_vec <- RColorBrewer::brewer.pal(11, name = "RdBu") |> rev()
    if (stringr::str_detect(measure, "prop")) col_vec <- RColorBrewer::brewer.pal(11, name = "OrRd")[-1]
    ind_vec <- seq_along(col_vec)
    for (i in 1:2) ind_vec <- ind_vec[-length(ind_vec)]
    break_vec <- c(-1, 0, ind_vec / length(ind_vec))
    lim_vec <- c(0, 1)
  }

  if (stringr::str_detect(measure, "std_diff")) {
    col_vec <- RColorBrewer::brewer.pal(11, name = "RdBu") |> rev()
    col_vec <- col_vec[-c(1, length(col_vec))]
    break_vec <- c(-3, -2, 1, 0.5, 0, 0.5, 1, 2, 3)
    lim_vec <- c(-3, 3)
    data <- data |>
      dplyr::mutate(
        var_fill = pmax(-3, var_fill),
        var_fill = pmin(3, var_fill)
      )
  }


  ggplot(
    data,
    aes(x = var_x, y = var_y, fill = var_fill)
  ) +
    geom_raster(show.legend = FALSE) +
    geom_text(aes(label = var_fill_lab), col = "black") +
    scale_fill_gradientn(
      colors = col_vec,
      breaks = break_vec,
      limits = lim_vec
    ) +
    theme(legend.key.height = unit(12, units = "mm"))
}


#' @title Plot and save correlation plots
#' @export
plot_heat_map <- function(facet, var_fill,
                          data, measure = NULL,
                          stat = "freq_bs",
                          data_name = params$data_name,
                          gate_lab = gate_name_lab_vec,
                          marker_lab = marker_lab_vec,
                          pop_lab = pop_lab_vec,
                          measure_lab = measure_lab_vec,
                          stat_lab = stat_lab_vec,
                          stim_lab = stim_lab_vec,
                          stim = NULL,
                          extra_lab_col = NULL,
                          extra_lab_col_sep = "/") {
  # Preparation
  # --------------------------

  # order pops and markers for gs_proto
  if (str_detect_any(data_name, c(
    "gs_cytof", "gs_proto", "gs_cd8_base",
    "gs_cytof_acs"
  ))) {
    data <- data |>
      dplyr::mutate(
        pop = factor(.data$pop,
          levels = pop_lab_vec[c(
            4, 5, 6, 7,
            1, 2, 3, 8
          )]
        ),
        marker = factor(.data$marker,
          levels = marker_lab_vec[c(
            1, 7, 3, 2,
            4, 5, 6
          )]
        )
      )
  }
  # rename table for convenience
  if (facet == "marker") {
    data <- data |>
      dplyr::rename(
        facet = marker,
        var_x = pop,
        var_y = gate_name,
        var_fill = !!var_fill
      )
    title_x <- "Cell population"
    title_y <- "Gating method"
    lab_x_vec <- pop_lab
    lab_y_vec <- gate_lab
  } else if (facet == "pop") {
    data <- data |>
      dplyr::rename(
        facet = pop,
        var_x = marker,
        var_y = gate_name,
        var_fill = !!var_fill
      )
    title_x <- "Marker"
    title_y <- "Gating method"
    lab_x_vec <- marker_lab
    lab_y_vec <- gate_lab
  } else if (facet == "gate_name") {
    data <- data |>
      dplyr::rename(
        facet = gate_name,
        var_x = marker,
        var_y = pop,
        var_fill = !!var_fill
      )
    title_x <- "Marker"
    title_y <- "Cell population"
    lab_x_vec <- marker_lab
    lab_y_vec <- pop_lab
  }

  # Plot and save
  # -------------------------

  purrr::walk(unique(data$facet), function(facet_curr) {
    plot_tbl <- data |> dplyr::filter(facet == facet_curr)
    if (!is.null(measure)) plot_tbl <- plot_tbl |> dplyr::filter(measure == .env$measure)
    if (!is.null(stat)) plot_tbl <- plot_tbl |> dplyr::filter(stat == .env$stat)
    if (!is.null(stim)) plot_tbl <- plot_tbl |> dplyr::filter(stim == .env$stim)

    var_fill_lab_vec <- round(plot_tbl$var_fill, 2)
    if (!is.null(extra_lab_col)) {
      stringr::str_init <- paste0(
        var_fill_lab_vec, " (",
        plot_tbl[[extra_lab_col[1]]]
      )
      if (length(extra_lab_col) == 2) {
        stringr::str_end <- paste0("/", plot_tbl[[extra_lab_col[2]]], ")")
      } else {
        stringr::str_end <- ")"
      }
      var_fill_lab_vec <- paste0(stringr::str_init, stringr::str_end)
    }

    plot_tbl <- plot_tbl |>
      dplyr::mutate(var_fill_lab = var_fill_lab_vec)

    plot_tbl <- plot_tbl |>
      dplyr::select(var_x, var_y, var_fill, var_fill_lab)



    if (stringr::str_detect(facet, "gate_name")) facet_curr <- gate_lab[[facet_curr]]
    # if(facet == "gate_name") facet_curr <- gate_name_lab_vec[facet_curr]

    title <- facet_curr
    title_desc_vec <- c(
      stat,
      measure,
      stim
    )
    title_lab_list <- list()
    if (!is.null(stat)) {
      title_lab_list <- title_lab_list |> append(list(stat = stat_lab))
    } else {
      title_lab_list <- title_lab_list |> append(list(stat = NULL))
    }
    if (!is.null(measure)) {
      title_lab_list <- title_lab_list |> append(list(measure = measure_lab))
    } else {
      title_lab_list <- title_lab_list |> append(list(measure = NULL))
    }
    if (!is.null(stim)) {
      title_lab_list <- title_lab_list |> append(list(stim = stim_lab))
    } else {
      title_lab_list <- title_lab_list |> append(list(stim = NULL))
    }

    for (i in seq_along(title_lab_list)) {
      if (is.null(title_lab_list[[i]])) next
      title <- paste0(title, " - ", title_lab_list[[i]][title_desc_vec[i]])
    }


    p <- .plot_heat_map_ind(
      data = plot_tbl,
      measure = measure
    ) +
      labs(
        y = title_y, x = title_x,
        title = title
      ) +
      scale_y_discrete(labels = lab_y_vec) +
      scale_x_discrete(labels = lab_x_vec)

    dir_save <- file.path(dir_base, "overall", facet_curr)
    if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)

    dir_save_final <- dir_save
    inner_dirs <- ""
    if (is.null(stat)) {
      inner_dirs <- "stat-any"
    } else {
      inner_dirs <- stat_lab[stat]
    }
    if (is.null(measure)) {
      inner_dirs <- paste0(inner_dirs, "/measure-any")
    } else {
      inner_dirs <- paste0(inner_dirs, "/", measure_lab[measure])
    }
    dir_save_final <- file.path(dir_save, inner_dirs)
    if (is.null(stim)) {
      fn <- "All stimulations.png"
    } else {
      fn <- paste0(stim_lab[stim], ".png")
    }

    if (!dir.exists(dir_save_final)) dir.create(dir_save_final, recursive = TRUE)
    cowplot::ggsave2(file.path(dir_save_final, fn), p,
      height = 20, width = 37, units = "cm"
    )
  })

  invisible(TRUE)
}

#' @title Create a proportion above a cutoff plot
#' @export
.plot_prop_plot <- function(.s2n_tbl = s2n_tbl, level = 1, gate_name = "locb5_no") {
  col_vec <- RColorBrewer::brewer.pal(11, name = "RdBu") |> rev()
  break_vec <- c(-1, 0, (1:9) / 9)

  plot_tbl <- .s2n_tbl |>
    dplyr::filter(gate_name == "locb5_no") |>
    dplyr::filter(!is.na(freq_bs_s2n)) |>
    dplyr::group_by(
      pop, marker, gate_name, gate_type,
      gate_combn
    ) |>
    dplyr::summarise(prop_g = sum(freq_bs_s2n > level) / dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::arrange(gate_name, pop, desc(prop_g))

  ggplot(
    plot_tbl,
    aes(x = marker, y = pop, fill = prop_g)
  ) +
    geom_raster(show.legend = FALSE) +
    geom_text(aes(label = round(prop_g, 2)), col = "black") +
    scale_fill_gradientn(
      colors = col_vec,
      breaks = break_vec,
      limits = c(0, 1)
    ) +
    theme(legend.key.height = unit(12, units = "mm"))
}

# get freq_bs for non-bootstrap stats for the selected gate and man_no
.plot_auto_vs_man_perf <- function(marker, stats = .get_results("stats", FALSE),
                                   stat = "freq_bs", gate_name = params$gate_name,
                                   dir_base_init = dir_base,
                                   gate_lab = gate_name_lab_vec,
                                   stat_lab = stat_lab_vec,
                                   marker_lab = marker_lab_vec) {
  plot_tbl <- stats |>
    dplyr::filter(.data$gate_name %in% c(.env$gate_name, c("man_no"))) |>
    dplyr::filter(.data$marker == .env$marker)

  # adjust column name for programming sake
  col_name_vec <- colnames(plot_tbl)
  stat_ind <- which(col_name_vec == stat)
  col_name_vec[stat_ind] <- "stat"
  colnames(plot_tbl) <- col_name_vec

  # remove NA entries
  plot_tbl <- plot_tbl |> dplyr::filter(!is.na(stat))

  # get manual cut values
  man_stat_temp_tbl <- plot_tbl |> dplyr::filter(gate_name == "man_no")
  man_stat_vec <- stats::setNames(
    man_stat_temp_tbl[["stat"]],
    paste0(
      man_stat_temp_tbl[["fcs"]], "_",
      man_stat_temp_tbl[["pop"]]
    )
  )

  # get plot tbl
  plot_tbl <- plot_tbl |>
    dplyr::filter(!gate_name == "man_no") |>
    dplyr::mutate(man_stat = man_stat_vec[paste0(fcs, "_", pop)])

  # get shape label vector
  batch_vec <- unique(plot_tbl$batch_sh)
  shape_vec <- stats::setNames(
    c(15, 16, 17, 18, 14, 10, 6, 7, 8, 12, 23)[1:length(batch_vec)],
    batch_vec
  )

  # plot these
  if (stringr::str_detect(stat, "count")) {
    stat_min <- -5
  } else if (stringr::str_detect(stat, "freq")) {
    stat_min <- -0.05
  }
  p <- ggplot(
    plot_tbl |>
      dplyr::mutate(
        stat = pmax(stat_min, stat),
        man_stat = pmax(stat_min, man_stat)
      ),
    aes(x = man_stat, y = stat)
  ) +
    cowplot::theme_cowplot(font_size = 14) +
    cowplot::background_grid(major = "xy") +
    geom_abline(slope = 1, intercept = 0) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 0.1, col = "red", linetype = "dotted", size = 1) +
    geom_vline(xintercept = 0.1, col = "red", linetype = "dotted", size = 1) +
    geom_vline(xintercept = 0) +
    geom_smooth(se = FALSE, span = 1.5) +
    geom_point(aes(shape = batch_sh, col = stim), size = 4, alpha = 0.4) +
    facet_wrap(~pop, ncol = 3, scales = "free") +
    labs(
      x = paste0("Manual ", stringr::str_to_lower(stat_lab[stat])),
      y = paste0("Auto ", stringr::str_to_lower(stat_lab[stat])),
      title = marker
    )

  dir_save <- file.path(
    dir_base_init, "overall", gate_lab[[gate_name]],
    stat_lab[[stat]], "Auto vs manual"
  )
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  fn <- paste0(marker, ".png")
  cowplot::ggsave2(file.path(dir_save, fn), p,
    height = 25, width = 25, units = "cm"
  )
}
