.get_gate_perf <- function(gate_stats, set_to_zero = TRUE, params, path_project) {
  # ================================
  # Preparation
  # ================================\

  if (set_to_zero) {
    gate_stats <- gate_stats |>
      dplyr::mutate(
        count_bs = pmax(0, count_bs),
        freq_bs = pmax(0, freq_bs)
      )
  }

  # get manual cut values
  man_stat_temp_tbl <- gate_stats |> dplyr::filter(gate_name == "man_no")

  # merged tbl
  join_tbl <- gate_stats |>
    dplyr::filter(gate_name != "man_no") |>
    dplyr::left_join(man_stat_temp_tbl, by = c("fcs", "pop"))

  # ================================
  # Calculate measure scores
  # ================================

  gate_perf_tbl <- join_tbl |>
    dplyr::group_by(pop, gate_name.x) |>
    dplyr::summarise(
      pcc.freq = .calc_pcc(freq.x, freq.y),
      pcc.freq_bs = .calc_pcc(freq_bs.x, freq_bs.y),
      pcc.count = .calc_pcc(count.x, count.y),
      pcc.count_bs = .calc_pcc(count_bs.x, count_bs.y),
      ccc.freq = .calc_ccc(freq.x, freq.y),
      ccc.freq_bs = .calc_ccc(freq_bs.x, freq_bs.y),
      ccc.count = .calc_ccc(count.x, count.y),
      ccc.count_bs = .calc_ccc(count_bs.x, count_bs.y),
      std_diff_mean.freq_bs =
        .get_std_diff(
          freq_bs_auto = freq_bs.x,
          freq_bs_man = freq_bs.y,
          n_stim = n_cell_stim.y, n_uns = n_cell_uns.y,
          n_stim_pos_man = count.y,
          n_uns_pos_man = count.y - count_bs.y
        ) |>
          abs() |>
          mean(na.rm = TRUE),
      std_diff_mean_trim10.freq_bs =
        .get_std_diff(
          freq_bs_auto = freq_bs.x, freq_bs_man = freq_bs.y,
          n_stim = n_cell_stim.y, n_uns = n_cell_uns.y,
          n_stim_pos_man = count.y, n_uns_pos_man = count.y - count_bs.y
        ) |>
          abs() |>
          mean(trim = 0.1, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(
      pcc.freq:std_diff_mean_trim10.freq_bs,
      names_to = "measure",
      values_to = "score"
    ) |>
    tidyr::separate(
      col = "measure", into = c("measure", "stat"),
      sep = "\\."
    ) |>
    dplyr::rename(gate_name = gate_name.x) |>
    dplyr::mutate(level = "method", ind = NA, batch = NA, stim = NA)

  gate_perf_tbl_ind <- join_tbl |>
    dplyr::group_by(pop, gate_name.x) |>
    dplyr::mutate(
      std_diff_ind.freq_bs = .get_std_diff(
        freq_bs_auto = freq_bs.x, freq_bs_man = freq_bs.y,
        n_stim = n_cell_stim.y, n_uns = n_cell_uns.y,
        n_stim_pos_man = count.y, n_uns_pos_man = count.y - count_bs.y
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!stim.x == "uns") |>
    dplyr::select(
      pop, gate_name.x, ind.x, batch.x, stim.x, std_diff_ind.freq_bs
    ) |>
    tidyr::pivot_longer(
      std_diff_ind.freq_bs,
      names_to = "measure", values_to = "score"
    ) |>
    tidyr::separate(
      col = "measure", into = c("measure", "stat"), sep = "\\."
    ) |>
    dplyr::rename(
      gate_name = gate_name.x, ind = ind.x, batch = batch.x, stim = stim.x
    ) |>
    dplyr::mutate(level = "ind")

  gate_perf_tbl <- gate_perf_tbl |> dplyr::bind_rows(gate_perf_tbl_ind)

  # ==========================
  # Save table
  # ==========================

  # create base directory if need be
  dir_base <- stim_gate_dir_base_create(params = params, dir_base_init = path_project)
  dir_save <- file.path(dir_base, "perf")
  # if(dir.exists(dir_save)) unlink(dir_save, recursive = FALSE)
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = FALSE)

  # save stats tbl
  saveRDS(gate_perf_tbl, file = file.path(dir_save, "gate_perf_tbl"))
  write.csv(gate_perf_tbl, file = file.path(dir_save, "gate_perf_tbl.csv"))

  gate_perf_tbl
}

#' @title Calculate standardised difference to manual
.get_std_diff_ind <- function(freq_bs_auto,
                              freq_bs_man,
                              n_stim,
                              n_uns,
                              n_stim_pos_man,
                              n_uns_pos_man) {
  diff_prop <- (freq_bs_auto - freq_bs_man) / 100
  n_stim_pos_man <- max(1, n_stim_pos_man)
  prop_stim_man <- n_stim_pos_man / n_stim
  var_stim_man <- prop_stim_man * (1 - prop_stim_man) / n_stim
  n_uns_pos_man <- max(1, n_uns_pos_man)
  prop_uns_man <- n_uns_pos_man / n_uns
  var_uns <- prop_uns_man * (1 - prop_uns_man) / n_uns
  diff_prop / sqrt(var_stim_man + var_uns)
}


.get_std_diff <- function(freq_bs_auto,
                          freq_bs_man,
                          n_stim,
                          n_uns,
                          n_stim_pos_man,
                          n_uns_pos_man) {
  purrr::map_dbl(seq_along(freq_bs_auto), function(i) {
    .get_std_diff_ind(
      freq_bs_auto = freq_bs_auto[i],
      freq_bs_man = freq_bs_man[i],
      n_stim = n_stim[i],
      n_uns = n_uns[i],
      n_stim_pos_man = n_stim_pos_man[i],
      n_uns_pos_man[i]
    )
  })
}

#' @title Plot PCC and CCC for each pop, stat and cutpoint method
.plot_gate_perf <- function(gate_perf, gate_lab, params, path_project) {
  # ================================
  # Preparation
  # ================================

  # col_vec <- RColorBrewer::brewer.pal(9, name = "OrRd")[-9]
  col_vec <- RColorBrewer::brewer.pal(11, name = "RdBu") |> rev()

  break_vec <- c(-1, 0, (1:9) / 9)
  # ==========================
  # Get plots
  # ==========================


  plot_list <- purrr::map(unique(gate_perf$pop), function(pop) {
    title <- ifelse(pop == "gate", params$pop_gate, pop)
    purrr::map(unique(gate_perf$measure), function(measure) {
      if (measure %in% c("ccc", "pcc")) {
        p <- ggplot(
          gate_perf |> dplyr::filter(
            .data$pop == .env$pop,
            .data$measure == .env$measure
          ),
          aes(x = stat, y = gate_name, fill = score)
        ) +
          geom_raster(show.legend = FALSE) +
          geom_text(aes(label = round(score, 2)), col = "black") +
          scale_fill_gradientn(
            colors = col_vec,
            breaks = break_vec,
            limits = c(0, 1)
          ) +
          theme(legend.key.height = unit(12, units = "mm")) +
          scale_y_discrete(labels = gate_lab) +
          labs(
            y = "Cutpoint method", x = "Statistic",
            title = title
          ) +
          scale_x_discrete(
            labels = c(
              "freq_bs" = "Freq (BS)",
              "freq" = "Freq",
              "count" = "Count",
              "count_bs" = "Count (BS)"
            )
          )
      } else if (stringr::str_detect(measure, "std_diff_mean")) {
        return(NULL)
        p <- ggplot(
          gate_perf |> dplyr::filter(
            .data$pop == .env$pop,
            .data$measure == .env$measure
          ),
          aes(x = stat, y = gate_name, fill = score)
        ) +
          geom_raster(show.legend = FALSE) +
          geom_text(aes(label = round(score, 2)), col = "black") +
          theme(legend.key.height = unit(12, units = "mm")) +
          scale_y_discrete(labels = gate_lab) +
          labs(
            y = "Cutpoint method", x = "Statistic",
            title = title
          ) +
          scale_x_discrete(
            labels = c(
              "freq_bs" = "Freq (BS)",
              "freq" = "Freq",
              "count" = "Count",
              "count_bs" = "Count (BS)"
            )
          )
      } else if (stringr::str_detect(measure, "std_diff_ind")) {
        plot_tbl <- gate_perf |>
          dplyr::filter(
            .data$pop == .env$pop,
            .data$measure == .env$measure
          ) |>
          dplyr::mutate(batch = stringr::str_sub(batch, end = 6))
        batch_vec <- unique(plot_tbl$batch)
        shape_vec <- stats::setNames(
          c(15, 16, 17, 18, 14, 10, 6, 7, 8, 12, 23)[seq_along(batch_vec)],
          batch_vec
        )
        p <- ggplot(
          plot_tbl,
          aes(x = gate_name, y = score)
        ) +
          geom_hline(yintercept = 0, size = 1, linetype = "dotted", col = "gray50") +
          geom_jitter(aes(col = stim, shape = batch), width = 0.25, size = 3) +
          cowplot::background_grid(major = "y", minor = "y") +
          geom_boxplot(alpha = 0.1) +
          scale_shape_manual(
            values = shape_vec,
            name = "Batch"
          ) +
          scale_colour_brewer(
            palette = "Set1",
            name = "Stimulation"
          ) +
          labs(
            x = "Gate method",
            y = "Standardised difference from manual in\n background-subtracted frequency",
            title = title
          ) +
          theme(legend.key.height = unit(12, units = "mm")) +
          scale_x_discrete(labels = gate_lab)
      }
      p
    }) |>
      stats::setNames(unique(gate_perf$measure)) |>
      purrr::compact()
  }) |>
    stats::setNames(unique(gate_perf$pop))

  plot_list_add <- purrr::map(unique(gate_perf$pop), function(pop) {
    title <- ifelse(pop == "gate", params$pop_gate, pop)
    purrr::map(unique(gate_perf$measure), function(measure) {
      if (stringr::str_detect(measure, "std_diff_ind")) {
        plot_tbl <- gate_perf |>
          dplyr::filter(
            .data$pop == .env$pop,
            .data$measure == .env$measure
          ) |>
          dplyr::mutate(batch = stringr::str_sub(batch, end = 6))
        batch_vec <- unique(plot_tbl$batch)
        shape_vec <- stats::setNames(
          c(15, 16, 17, 18, 14, 10, 6, 7, 8, 12, 23)[seq_along(batch_vec)],
          batch_vec
        )
        p <- ggplot(
          plot_tbl,
          aes(x = gate_name, y = abs(score))
        ) +
          geom_hline(
            yintercept = 0, size = 1, linetype = "dotted", col = "gray50"
          ) +
          geom_jitter(
            aes(col = stim, shape = batch),
            width = 0.25, size = 3
          ) +
          cowplot::background_grid(major = "y", minor = "y") +
          geom_boxplot(alpha = 0.1) +
          scale_shape_manual(
            values = shape_vec,
            name = "Batch"
          ) +
          scale_colour_brewer(
            palette = "Set1",
            name = "Stimulation"
          ) +
          labs(
            x = "Gate method",
            y = "Absolute standardised difference from manual in\n background-subtracted frequency", # nolint
            title = title
          ) +
          theme(legend.key.height = unit(12, units = "mm")) +
          scale_x_discrete(labels = gate_lab)
      } else {
        return(NULL)
      }
      p
    }) |>
      stats::setNames(paste0(unique(gate_perf$measure), "-abs")) |>
      purrr::compact()
  }) |>
    stats::setNames(unique(gate_perf$pop))

  plot_list <- plot_list |> append(plot_list_add)

  plot_list <- purrr::map(unique(names(plot_list)), function(x) {
    plot_list[[x]] |>
      append(plot_list_add[[x]])
  }) |>
    stats::setNames(unique(names(plot_list)))
  # ==========================
  # Save plots
  # ==========================

  print("saving perf plots")
  .save_plot_gate_stats_or_perf(
    plot_list = plot_list,
    params = params,
    plot_type = "perf",
    path_project = path_project
  )


  invisible(TRUE)
}

#' @title Calculate CCC
#'
#' @description
#' A wrapper around \code{cccrm::cccUst}.
#'
#' @param x,y numeric vectors of equal length. Measurements.
#'
#' @return Numeric.
.calc_ccc <- function(x, y) {
  x_vec <- x[!is.na(x)]
  y_vec <- y[!is.na(y)]
  if (sd(x_vec) == 0 || sd(y_vec) == 0) {
    return(1.1)
  }
  data_tbl <- tibble::tibble(x = x, y = y) |>
    dplyr::filter(!is.na(x), !is.na(y)) |>
    tidyr::pivot_longer(
      cols = c("x", "y"),
      names_to = "rmet",
      values_to = "ry"
    )

  suppressWarnings(cccrm::cccUst(
    dataset = data_tbl,
    ry = "ry",
    rmet = "rmet"
  )[["CCC"]])
}

#' @title Calculate PCC
#'
#' @description
#' A wrapper around \code{stats::cor}.
#'
#' @param x,y numeric vectors of equal length. Measurements.
#'
#' @return Numeric.
.calc_pcc <- function(x, y) {
  x_vec <- x[!is.na(x)]
  y_vec <- y[!is.na(y)]
  if (sd(x_vec) == 0 || sd(y_vec) == 0) {
    return(1.1)
  }
  cor(x[!is.na(x)], y[!is.na(y)])
}

.get_s2n <- function(boot, est = NULL, rep_disp = 0.5, type = "mean") {
  if (type == "median") {
    if (is.null(est)) est <- median(boot, na.rm = TRUE)
    stat_mad <- mad(stat, na.rm = TRUE)
    if (est == 0) {
      return(0)
    }
    if (stat_mad == 0) stat_mad <- rep_disp
    return(est / stat_mad)
  } else if (type == "mix") {
    if (is.null(est)) est <- median(stat, na.rm = TRUE)
    stat_se <- sd(boot, na.rm = TRUE)
    if (est == 0) {
      return(0)
    }
    if (is.na(stat_se)) {
      return(NA)
    }
    return(est / stat_se)
  } else if (type == "mean") {
    if (is.null(est)) est <- mean(stat, na.rm = TRUE)
    stat_se <- sd(boot, na.rm = TRUE)
    if (est == 0) {
      return(0)
    }
    if (is.na(stat_se)) stat_se <- rep_disp
    return(est / stat_se)
  }
}
