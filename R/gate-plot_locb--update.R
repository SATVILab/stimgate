.update_and_save_locb_plots_with_lines <- function(params, gate_tbl, col, gate_lab) {
  purrr::walk(params$bias_uns, function(bias) {
    plot_list <- .get_locb_plots_with_lines(
      bias = bias, params = params,
      gate_tbl = gate_tbl, col = col,
      gate_lab = gate_lab
    )

    .save_plot_gate(params = params, plot_list = plot_list)
  })
}

.get_locb_plots_with_lines <- function(bias, params, gate_tbl, col, gate_lab) {
  # Load plots
  # ----------------------------------------

  dir_load <- file.path(tempdir(), params$data_name, paste0("cp_locb", bias, "_plots"))
  fn_vec <- list.files(dir_load, full.names = TRUE)

  p_list <- list("gate" = purrr::map(fn_vec, function(fn) readRDS(fn)) |> flatten())

  p_list[[1]] <- .add_cp_lines_to_plot_batch(
    p_list = p_list[[1]],
    gate_tbl = gate_tbl,
    col = col,
    gate_lab = gate_lab
  )

  p_list
}

.add_cp_lines_to_plot_batch <- function(p_list, gate_tbl, col, gate_lab) {
  ind_list <- stringr::str_split(names(p_list), "_")
  purrr::map(seq_along(p_list), function(i) { # loop over batches
    .add_cp_lines_to_plot_sample(
      p_list = p_list[[i]], # plots in batch
      gate = gate_tbl,
      ind = ind_list[[i]],
      col = col,
      gate_lab = gate_lab
    )
  }) |>
    stats::setNames(names(p_list))
}


#' @title Add gate lines to individual samples
#'
#' @param p_list named list. Named list where names are names of samples
#' and the values are themselves named lists of plots.
#' @param gate dataframe
.add_cp_lines_to_plot_sample <- function(p_list, gate, ind, col, gate_lab) {
  purrr::map(seq_along(p_list), function(i) { # loop over plots
    .add_cp_lines_to_plot_plot(
      p_list = p_list[[i]], # plots for just this sample
      gate = gate <- gate |> dplyr::filter(.data$ind == .env$ind[i]), # gate_tbl for just this sample
      col = col,
      gate_lab = gate_lab
    )
  }) |>
    stats::setNames(names(p_list))
}


#' @title Add gate lines to multiple plots
#'
#' @param p_list named list. Named list where each name is the name of plot and
#' each element is a plot.
#' @param gate dataframe. Dataframe with columns gate_name and gate, giving
#' the name of the gate and its value, respectively. All rows in this
#' dataframe will be plotted on each plot.
#' @param col named character vector. Colours of gate lines.
#' @param gate_lab named character vector. Specifies names of gates.
.add_cp_lines_to_plot_plot <- function(p_list, gate, col, gate_lab) {
  purrr::map(p_list, function(p) {
    .add_cp_lines_to_plot_ind(p,
      gate = gate,
      col = col, gate_lab = gate_lab
    )
  }) |>
    stats::setNames(names(p_list))
}


#' @title Add gate lines to a single plot
#'
#' @param p ggplot object. Plot to add gates to.
.add_cp_lines_to_plot_ind <- function(p, gate, col, gate_lab) {
  p +
    geom_vline(
      aes(
        xintercept = gate,
        col = gate_name
      ),
      data = gate, size = 1.5
    ) +
    scale_colour_manual(
      values = col,
      labels = gate_lab
    )
}
