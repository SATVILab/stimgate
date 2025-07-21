#' Plot stimulation gate
#'
#' Plot bivariate hex and univariate density plots for batches of samples, along
#' with their gates.
#'
#' @param ind numeric vector. Specifies indices in `.data` to plot.
#' @param .data GatingSet. Same GatingSet passed to `stimgate_gate`.
#' @param path_project character. Path to the project directory used for `stimgate_gate`.
#' @param marker character vector of length one or two. Specifies markers (channels, really)
#' to be plotted. If only one is passed, then only univariate plots are created.
#' @param ind_lab named character vector.
#' Labels for `ind` used in plot.
#' Optional.
#' @param marker_lab named character vector.
#' Labels for `marker` used in plot.
#' Optional.
#' @param exc_min Logical.
#' If `TRUE`, excludes the minimum expression values when processing the data.
#' Default is `TRUE`.
#' @param limits_expand list.
#' Expand the limits of the plot axes.
#' Default is `NULL`.
#' @param limits_equal Logical.
#' If TRUE, forces equal lengths of the limits.
#' @param grid Logical.
#' If TRUE, arranges the resulting plots in a grid format
#' using `cowplot::plot_grid`.
#' Default is `TRUE`.
#' @param grid_n_col Integer.
#' Number of columns in grid layout.
#' @param show_gate Logical.
#' If `TRUE`, overlays gate lines on the plots.|>
#' Default is `TRUE`.
#' @param min_cell integer.
#' Minimum number of cells to be plotted.
#' Will skip plots with fewer cells.
#' Default is 10.
#'
#' @return A grid of plots if `grid` is TRUE, otherwise a list of ggplot objects.
#'
#' @examples
#' # Create example data and run gating
#' example_data <- get_example_data()
#' gs <- flowWorkspace::load_gs(example_data$path_gs)
#' path_project <- file.path(dirname(example_data$path_gs), "stimgate")
#'
#' # Run gating
#' stimgate::stimgate_gate(
#'   .data = gs,
#'   path_project = path_project,
#'   pop_gate = "root",
#'   batch_list = example_data$batch_list,
#'   marker = example_data$marker
#' )
#'
#' # Create plots
#' plots <- stimgate_plot(
#'   ind = example_data$batch_list[[1]], # indices in `gs` to plot
#'   .data = gs, # GatingSet
#'   path_project = path_project,
#'   marker = example_data$marker,
#'   grid = TRUE
#' )
#' @export
stimgate_plot <- function(ind,
                          .data,
                          path_project,
                          marker,
                          ind_lab = NULL,
                          marker_lab = NULL,
                          exc_min = TRUE,
                          limits_expand = NULL,
                          limits_equal = FALSE,
                          grid = TRUE,
                          grid_n_col = 2,
                          show_gate = TRUE,
                          min_cell = 10) {
  p_list <- .plot_gate(
    ind = ind, ind_lab = ind_lab, .data = .data,
    marker = marker, marker_lab = marker_lab,
    path_project = path_project, exc_min = exc_min,
    limits_expand = limits_expand, limits_equal = limits_equal,
    show_gate = show_gate, min_cell = min_cell
  )
  if (length(p_list) == 0L) {
    return(NULL)
  }
  .plot_grid(plot = grid, p_list = p_list, n_col = grid_n_col)
}

.plot_gate <- function(marker,
                       ind,
                       ind_lab,
                       .data,
                       marker_lab,
                       path_project,
                       exc_min,
                       limits_expand,
                       limits_equal,
                       show_gate,
                       min_cell) {
  # bv
  p_list_bv <- .plot_gate_bv(
    marker = marker, ind = ind, ind_lab = ind_lab,
    .data = .data, marker_lab = marker_lab,
    path_project = path_project, exc_min = exc_min,
    limits_expand = limits_expand, limits_equal = limits_equal,
    show_gate = show_gate, min_cell = min_cell
  )

  # uv
  p_list_uv <- .plot_gate_uv(
    ind = ind, ind_lab = ind_lab, .data = .data,
    marker = marker, exc_min = exc_min, marker_lab = marker_lab,
    show_gate = show_gate, path_project = path_project,
    min_cell = min_cell
  )

  p_list_bv |> append(p_list_uv)
}

.plot_gate_bv <- function(marker,
                          ind,
                          ind_lab,
                          .data,
                          marker_lab,
                          path_project,
                          exc_min,
                          limits_expand,
                          limits_equal,
                          show_gate,
                          min_cell) {
  if (length(marker) == 1L) {
    return(NULL)
  }
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    install.packages("hexbin")
  }
  p_list <- lapply(seq_along(ind), function(i) {
    ind_curr <- ind[[i]]
    ex_tbl <- .plot_get_ex_tbl(ind_curr, .data, marker, exc_min = FALSE)
    if (nrow(ex_tbl) < min_cell) {
      return(NULL)
    }
    p <- plot_cyto(
      data = ex_tbl,
      marker = marker,
      exc_min = exc_min,
      limits_expand = limits_expand,
      limits_equal = limits_equal
    ) +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")
      )
    p <- .plot_add_axis_title(p, marker, marker_lab)
    p <- .plot_add_title(p, ind_curr, i, ind_lab)
    p <- .plot_add_gate(p, .data, ind_curr, marker, path_project, show_gate)
    p
  }) |>
    stats::setNames(.plot_get_lab(ind, ind_lab))
  p_list <- p_list[vapply(p_list, Negate(is.null), logical(1))]
  if (length(p_list) == 0L) {
    return(NULL)
  }
  p_list
}

.plot_get_ex_tbl <- function(ind, .data, marker, exc_min) {
  lapply(ind, function(x) {
    .plot_get_ex_tbl_ind(x, .data, marker, exc_min) |>
      dplyr::mutate(ind = x)
  }) |>
    Reduce(rbind, x = _)
}


.plot_get_ex_tbl_ind <- function(ind, .data, marker, exc_min) {
  fr <- flowWorkspace::gh_pop_get_data(.data[[ind]])
  ex_tbl <- flowCore::exprs(fr) |> tibble::as_tibble()
  ex_tbl <- ex_tbl[, marker, drop = FALSE]
  ex_tbl <- .plot_get_ex_tbl_ind_exc_min(ex_tbl, exc_min, marker)
  ex_tbl
}

.plot_get_ex_tbl_ind_exc_min <- function(ex_tbl, exc_min, marker) {
  if (!exc_min) {
    return(ex_tbl)
  }
  n_row_init <- nrow(ex_tbl)
  attr(ex_tbl, "n_row_init") <- n_row_init
  min_val_vec <- vapply(
    marker, function(x) min(ex_tbl[[x]], na.rm = TRUE), numeric(1)
  )
  for (i in seq_along(marker)) {
    ex_tbl <- ex_tbl[ex_tbl[[marker[i]]] > min_val_vec[i], ]
  }
  n_row_final <- nrow(ex_tbl)
  attr(ex_tbl, "prob_g_min") <- n_row_final / n_row_init

  ex_tbl
}

.plot_add_axis_title <- function(p, val, val_lab) {
  lab <- .plot_get_lab(val, val_lab)
  p <- p + labs(x = lab[[1]])
  if (length(lab) > 1L) {
    p <- p + labs(y = lab[[2]])
  }
  p
}

.plot_add_title <- function(p, ind, i, ind_lab) {
  p + ggtitle(.plot_get_lab(ind, ind_lab, i))
}

.plot_get_lab <- function(val, val_lab, i = NULL) {
  if (is.null(val_lab)) {
    return(val)
  }
  lab <- if (!is.null(names(val_lab))) {
    val_lab[val]
  } else {
    if (!is.null(i)) val_lab[i] else val_lab
  }
  lab |> stats::setNames(NULL)
}

.plot_add_gate <- function(p, .data, ind, marker, path_project, show_gate) {
  if (!show_gate) {
    return(p)
  }
  gate_tbl <- .plot_get_gate_tbl(ind, marker, path_project)
  for (i in seq_along(marker)) {
    gate_vec <- gate_tbl[["gate"]][gate_tbl[["chnl"]] == marker[i]]
    for (j in seq_along(gate_vec)) {
      gate <- gate_vec[[j]]
      if (i == 1) {
        p <- p +
          geom_vline(
            xintercept = gate,
            color = "red",
            alpha = 0.5
          ) +
          expand_limits(
            x = gate * 1.1
          )
      } else if (i == 2) {
        p <- p +
          geom_hline(
            yintercept = gate,
            color = "red",
            alpha = 0.5
          ) +
          expand_limits(y = gate * 1.1)
      }
    }
  }
  p
}

.plot_get_gate_tbl <- function(ind, marker, path_project) {
  gate_tbl <- stimgate_gate_get(path_project) |>
    dplyr::group_by(gate_name, chnl, marker, ind, batch) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  gate_tbl <- gate_tbl[gate_tbl[["ind"]] %in% ind, ]
  ind_vec <- NULL
  for (i in seq_along(marker)) {
    ind_vec[[i]] <- which(gate_tbl[["chnl"]] == marker[i])
  }
  ind_vec <- ind_vec |> unlist()
  gate_tbl[ind_vec, ]
}

.plot_gate_uv <- function(ind,
                          ind_lab,
                          .data,
                          marker,
                          exc_min,
                          marker_lab,
                          show_gate,
                          path_project,
                          min_cell) {
  p_list <- lapply(marker, function(m) {
    .plot_gate_uv_marker(
      marker = m, ind = ind, .data = .data,
      exc_min = exc_min, ind_lab = ind_lab,
      marker_lab = marker_lab, show_gate = show_gate,
      path_project = path_project, min_cell = min_cell
    )
  }) |>
    stats::setNames(.plot_get_lab(marker, marker_lab))
  p_list <- p_list[vapply(p_list, Negate(is.null), logical(1))]
  if (length(p_list) == 0L) {
    return(NULL)
  }
  p_list
}

.plot_gate_uv_marker <- function(marker,
                                 ind,
                                 .data,
                                 exc_min,
                                 ind_lab,
                                 marker_lab,
                                 show_gate,
                                 path_project,
                                 min_cell) {
  plot_tbl <- .plot_gate_uv_marker_get_plot_tbl(
    ind = ind, .data = .data, marker = marker, exc_min = exc_min,
    ind_lab = ind_lab, min_cell = min_cell
  )
  if (is.null(plot_tbl)) {
    return(NULL)
  }
  .plot_gate_uv_marker_plot(
    plot_tbl = plot_tbl, exc_min = exc_min,
    ind = ind, ind_lab = ind_lab,
    marker = marker, marker_lab = marker_lab,
    show_gate = show_gate, path_project = path_project
  )
}

.plot_gate_uv_marker_get_plot_tbl <- function(ind,
                                              .data,
                                              marker,
                                              exc_min,
                                              ind_lab,
                                              min_cell) {
  plot_tbl_list <- lapply(seq_along(ind), function(i) {
    plot_tbl <- .plot_gate_uv_marker_get_plot_tbl_ind(
      ind = ind[[i]], .data = .data, marker = marker, exc_min = exc_min,
      min_cell = min_cell
    )
    if (is.null(plot_tbl)) {
      return(NULL)
    }
    plot_tbl[, "ind"] <- ind[[i]]
    plot_tbl[, "ind_lab"] <- .plot_get_lab(ind[[i]], ind_lab, i)
    plot_tbl
  })
  plot_tbl_list <- plot_tbl_list[
    vapply(plot_tbl_list, Negate(is.null), logical(1))
  ]
  if (length(plot_tbl_list) == 0L) {
    return(NULL)
  }
  Reduce(rbind, plot_tbl_list)
}

.plot_gate_uv_marker_get_plot_tbl_ind <- function(ind,
                                                  .data,
                                                  marker,
                                                  exc_min,
                                                  min_cell) {
  ex_tbl <- .plot_get_ex_tbl(ind, .data, marker, exc_min)
  if (nrow(ex_tbl) < min_cell) {
    return(NULL)
  }
  dens_obj_raw <- density(ex_tbl[[marker]], na.rm = TRUE)
  plot_tbl <- tibble::tibble(x = dens_obj_raw$x, y = dens_obj_raw$y)
  .plot_gate_uv_marker_add_adj(
    exc_min, plot_tbl, dens_obj_raw, attr(ex_tbl, "prob_g_min")
  )
}

.plot_gate_uv_marker_add_adj <- function(exc_min,
                                         plot_tbl,
                                         dens_obj_raw,
                                         prob_g_min) {
  if (!exc_min) {
    return(plot_tbl)
  }
  plot_tbl[, "type"] <- "raw"
  dens_obj_adj <- dens_obj_raw
  dens_obj_adj$y <- dens_obj_adj$y * prob_g_min
  plot_tbl_adj <- tibble::tibble(
    x = dens_obj_adj$x, y = dens_obj_adj$y, type = "adj"
  )
  plot_tbl |>
    dplyr::bind_rows(plot_tbl_adj)
}

.plot_gate_uv_marker_plot <- function(plot_tbl,
                                      exc_min,
                                      ind,
                                      ind_lab,
                                      marker,
                                      marker_lab,
                                      show_gate,
                                      path_project) {
  p <- .plot_gate_uv_marker_plot_init(plot_tbl, exc_min, ind, ind_lab)
  p <- .plot_add_axis_title(p, marker, marker_lab)
  p <- p + ggplot2::labs(y = "Density")
  p <- .plot_add_title(p, marker, NULL, marker_lab)
  p <- .plot_add_gate(p, .data, ind, marker, path_project, show_gate)
  p
}

.plot_gate_uv_marker_plot_init <- function(plot_tbl,
                                           exc_min,
                                           ind,
                                           ind_lab) {
  alpha_lab_vec <- c("raw" = 0.5, "adj" = 1)
  p <- if (exc_min) {
    if (length(ind) > 1L) {
      ggplot(plot_tbl, aes(x = x, y = y, alpha = type, color = ind_lab)) +
        scale_alpha_manual(values = alpha_lab_vec)
    } else {
      ggplot(plot_tbl, aes(x = x, y = y, alpha = type)) +
        scale_alpha_manual(values = alpha_lab_vec)
    }
  } else {
    if (length(ind) > 1L) {
      ggplot(plot_tbl, aes(x = x, y = y, colour = ind_lab)) +
        scale_alpha_manual(values = alpha_lab_vec)
    } else {
      ggplot(plot_tbl, aes(x = x, y = y)) +
        scale_alpha_manual(values = alpha_lab_vec)
    }
  }
  p +
    geom_line() +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "x") +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white")
    )
}

.plot_grid <- function(plot,
                       p_list,
                       n_col) {
  if (!plot) {
    return(p_list)
  }
  cowplot::plot_grid(
    plotlist = p_list, ncol = n_col,
    align = "hv"
  ) +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white")
    )
}
