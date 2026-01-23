#' Plot stimulation gate
#'
#' Plot bivariate hex and univariate density plots for batches of samples, along
#' with their gates.
#'
#' @param ind numeric vector. Specifies indices in `.data` to plot.
#' @param .data GatingSet. Same GatingSet passed to `stimgate_gate`.
#' @param path_project character.
#' Path to the project directory used for `stimgate_gate`.
#' @param marker character vector of length one or two. Specifies markers
#' to be plotted. If only one is passed, then only univariate plots are created.
#' @param chnl character vector of length one or two. Specifies channels
#' to be plotted. Ignored if `marker` is provided.
#' @param pop character. Specifies population within GatingSet that
#' gates were calculated on. If `NULL`, defaults to population specified
#' by folder name in `project_path/gates/pop_<pop>`, but throws
#' an error if more than one population is detected (i.e. more
#' than one directory in `gates/`). Default is `NULL`.
#' @param ind_lab named character vector.
#' Labels for `ind` used in plot.
#' Optional.
#' @param axis_lab named character vector.
#' Labels for axis titles, applied to `marker` or `chnl`.
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
#' @inheritParams stimgate_data_get_ex
#'
#' @return A grid of plots if `grid` is TRUE, otherwise a list of ggplot objects.
#'
#' @import ggplot2
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
                          marker = NULL,
                          chnl = NULL,
                          pop = NULL,
                          ind_lab = NULL,
                          axis_lab = NULL,
                          exc_min = TRUE,
                          limits_expand = NULL,
                          limits_equal = FALSE,
                          grid = TRUE,
                          grid_n_col = 2,
                          show_gate = TRUE,
                          min_cell = 10,
                          bias = FALSE,
                          combn_exc = NULL,
                          chnl_gate = NULL,
                          marker_gate = NULL,
                          gate_type_cyt_pos = "cyt",
                          gate_type_single_pos = "single",
                          mult = FALSE,
                          gate_uns_method = "min") {
  if (is.null(marker) && is.null(chnl)) {
    stop("Must specify one of marker or chnl")
  }
  pop <- pop %||% .gate_get_pop(path_project)
  if (length(pop) > 1L) {
    stop("Cannot plot gates for multiple populations")
  }
  if (length(pop) == 0L || !nzchar(pop)) {
    stop("No population found for plotting gates")
  }
  p_list <- .plot_gate(
    ind = ind, ind_lab = ind_lab, .data = .data,
    marker = marker, chnl = chnl, pop = pop, axis_lab = axis_lab,
    path_project = path_project, exc_min = exc_min,
    limits_expand = limits_expand, limits_equal = limits_equal,
    show_gate = show_gate, min_cell = min_cell,
    bias = bias, combn_exc = combn_exc, chnl_gate = chnl_gate,
    marker_gate = marker_gate, gate_type_cyt_pos = gate_type_cyt_pos,
    gate_type_single_pos = gate_type_single_pos, mult = mult,
    gate_uns_method = gate_uns_method
  )
  if (length(p_list) == 0L) {
    return(NULL)
  }
  .plot_grid(plot = grid, p_list = p_list, n_col = grid_n_col)
}

#' @keywords internal
.plot_gate <- function(marker,
                       chnl,
                       pop,
                       ind,
                       ind_lab,
                       .data,
                       axis_lab,
                       path_project,
                       exc_min,
                       limits_expand,
                       limits_equal,
                       show_gate,
                       min_cell,
                       bias,
                       combn_exc,
                       chnl_gate,
                       marker_gate,
                       gate_type_cyt_pos,
                       gate_type_single_pos,
                       mult,
                       gate_uns_method) {
  # bv
  p_list_bv <- .plot_gate_bv(
    marker = marker, chnl = chnl, pop = pop,
    ind = ind, ind_lab = ind_lab,
    .data = .data, axis_lab = axis_lab,
    path_project = path_project, exc_min = exc_min,
    limits_expand = limits_expand, limits_equal = limits_equal,
    show_gate = show_gate, min_cell = min_cell,
    bias = bias, combn_exc = combn_exc, chnl_gate = chnl_gate,
    marker_gate = marker_gate, gate_type_cyt_pos = gate_type_cyt_pos,
    gate_type_single_pos = gate_type_single_pos, mult = mult,
    gate_uns_method = gate_uns_method
  )

  # uv
  p_list_uv <- .plot_gate_uv(
    ind = ind, ind_lab = ind_lab, .data = .data,
    marker = marker, chnl = chnl, pop = pop,
    exc_min = exc_min, axis_lab = axis_lab,
    show_gate = show_gate, path_project = path_project,
    min_cell = min_cell,
    bias = bias, combn_exc = combn_exc, chnl_gate = chnl_gate,
    marker_gate = marker_gate, gate_type_cyt_pos = gate_type_cyt_pos,
    gate_type_single_pos = gate_type_single_pos, mult = mult,
    gate_uns_method = gate_uns_method
  )

  p_list_bv |> append(p_list_uv)
}

#' @keywords internal
.plot_gate_bv <- function(marker,
                          chnl,
                          pop,
                          ind,
                          ind_lab,
                          .data,
                          axis_lab,
                          path_project,
                          exc_min,
                          limits_expand,
                          limits_equal,
                          show_gate,
                          min_cell,
                          bias,
                          combn_exc,
                          chnl_gate,
                          marker_gate,
                          gate_type_cyt_pos,
                          gate_type_single_pos,
                          mult,
                          gate_uns_method) {
  one_chnl <- is.null(marker) && !is.null(chnl) && length(chnl) == 1L
  one_marker <- is.null(chnl) && !is.null(marker) && length(marker) == 1L
  one_var <- one_chnl || one_marker
  if (one_var) {
    return(NULL)
  }
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    if (interactive()) {
      prompt_answer <- readline(
        prompt = paste0(
          "The 'hexbin' package is required for bivariate plots. ", # nolint linter_line_length_linter
          "Do you want to install it now? [y/n]: "
        )
      )
      if (tolower(prompt_answer) != "y") {
        stop("Cannot proceed without installing 'hexbin' package.")
      }
      utils::install.packages("hexbin")
    } else {
      stop("The 'hexbin' package is required but not installed.")
    }
  }
  p_list <- lapply(seq_along(ind), function(i) {
    ind_curr <- ind[[i]]
    ex_tbl <- .plot_get_ex_tbl(
      ind_curr,
      .data, pop, marker, chnl,
      exc_min, path_project, bias,
      combn_exc, chnl_gate, marker_gate,
      gate_type_cyt_pos, gate_type_single_pos,
      mult, gate_uns_method
    )
    if (nrow(ex_tbl) < min_cell) {
      return(NULL)
    }
    p <- plot_cyto(
      data = ex_tbl,
      marker = chnl %||% marker,
      exc_min = FALSE,
      limits_expand = limits_expand,
      limits_equal = limits_equal
    ) +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")
      )
    p <- .plot_add_axis_title(p, marker, chnl, axis_lab)
    p <- .plot_add_title(p, ind_curr, i, ind_lab)
    p <- .plot_add_gate(
      p, ind_curr, marker, chnl, pop, path_project, show_gate
    )
    p
  }) |>
    stats::setNames(.plot_get_lab(ind, ind_lab))
  p_list <- p_list[vapply(p_list, Negate(is.null), logical(1))]
  if (length(p_list) == 0L) {
    return(NULL)
  }
  p_list
}

#' @keywords internal
.plot_get_ex_tbl <- function(ind,
                             .data,
                             pop,
                             marker,
                             chnl,
                             exc_min,
                             path_project,
                             bias,
                             combn_exc,
                             chnl_gate,
                             marker_gate,
                             gate_type_cyt_pos,
                             gate_type_single_pos,
                             mult,
                             gate_uns_method) {
  lapply(ind, function(ind_curr) {
    stimgate_data_get_ex(
      path_project, .data, pop, ind_curr, chnl,
      marker, bias, exc_min, combn_exc, chnl_gate,
      marker_gate, gate_type_cyt_pos,
      gate_type_single_pos, mult, gate_uns_method
    )
  }) |>
    Reduce(rbind, x = _)
}


#' @keywords internal
.plot_get_ex_tbl_ind <- function(ind,
                                 .data,
                                 pop,
                                 chnl,
                                 chnl_lab,
                                 exc_min,
                                 path_project) {
  ex <- if (!is.null(.data)) {
    .get_ex_new(.data, pop, chnl, ind, path_project, FALSE)
  } else {
    .get_ex_old(pop, chnl, ind, path_project)
  }
  ex <- .plot_get_ex_tbl_ind_exc_min(ex, exc_min, chnl)
  if (!is.null(chnl_lab)) {
    colnames(ex) <- chnl_lab[colnames(ex)]
  }
  ex
}

#' @keywords internal
.plot_get_ex_tbl_ind_new <- function(ind,
                                     .data,
                                     pop,
                                     chnl) {
  fr <- flowWorkspace::gh_pop_get_data(.data[[ind]], y = pop)
  ex_tbl <- flowCore::exprs(fr) |> tibble::as_tibble()
  ex_tbl[, chnl, drop = FALSE]
}

#' @keywords internal
.plot_get_ex_tbl_ind_old <- function(ind,
                                     pop,
                                     chnl) {}

#' @keywords internal
.plot_get_ex_tbl_ind_exc_min <- function(ex_tbl, exc_min, chnl) {
  if (!exc_min) {
    return(ex_tbl)
  }
  n_row_init <- nrow(ex_tbl)
  attr(ex_tbl, "n_row_init") <- n_row_init
  min_val_vec <- vapply(
    chnl, function(x) min(ex_tbl[[x]], na.rm = TRUE), numeric(1)
  )
  for (i in seq_along(chnl)) {
    ex_tbl <- ex_tbl[ex_tbl[[chnl[i]]] > min_val_vec[i], ]
  }
  n_row_final <- nrow(ex_tbl)
  attr(ex_tbl, "prob_g_min") <- n_row_final / n_row_init

  ex_tbl
}

#' @keywords internal
.plot_add_axis_title <- function(p, val1, val2, val_lab) {
  val <- if (!is.null(val1)) {
    val1
  } else {
    val2
  }
  lab <- .plot_get_lab(val, val_lab)
  p <- p + labs(x = lab[[1]])
  if (length(lab) > 1L) {
    p <- p + labs(y = lab[[2]])
  }
  p
}

#' @keywords internal
.plot_add_title <- function(p, ind, i, ind_lab) {
  p + ggtitle(.plot_get_lab(ind, ind_lab, i))
}

#' @keywords internal
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

#' @keywords internal
.plot_add_gate <- function(p,
                           ind,
                           marker,
                           chnl,
                           pop,
                           path_project,
                           show_gate) {
  if (!show_gate) {
    return(p)
  }
  marker <- marker %||% stimgate_meta_read_chnl_lab(path_project)[chnl]
  chnl <- chnl %||% stimgate_meta_read_marker_lab(path_project)[marker]
  pop <- pop %||% .gate_get_pop(path_project)
  if (length(pop) > 1L) {
    stop("Cannot plot gates for multiple populations")
  }
  if (length(pop) == 0L || !nzchar(pop)) {
    stop("No population found for plotting gates")
  }
  gate_tbl <- .plot_get_gate_tbl(ind, pop, marker, chnl, path_project)
  chnl_gate <- chnl[chnl %in% .gate_get_chnl(path_project, pop)]
  for (i in seq_along(chnl_gate)) {
    gate_vec <- gate_tbl[["gate"]][gate_tbl[["chnl"]] == chnl_gate[i]]
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

#' @keywords internal
.plot_get_gate_tbl <- function(ind, pop, marker, chnl, path_project) {
  gate_tbl <- stimgate_gate_get(
    path_project,
    pop,
    chnl = chnl,
    marker = marker
  ) |>
    dplyr::group_by(gate_name, chnl, marker, ind, batch) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  gate_tbl <- gate_tbl[gate_tbl[["ind"]] %in% ind, ]
  ind_vec <- NULL
  if (!is.null(marker)) {
    for (i in seq_along(marker)) {
      ind_vec[[i]] <- which(gate_tbl[["marker"]] == marker[i])
    }
  } else {
    for (i in seq_along(chnl)) {
      ind_vec[[i]] <- which(gate_tbl[["chnl"]] == chnl[i])
    }
  }
  ind_vec <- ind_vec |> unlist()
  gate_tbl[ind_vec, ]
}

#' @keywords internal
.plot_gate_uv <- function(ind,
                          ind_lab,
                          .data,
                          marker,
                          chnl,
                          pop,
                          exc_min,
                          axis_lab,
                          show_gate,
                          path_project,
                          min_cell,
                          bias,
                          combn_exc,
                          chnl_gate,
                          marker_gate,
                          gate_type_cyt_pos,
                          gate_type_single_pos,
                          mult,
                          gate_uns_method) {
  var_loop <- if (!is.null(marker)) marker else chnl
  p_list <- lapply(var_loop, function(v) {
    marker_curr <- if (!is.null(marker)) v else NULL
    chnl_curr <- if (!is.null(chnl)) v else NULL
    .plot_gate_uv_marker(
      marker = marker_curr, chnl = chnl_curr,
      ind = ind, .data = .data, pop = pop,
      exc_min = exc_min, ind_lab = ind_lab,
      axis_lab = axis_lab, show_gate = show_gate,
      path_project = path_project, min_cell = min_cell,
      bias = bias, combn_exc = combn_exc, chnl_gate = chnl_gate,
      marker_gate = marker_gate, gate_type_cyt_pos = gate_type_cyt_pos,
      gate_type_single_pos = gate_type_single_pos, mult = mult,
      gate_uns_method = gate_uns_method
    )
  }) |>
    stats::setNames(.plot_get_lab(var_loop, axis_lab))
  p_list <- p_list[vapply(p_list, Negate(is.null), logical(1))]
  if (length(p_list) == 0L) {
    return(NULL)
  }
  p_list
}

#' @keywords internal
.plot_gate_uv_marker <- function(marker,
                                 chnl,
                                 pop,
                                 ind,
                                 .data,
                                 exc_min,
                                 ind_lab,
                                 axis_lab,
                                 show_gate,
                                 path_project,
                                 min_cell,
                                 bias,
                                 combn_exc,
                                 chnl_gate,
                                 marker_gate,
                                 gate_type_cyt_pos,
                                 gate_type_single_pos,
                                 mult,
                                 gate_uns_method) {
  plot_tbl <- .plot_gate_uv_marker_get_plot_tbl(
    marker = marker, chnl = chnl, pop = pop,
    ind = ind, .data = .data, exc_min = exc_min,
    ind_lab = ind_lab, min_cell = min_cell,
    path_project = path_project,
    bias = bias, combn_exc = combn_exc, chnl_gate = chnl_gate,
    marker_gate = marker_gate, gate_type_cyt_pos = gate_type_cyt_pos,
    gate_type_single_pos = gate_type_single_pos, mult = mult,
    gate_uns_method = gate_uns_method
  )
  if (is.null(plot_tbl)) {
    return(NULL)
  }
  .plot_gate_uv_marker_plot(
    plot_tbl = plot_tbl, exc_min = exc_min,
    ind = ind, ind_lab = ind_lab, pop = pop,
    marker = marker, chnl = chnl, axis_lab = axis_lab,
    show_gate = show_gate, path_project = path_project
  )
}

#' @keywords internal
.plot_gate_uv_marker_get_plot_tbl <- function(ind,
                                              .data,
                                              marker,
                                              chnl,
                                              pop,
                                              exc_min,
                                              ind_lab,
                                              min_cell,
                                              path_project,
                                              bias,
                                              combn_exc,
                                              chnl_gate,
                                              marker_gate,
                                              gate_type_cyt_pos,
                                              gate_type_single_pos,
                                              mult,
                                              gate_uns_method) {
  plot_tbl_list <- lapply(seq_along(ind), function(i) {
    plot_tbl <- .plot_gate_uv_marker_get_plot_tbl_ind(
      ind = ind[[i]], .data = .data, pop = pop,
      marker = marker, chnl = chnl, exc_min = exc_min,
      min_cell = min_cell, path_project = path_project,
      bias = bias, combn_exc = combn_exc, chnl_gate = chnl_gate,
      marker_gate = marker_gate, gate_type_cyt_pos = gate_type_cyt_pos,
      gate_type_single_pos = gate_type_single_pos, mult = mult,
      gate_uns_method = gate_uns_method
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

#' @keywords internal
.plot_gate_uv_marker_get_plot_tbl_ind <- function(ind,
                                                  .data,
                                                  marker,
                                                  chnl,
                                                  pop,
                                                  exc_min,
                                                  min_cell,
                                                  path_project,
                                                  bias,
                                                  combn_exc,
                                                  chnl_gate,
                                                  marker_gate,
                                                  gate_type_cyt_pos,
                                                  gate_type_single_pos,
                                                  mult,
                                                  gate_uns_method) {
  ex_tbl <- stimgate_data_get_ex(
    path_project, .data, pop, ind, chnl,
    marker, bias, exc_min, combn_exc, chnl_gate,
    marker_gate, gate_type_cyt_pos,
    gate_type_single_pos, mult, gate_uns_method
  )
  if (nrow(ex_tbl) < min_cell) {
    return(NULL)
  }
  .var <- if (!is.null(marker)) marker else chnl
  dens_obj_raw <- density(ex_tbl[[.var]], na.rm = TRUE)
  plot_tbl <- tibble::tibble(x = dens_obj_raw$x, y = dens_obj_raw$y)
  .plot_gate_uv_marker_add_adj(
    exc_min, plot_tbl, dens_obj_raw, ex_tbl
  )
}

#' @keywords internal
.plot_gate_uv_marker_add_adj <- function(exc_min,
                                         plot_tbl,
                                         dens_obj_raw,
                                         ex_tbl) {
  if (!exc_min) {
    return(NULL)
  }
  prob_g_min <- attr(ex_tbl, "prob_g_min")[[1]][[1]][[1]]
  plot_tbl[, "type"] <- "raw"
  dens_obj_adj <- dens_obj_raw
  dens_obj_adj$y <- dens_obj_adj$y * prob_g_min
  plot_tbl_adj <- tibble::tibble(
    x = dens_obj_adj$x, y = dens_obj_adj$y, type = "adj"
  )
  plot_tbl |>
    dplyr::bind_rows(plot_tbl_adj)
}

#' @keywords internal
.plot_gate_uv_marker_plot <- function(plot_tbl,
                                      exc_min,
                                      ind,
                                      ind_lab,
                                      marker,
                                      chnl,
                                      pop,
                                      axis_lab,
                                      show_gate,
                                      path_project) {
  p <- .plot_gate_uv_marker_plot_init(plot_tbl, exc_min, ind, ind_lab)
  p <- .plot_add_axis_title(p, marker, chnl, axis_lab)
  p <- p + ggplot2::labs(y = "Density")
  .var <- if (!is.null(marker)) marker else chnl
  p <- .plot_add_title(p, .var, NULL, axis_lab)
  p <- .plot_add_gate(p, ind, marker, chnl, pop, path_project, show_gate)
  p
}

#' @keywords internal
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

#' @keywords internal
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
