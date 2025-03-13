plot_gate <- function(ind,
                      ind_lab = NULL,
                      .data,
                      path_project,
                      marker,
                      marker_lab = NULL,
                      exc_min = TRUE,
                      limits_expand = NULL,
                      limits_equal = NULL,
                      grid = TRUE,
                      grid_n_col = 2) {
  p_list <- .plot_gate(
    ind = ind, ind_lab = ind_lab, .data = .data,
    marker = marker, marker_lab = marker_lab,
    path_project = path_project, exc_min = exc_min,
    limits_expand = limits_expand, limits_equal = limits_equal
  )
  if (length(p_list) == 0L) {
    return(NULL)
  }
  .plot_grid(plot_grid = grid, p_list = p_list, n_col = grid_n_col)
}

.plot_gate <- function(marker,
                       ind,
                       ind_lab,
                       .data,
                       marker_lab,
                       path_project,
                       exc_min,
                       limits_expand,
                       limits_equal) {
  # bv
  p_list_bv <- .plot_gate_bv(
    marker = marker, ind = ind, ind_lab = ind_lab,
    .data = .data, marker_lab = marker_lab,
    path_project = path_project, exc_min = exc_min,
    limits_expand = limits_expand, limits_equal = limits_equal
  )

  # uv
  p_list_uv <- .plot_gate_uv(
    ind = ind, ind_lab = ind_lab, .data = .data,
    marker = marker, exc_min = exc_min
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
                          limits_equal) {
  if (length(marker) == 1L) {
    return(NULL)
  }
  p_list <- lapply(seq_along(ind), function(i) {
    ind_curr <- ind[[i]]
    ex_tbl <- .plot_get_ex_tbl(ind_curr, .data, marker, exc_min = FALSE)
    if (nrow(ex_tbl) == 0L) {
      return(NULL)
    }
    p <- UtilsCytoRSV::plot_cyto(
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
    p <- .plot_bv_add_axis_title(p, marker, marker_lab)
    p <- .plot_bv_add_title(p, ind_curr, i, ind_lab)
    p
  })
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
    Reduce(rbind, x = _) |>
}


.plot_get_ex_tbl_ind <- function(ind, .data, marker, exc_min) {
  fr <- flowWorkspace::gh_pop_get_data(.data[[ind]])
  ex_tbl <- flowCore::exprs(fr) |> tibble::as_tibble()
  ex_tbl <- ex_tbl[, marker_vec, drop = FALSE]
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

.plot_bv_add_axis_title <- function(p, marker, marker_lab) {
  marker_lab <- .plot_get_marker_lab(marker, marker_lab)
  p <- p + labs(x = marker_lab[[1]])
  if (length(marker_lab) > 1L) {
    p <- p + labs(y = marker_lab[[2]])
  }
  p
}

.plot_get_marker_lab <- function(marker, marker_lab) {
  if (is.null(marker_lab)) {
    return(marker)
  }
  if (!is.null(names(marker_lab))) {
    marker_lab <- marker_lab[marker]
  }
  marker_lab
}

.plot_bv_add_title <- function(p, ind, i, ind_lab) {
  p + ggtitle(.plot_get_ind_lab(ind, i, ind_lab))
}

.plot_get_ind_lab <- function(ind, ind_lab, i) {
  if (is.null(ind_lab)) {
    return(ind)
  }
  if (!is.null(names(ind_lab))) ind_lab[[ind_curr]] else ind_lab[[i]]
}

.plot_add_gate <- function(p, ind, marker, path_project) {
  gate_tbl <- .plot_get_gate_tbl(ind, marker, path_project)

}

.plot_get_gate_tbl <- function(ind, marker, path_project) {
  gate_tbl <- get_gates(path_project)
  gate_tbl[gate_tbl[["ind"]] == ind & gate_tbl[["chnl"]] %in% marker,]
}

.plot_gate_uv <- function(ind,
                          ind_lab,
                          .data,
                          marker,
                          exc_min) {
  P_list <- lapply(marker, function(i) {
    .plot_gate_uv_marker(
      marker = m, ind = ind, .data = .data,
      exc_min = exc_min, ind_lab = ind_lab, i = i
    )
  })
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
                                 i) {
  plot_tbl <- .plot_gate_uv_marker_get_plot_tbl(
    ind = ind, .data = .data, marker = marker, exc_min = exc_min,
    ind_lab = ind_lab
  )
  .plot_gate_uv_marker_plot(
    plot_tbl = plot_tbl, exc_min = exc_min,
    ind = ind, i = i, ind_lab = ind_lab
  )

}

.plot_gate_uv_marker_get_plot_tbl <- function(ind,
                                              .data,
                                              marker,
                                              exc_min,
                                              ind_lab) {
  lapply(ind, function(x) {
    .plot_gate_uv_marker_get_plot_tbl_ind(
      ind = x, .data = .data, marker = marker, exc_min = exc_min
    )
    ind_lab <- .plot_get_ind_lab(x, ind_lab)
    plot_tbl[, "ind"] <- ind_lab
    plot_tbl
  }) |>
    Reduce(rbind, x = _)
}

.plot_gate_uv_marker_get_plot_tbl_ind <- function() {
  ex_tbl <- .plot_get_ex_tbl(ind, .data, marker, exc_min)
  dens_obj_raw <- density(ex_tbl[[marker]], na.rm = TRUE)
  plot_tbl <- tibble::tibble(x = dens_obj_raw$x, y = dens_obj_raw$y)
  .plot_gate_uv_marker_add_adj(exc_min, plot_tbl, dens_obj_raw)
}

.plot_gate_uv_marker_add_adj <- function(exc_min, plot_tbl, dens_obj_raw) {
  if (!exc_min) {
    return(plot_tbl)
  }
  plot_tbl <- plot_tbl |> dplyr::mutate(type = raw)
  dens_obj_adj <- dens_obj_raw
  dens_obj_adj$y <- dens_obj_adj$y * attr(ex_tbl, "prob_g_min")
  plot_tbl_adj <- tibble::tibble(
    x = dens_obj_adj$x, y = dens_obj_adj$y, type = "adj"
  )
  plot_tbl |>
    dplyr::bind_rows(plot_tbl_adj)
}

.plot_gate_uv_marker_plot <- function(plot_tbl, exc_min, ind, i, ind_lab) {
  p <- .plot_gate_uv_marker_plot_init(plot_tbl, exc_min)
  p <- .plot_add_axis_title(p, marker, marker_lab)
  p <- .plot_add_title(p, ind, i, ind_lab)
  p
}

.plot_gate_uv_marker_plot_init <- function(plot_tbl, exc_min) {
  p <- if (exc_min) {
    ggplot(plot_tbl, aes(x = x, y = y, color = type))
  } else {
    ggplot(plot_tbl, aes(x = x, y = y))
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

.plot_grid <- function(plot_grid,
                       p_list,
                       n_col) {
  if (!plot_grid) {
    return(p_list)
  }
  cowplot::plot_grid(
    plotlist = plot_grid, ncol = n_col
  ) +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white")
    )
}