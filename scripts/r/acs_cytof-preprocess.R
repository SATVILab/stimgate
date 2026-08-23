create_gatingset <- function(
  path_fcs,
  path_gs,
  n_sample = NULL,
  ind_sample = NULL
) {
  if (!is.null(n_sample) && !is.null(ind_sample)) {
    stop("Specify only one of n_sample and ind_sample.")
  }

  fcs_vec <- list.files(
    path_fcs,
    pattern = "\\.fcs$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  ) |>
    sort()
  if (length(fcs_vec) == 0L) {
    stop("No FCS files found at: ", path_fcs)
  }

  if (!is.null(ind_sample)) {
    if (anyNA(ind_sample) ||
      any(ind_sample < 1L) ||
      any(ind_sample > length(fcs_vec))) {
      stop("ind_sample contains an index outside the available FCS files.")
    }
    fcs_vec <- fcs_vec[ind_sample]
  } else if (!is.null(n_sample)) {
    if (length(n_sample) != 1L ||
      is.na(n_sample) ||
      n_sample < 1L ||
      n_sample > length(fcs_vec)) {
      stop("n_sample must be between 1 and the number of available FCS files.")
    }
    fcs_vec <- fcs_vec[seq_len(n_sample)]
  }

  cs <- flowWorkspace::load_cytoset_from_fcs(fcs_vec)

  gs <- flowWorkspace::GatingSet(cs)
  forwardTransform <- function(x) {
    asinh(x / 5)
  }
  backTransform <- function(x) {
    5 * sinh(x)
  }
  round(2.3, 10) == round(backTransform(forwardTransform(2.3)), 10)
  round(2.3, 10) == round(forwardTransform(backTransform(2.3)), 10)
  trans.obj <- flowWorkspace::flow_trans(
    "asinh",
    trans.fun = forwardTransform,
    inverse.fun = backTransform
  )
  trans <- flowWorkspace::transformerList(
    c(
      "Dy161Di",
      "Dy162Di",
      "Dy163Di",
      "Dy164Di",
      "Er166Di",
      "Er167Di",
      "Er168Di",
      "Er170Di",
      "Eu151Di",
      "Eu153Di",
      "Gd155Di",
      "Gd156Di",
      "Gd158Di",
      "Gd160Di",
      "Ho165Di",
      "Lu175Di",
      "Lu176Di",
      "Nd142Di",
      "Nd143Di",
      "Nd144Di",
      "Nd145Di",
      "Nd146Di",
      "Nd148Di",
      "Nd150Di",
      "Pr141Di",
      "Sm147Di",
      "Sm149Di",
      "Sm152Di",
      "Sm154Di",
      "Tb159Di",
      "Tm169Di",
      "Yb171Di",
      "Yb172Di",
      "Yb173Di",
      "Yb174Di"
    ),
    trans.obj
  )
  gs_trans <- flowWorkspace::transform(gs, trans)
  if (dir.exists(path_gs)) {
    unlink(path_gs, recursive = TRUE)
  }
  if (!dir.exists(dirname(path_gs))) {
    dir.create(dirname(path_gs), recursive = TRUE)
  }
  flowWorkspace::save_gs(
    gs = gs_trans,
    path = path_gs
  )
  path_gs
}

plot_gatingset_check <- function(path_gs, path_plot_dir) {
  gs <- flowWorkspace::load_gs(path_gs)
  cf <- flowWorkspace::gh_pop_get_data(gs[[1]])
  chnl_to_marker <- UtilsCytoRSV::chnl_to_marker(cf)
  expr_tbl <- flowCore::exprs(cf) |>
    tibble::as_tibble()
  forwardTransform <- function(x) asinh(x / 5)
  backTransform <- function(x) 5 * sinh(x)
  trans_obj_asinh <- flowWorkspace::flow_trans(
    "asinh",
    trans.fun = forwardTransform,
    inverse.fun = backTransform
  )
  expr_tbl_long <- expr_tbl |>
    dplyr::mutate(cell_id = seq_len(dplyr::n())) |>
    tidyr::pivot_longer(
      -cell_id,
      names_to = "chnl",
      values_to = "expr"
    ) |>
    dplyr::mutate(marker = chnl_to_marker[chnl] |> as.character()) |>
    dplyr::mutate(trans = "asinh")
  expr_tbl_long <- expr_tbl_long |>
    dplyr::bind_rows(
      expr_tbl_long |>
        dplyr::mutate(trans = "none") |>
        dplyr::mutate(expr = backTransform(expr))
    )
  try(
    lapply(unique(expr_tbl_long$trans), function(x) {
      plot_tbl <- expr_tbl_long |> dplyr::filter(trans == x)
      p <- ggplot(
        plot_tbl,
        aes(x = expr, fill = marker)
      ) +
        cowplot::theme_cowplot() +
        cowplot::background_grid(major = "x") +
        theme(
          plot.background = element_rect(fill = "white", colour = "white"),
          panel.background = element_rect(fill = "white", colour = "white")
        ) +
        geom_histogram(bins = 30) +
        facet_wrap(~marker, scales = "free", ncol = 8) +
        theme(legend.position = "none") +
        theme(
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()
        ) +
        labs(y = "Count", x = "Marker expression", title = x)

      path_plot <- file.path(
        path_plot_dir,
        paste0("all_markers-trans_", x, ".png")
      )

      if (file.exists(path_plot)) {
        invisible(file.remove(path_plot))
      }
      if (!dir.exists(dirname(path_plot))) {
        dir.create(dirname(path_plot), recursive = TRUE)
      }

      ggplot2::ggsave(
        filename = path_plot,
        plot = p,
        width = 20,
        height = 16,
        units = "cm"
      )

      path_plot
    }) |>
      unlist()
  )
}
