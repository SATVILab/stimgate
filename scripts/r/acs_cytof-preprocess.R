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

  sampleMap <- .acsCytofMapFiles(fcs_vec)
  .acsCytofBatchList(sampleMap)
  preprocessing <- list(
    settings = list(transform = "asinh(x / 5)", pairingVersion = 1L),
    sampleMap = sampleMap,
    inputFileListHash = .acsCytofHash(basename(fcs_vec)),
    inputContentHash = .acsCytofHash(unname(tools::md5sum(fcs_vec)))
  )
  cs <- flowWorkspace::load_cytoset_from_fcs(fcs_vec)
  flowWorkspace::sampleNames(cs) <- basename(fcs_vec)

  gs <- flowWorkspace::GatingSet(cs)
  forwardTransform <- function(x) {
    asinh(x / 5)
  }
  backTransform <- function(x) {
    5 * sinh(x)
  }
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
  # Remove the old manifest before the swap, so an interrupted run leaves it
  # missing rather than paired with the new GatingSet.
  path_manifest <- .acsCytofPreprocessingFile(path_gs)
  path_manifest_tmp <- paste0(path_manifest, ".tmp-", Sys.getpid())
  on.exit(unlink(path_manifest_tmp), add = TRUE)
  saveRDS(preprocessing, path_manifest_tmp)
  unlink(path_manifest)
  .acsCytofReplaceDir(path_gs, function(path_tmp) {
    flowWorkspace::save_gs(gs = gs_trans, path = path_tmp)
  })
  if (!file.rename(path_manifest_tmp, path_manifest)) {
    stop("Could not move the preprocessing manifest into place: ", path_manifest)
  }
  path_gs
}

.acsCytofPlotGatingSetCheck <- function(path_gs) {
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
  plots <- lapply(unique(expr_tbl_long$trans), function(x) {
    plot_tbl <- expr_tbl_long |> dplyr::filter(trans == x)
    ggplot2::ggplot(plot_tbl, ggplot2::aes(x = expr, fill = marker)) +
      .analysis_theme(grid = "x") +
      ggplot2::geom_histogram(bins = 30) +
      ggplot2::facet_wrap(~marker, scales = "free", ncol = 8) +
      ggplot2::scale_x_continuous(
        labels = .analysis_label_number,
        n.breaks = 3,
        guide = ggplot2::guide_axis(check.overlap = TRUE)
      ) +
      ggplot2::theme(
        legend.position = "none",
        strip.text = ggplot2::element_text(size = 7),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      ) +
      ggplot2::labs(y = "Count", x = "Marker expression")
  })
  stats::setNames(plots, unique(expr_tbl_long$trans))
}

plot_gatingset_check <- function(path_gs, path_plot_dir, plots = NULL) {
  if (is.null(plots)) plots <- .acsCytofPlotGatingSetCheck(path_gs)
  dir.create(path_plot_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- vapply(names(plots), function(trans) {
    path_plot <- file.path(
      path_plot_dir, paste0("all_markers-trans_", trans, ".png")
    )
    .analysis_save_fig(plots[[trans]], path_plot, height = 12)
    path_plot
  }, character(1))
  unname(paths)
}
