library(testthat)

test_that("stimgate_gate runs", {
  path_gs <- testthat::test_path(
    "testdata", "postmortem", "gs", "gs_cytof_acs_cd4"
  )
  gs <- flowWorkspace::load_gs(path_gs)
  fn_tbl_info <- get_fn_tbl_info_postmortem(gs)
  # .debugonce(get_batch_list_postmortem)
  # browser()
  min_cell <- TRUE
  filter_method <- if (min_cell) "min_cell" else "min_uns"
  batch_list <- get_batch_list_postmortem(
    fn_tbl_info,
    col_grp = c("pid", "location"),
    col_stim = "stim",
    uns_chr = "uns",
    min_cell_uns = 100,
    min_cell_stim = 100,
    col_n_cell = "n_cell_pop",
    filter_method = filter_method
  )

  marker_vec <- c("Er168Di", "Lu175Di")
  # marker_vec <- c("Er168Di", "Lu175Di", "Er166Di", "Yb172Di", "Nd150Di")
  path_project <- setup_project_postmortem(min_cell)
  # browser()
  # .debugonce(stimgate_gate)
  # .debugonce(.gate_stats)

  stimgate_gate(
    path_project = path_project,
    .data = gs,
    batch_list = batch_list,
    marker = marker_vec,
    .debug = TRUE,
    tol_clust = 1e-6,
    min_cell = 40
  )
  browser()
  browser()
  library(ggplot2)

  chnl_lab <- .get_chnl_lab(gs)
  ind_lab <- .get_ind_lab(fn_tbl_info)
  # .debugonce(.plot_get_gate_tbl)
  # START HERE:
  # Handle multiple gates for multiple indices,
  # all to be plotted, for the uv plots
  # (and add an alpha parameter, so that they can
  # all be seenk)
  # .debug(.plot_add_gate)

  p_grid <- plot_gate(
    ind = batch_list[[1]],
    ind_lab = ind_lab,
    .data = gs,
    path_project = path_project,
    marker = marker_vec,
    limits_expand = list(c(0, 4)),
    limits_equal = TRUE,
    marker_lab = chnl_lab
  )
  path_dir_fig <- here::here("_tmp", "fig", "p_gate")
  ggplot2::ggsave(
    file.path(path_dir_fig, "p_grid.png"),
    p_grid,
    units = "cm",
    width = 20,
    height = 20,
  )
})
