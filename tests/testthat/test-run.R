library(testthat)

test_that("stimgate_gate runs", {
  path_gs <- testthat::test_path(
    "testdata", "postmortem", "gs", "gs_cytof_acs_cd4"
  )
  gs <- flowWorkspace::load_gs(path_gs)
  fn_tbl_info <- get_fn_tbl_info_postmortem(gs)
  # debugonce(get_batch_list_postmortem)
  UtilsCytoRSV::chnl_lab(gs[[1]])
  browser()
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
  marker_vec <- c("Er168Di", "Lu175Di", "Er166Di", "Yb172Di", "Nd150Di")
  path_project <- setup_project_postmortem(min_cell)
  browser()
  browser()
  browser()
  debugonce(stimgate_gate)
  debugonce(.gate_stats)

  stimgate_gate(
    path_project = path_project,
    .data = gs,
    batch_list = batch_list,
    marker = marker_vec,
    debug = TRUE,
    tol_clust = 1e-6,
    min_cell = 40
  )
})
