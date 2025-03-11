library(testthat)

test_that("stimgate_gate runs", {
  path_gs <- testthat::test_path(
    "testdata", "postmortem", "gs", "gs_cytof_acs_cd4"
  )
  gs <- flowWorkspace::load_gs(path_gs)
  fn_tbl_info <- get_fn_tbl_info_postmortem(gs)
  # debugonce(get_batch_list_postmortem)
  batch_list <- get_batch_list_postmortem(
    fn_tbl_info,
    col_grp = c("pid", "location"),
    col_stim = "stim",
    uns_chr = "uns",
    min_cell_uns = 100,
    min_cell_stim = 100,
    col_n_cell = "n_cell_pop"
  )
  batch_list <- batch_list |>
    stats::setNames(paste0("batch_", seq_along(batch_list)))

  batch_vec_sort <- lapply(batch_list, function(x) {
    fn_tbl_info[["n_cell_pop"]][x] |>
      min()
  }) |>
    stats::setNames(seq_along(batch_list)) |>
    unlist() |>
    sort() |>
    rev()
  batch_vec_sel <- batch_vec_sort[1:5] |>
    names() |>
    as.numeric()
  batch_list <- batch_list[batch_vec_sel]
  
  marker_vec <- c("Er168Di", "Lu175Di")
  path_project <- setup_project_postmortem()

  browser()
  browser()
  browser()
  debugonce(.get_cp_cluster_prop_bs_by_cp_tbl_obj)
  debugonce(stimgate_gate)
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
