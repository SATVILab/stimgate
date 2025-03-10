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
  browser()
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
  
  marker_vec <- c("Er168Di", "Lu175Di")

  browser()
  browser()
  browser()
  # debugonce(.complete_marker_list)
  # debugonce(.gate_init)
  # debugonce(.gate_batch)
  debugonce(.get_cp_uns_loc_sample)
  stimgate_gate(
    path_project = path_project,
    .data = gs,
    batch_list = batch_list,
    marker = marker_vec,
    debug = TRUE,
    min_cell = 40
  )
})
