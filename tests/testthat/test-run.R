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
  batch_list <- batch_list[seq_len(5)]
  path_project <- file.path(tempdir(), "stimgate_gate_run")
  suppressWarnings(unlink(path_project, recursive = TRUE))
  dir.create(path_project, showWarnings = FALSE, recursive = TRUE)
  marker_vec <- c("Lu175Di", "Er168Di")
  debugonce(.gate_batch)
  debugonce(.gate_batch_all)
  stimgate_gate(
    path_project = path_project,
    .data = gs,
    batch_list = batch_list,
    marker = marker_vec,
    debug = TRUE
  )


})
