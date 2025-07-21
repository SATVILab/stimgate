library(testthat)

test_that("stimgate_fcs_write function exists and has correct signature", {
  # Test that the function exists and has the expected parameters
  expect_true(exists("stimgate_fcs_write", where = asNamespace("stimgate")))

  # Test function signature by checking for argument names
  args <- names(formals(stimgate::stimgate_fcs_write))
  expected_args <- c("path_project", ".data", "ind_batch_list", "path_dir_save",
                     "chnl", "gate_tbl", "trans_fn", "trans_chnl", "combn_exc",
                     "gate_type_cyt_pos", "gate_type_single_pos", "mult", "gate_uns_method")

  expect_true(all(expected_args %in% args))
})

test_that("stimgate_fcs_write handles invalid inputs appropriately", {
  # Test with NULL .data
  expect_error(
    stimgate::stimgate_fcs_write(
      path_project = tempfile(),
      .data = NULL,
      ind_batch_list = list(c(1, 2)),
      path_dir_save = tempfile()
    )
  )

  # Test with missing required parameters
  expect_error(
    stimgate::stimgate_fcs_write()
  )
})

test_that("stimgate_fcs_write runs", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))
  path_dir_save <- file.path(tempdir(), "fcs_output_test")

  # Function should create the directory before failing on missing gates
  # debugonce(.fcs_write_get_gate_tbl_add_uns)
  stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]]
  )
  browser()
  ))
})
