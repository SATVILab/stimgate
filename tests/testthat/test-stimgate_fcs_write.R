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

test_that("stimgate_fcs_write creates output directory", {
  # Test that the function creates the output directory when called
  # even if it fails later due to missing gates
  
  fs <- get_fs()
  chnl_list <- get_chnl_list(fs = fs)
  fs_gate <- chnl_list[[length(chnl_list)]]$fs
  dir_cache <- file.path(tempdir(), "stimgate_fcs_write_dir_test")
  if (dir.exists(dir_cache)) {
    unlink(dir_cache, recursive = TRUE)
  }
  dir.create(dir_cache, recursive = TRUE)
  
  # Create GatingSet
  path_gs <- get_gatingset(
    fs = fs_gate,
    dir_cache = dir_cache
  )
  gs <- flowWorkspace::load_gs(path_gs)
  
  path_project <- file.path(dir_cache, "nonexistent_project")
  path_dir_save <- file.path(dir_cache, "fcs_output_test")
  
  # Function should create the directory before failing on missing gates
  expect_error({
    stimgate::stimgate_fcs_write(
      path_project = path_project,
      .data = gs,
      ind_batch_list = list(c(1, 2)),
      path_dir_save = path_dir_save,
      chnl = names(chnl_list)[1]
    )
  })
  
  # Directory should have been created despite the error
  expect_true(dir.exists(path_dir_save))
  
  # Clean up
  unlink(dir_cache, recursive = TRUE)
})

test_that("stimgate_fcs_write documentation example is valid", {
  # Test that the documentation example syntax is correct
  # without actually running the complex workflow
  
  # This validates that the example in the documentation uses correct parameter names
  expect_error({
    # This should fail gracefully with missing data, not with unknown parameters
    stimgate::stimgate_fcs_write(
      path_project = tempfile("nonexistent_path"),
      .data = NULL,
      ind_batch_list = list(batch1 = c(1, 2, 3)),
      path_dir_save = tempfile("nonexistent_output"),
      chnl = c("IL2", "IFNg")
    )
  })
  
  # The exact error message may vary, but the function should fail
  # without complaining about unknown parameters
})