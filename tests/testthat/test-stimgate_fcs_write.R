library(testthat)

test_that("stimgate_fcs_write runs with complete workflow", {
  skip_on_cran()
  
  # Setup test data
  fs <- get_fs()
  chnl_list <- get_chnl_list(fs = fs)
  batch_list <- chnl_list[[1]]$batch_list
  fs_gate <- chnl_list[[length(chnl_list)]]$fs
  dir_cache <- file.path(tempdir(), "stimgate_fcs_write_test")
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
  
  # Create gates first using stimgate_gate
  path_project <- file.path(dir_cache, "stimgate")
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = batch_list,
    marker = names(chnl_list)
  ))
  
  # Create a simple gate table for testing
  chnl <- names(chnl_list)
  simple_gate_tbl <- data.frame(
    chnl = rep(chnl, each = length(gs)),
    marker = rep(chnl, each = length(gs)), 
    batch = rep(1, length(chnl) * length(gs)),
    ind = rep(seq_along(gs), length(chnl)),
    gate = rep(1.0, length(chnl) * length(gs)),
    gate_name = rep("gate", length(chnl) * length(gs)),
    stringsAsFactors = FALSE
  )
  
  # Setup for FCS writing
  path_dir_save <- file.path(dir_cache, "fcs_output")
  
  # Create ind_batch_list from batch_list  
  ind_batch_list <- batch_list
  
  # Test stimgate_fcs_write with simple gate table
  result <- stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = ind_batch_list,
    path_dir_save = path_dir_save,
    chnl = chnl,
    gate_tbl = simple_gate_tbl
  )
  
  # Verify the function ran successfully
  expect_true(is.character(result))
  expect_equal(result, path_dir_save)
  
  # Verify directory was created
  expect_true(dir.exists(path_dir_save))
  
  # Clean up
  unlink(dir_cache, recursive = TRUE)
})

test_that("stimgate_fcs_write handles missing gate table gracefully", {
  # Setup minimal test data
  fs <- get_fs()
  chnl_list <- get_chnl_list(fs = fs)
  batch_list <- chnl_list[[1]]$batch_list
  fs_gate <- chnl_list[[length(chnl_list)]]$fs
  dir_cache <- file.path(tempdir(), "stimgate_fcs_write_error_test")
  if (dir.exists(dir_cache)) {
    unlink(dir_cache, recursive = TRUE)
  }
  dir.create(dir_cache, recursive = TRUE)
  
  # Create GatingSet but don't run gating
  path_gs <- get_gatingset(
    fs = fs_gate,
    dir_cache = dir_cache
  )
  gs <- flowWorkspace::load_gs(path_gs)
  
  path_project <- file.path(dir_cache, "stimgate_missing")
  path_dir_save <- file.path(dir_cache, "fcs_output_missing")
  chnl <- names(chnl_list)
  
  # Test should fail gracefully when gate files don't exist
  expect_error(
    stimgate::stimgate_fcs_write(
      path_project = path_project,
      .data = gs,
      ind_batch_list = batch_list,
      path_dir_save = path_dir_save,
      chnl = chnl
    ),
    "cannot open"
  )
  
  # Clean up
  unlink(dir_cache, recursive = TRUE)
})