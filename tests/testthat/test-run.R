library(testthat)

test_that("stimgate_gate runs", {
  fs <- get_fs()
  chnl_list <- get_chnl_list(fs = fs)
  batch_list <- chnl_list[[1]]$batch_list
  fs_gate <- chnl_list[[length(chnl_list)]]$fs
  dir_cache <- file.path(tempdir(), "stimgate_gate")
  if (dir.exists(dir_cache)) {
    unlink(dir_cache, recursive = TRUE)
  }
  dir.create(dir_cache, recursive = TRUE)
  path_gs <- get_gatingset(
    fs = fs_gate,
    dir_cache = dir_cache
  )
  gs <- flowWorkspace::load_gs(path_gs)
  path_project <- file.path(dir_cache, "stimgate")
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = batch_list,
    marker = names(chnl_list)
  ))
})

test_that("stimgate_gate_get returns gate table after stimgate_gate", {
  fs <- get_fs()
  chnl_list <- get_chnl_list(fs = fs)
  batch_list <- chnl_list[[1]]$batch_list
  fs_gate <- chnl_list[[length(chnl_list)]]$fs
  dir_cache <- file.path(tempdir(), "stimgate_gate_get")
  if (dir.exists(dir_cache)) {
    unlink(dir_cache, recursive = TRUE)
  }
  dir.create(dir_cache, recursive = TRUE)
  path_gs <- get_gatingset(
    fs = fs_gate,
    dir_cache = dir_cache
  )
  gs <- flowWorkspace::load_gs(path_gs)
  path_project <- file.path(dir_cache, "stimgate")
  
  # First run stimgate_gate to generate the data
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = batch_list,
    marker = names(chnl_list)
  ))
  
  # Then test stimgate_gate_get
  gate_tbl <- stimgate::stimgate_gate_get(path_project)
  
  # Verify the result is a data frame
  expect_s3_class(gate_tbl, "data.frame")
  
  # Verify it has the expected columns based on the function implementation
  expected_cols <- c("chnl")
  expect_true(all(expected_cols %in% colnames(gate_tbl)))
  
  # Verify it has data (rows)
  expect_gt(nrow(gate_tbl), 0)
})
