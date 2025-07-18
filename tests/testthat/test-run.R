library(testthat)

test_that("basic test works", {
  expect_true(TRUE)
})

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
  marker <- names(chnl_list)
  result <- stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = batch_list,
    marker = marker
  )
  expect_type(result, "character")
})

test_that("stimgate_gate runs with debug = TRUE", {
  fs <- get_fs()
  chnl_list <- get_chnl_list(fs = fs)
  batch_list <- chnl_list[[1]]$batch_list
  fs_gate <- chnl_list[[length(chnl_list)]]$fs
  dir_cache <- file.path(tempdir(), "stimgate_gate_debug")
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
  marker <- names(chnl_list)
  result <- stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = batch_list,
    marker = marker,
    debug = TRUE
  )
  expect_type(result, "character")
})
