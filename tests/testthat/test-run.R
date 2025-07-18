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

test_that("stim_gate_plot function exists", {
  # Just test that the function exists and is callable
  expect_true(exists("stim_gate_plot", envir = asNamespace("stimgate")))
  expect_true(is.function(stimgate::stim_gate_plot))
})
