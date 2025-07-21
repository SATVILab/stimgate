library(testthat)

test_that("stimgate_gate runs", {
  example_data <- get_gatingset_example()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))
  expect_true(file.exists(file.path(path_project, "gate_stats.rds")))
})

test_that("stimgate_plot function exists", {
  # Just test that the function exists and is callable
  expect_true(exists("stimgate_plot", envir = asNamespace("stimgate")))
  expect_true(is.function(stimgate::stimgate_plot))
})
