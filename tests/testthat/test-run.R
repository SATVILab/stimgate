library(testthat)
library(stimgate)

test_that("stimgate_gate runs", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")
  invisible(stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))
  expect_true(file.exists(file.path(path_project, "gate_stats.rds")))
})

test_that("stimgate_gate runs with debug = TRUE", {
  # Run stimgate_gate with debug = TRUE to ensure it works without error
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate_debug")
  
  # Test that the function can be called with debug = TRUE without error
  # and that it returns the expected path
  result_path <- stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker,
    debug = TRUE
  )
  
  # Verify the function returns the expected path
  expect_equal(result_path, path_project)

  # Verify basic output files still exist (same as non-debug version)
  expect_true(file.exists(file.path(path_project, "gate_stats.rds")))

  # Clean up debug test directory
  if (dir.exists(path_project)) {
    unlink(path_project, recursive = TRUE)
  }
})
