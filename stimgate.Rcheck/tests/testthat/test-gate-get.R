library(testthat)

test_that("stimgate_gate_get returns gate table after stimgate_gate", {
  # Skip if we can't load the required package
  skip_if_not_installed("stimgate")
  
  # Load the stimgate package
  library(stimgate)
  
  # Get example data
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate_gate_get_test")
  
  # Run stimgate_gate first to create the project structure
  invisible(stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))
  
  # Test that stimgate_gate_get function works
  gate_tbl <- stimgate_gate_get(path_project)
  
  # Verify the gate table has expected structure
  expect_true(is.data.frame(gate_tbl))
  expect_true(nrow(gate_tbl) > 0)
  
  # Verify expected columns exist
  expected_cols <- c("chnl", "marker", "batch", "ind", "gate")
  expect_true(all(expected_cols %in% colnames(gate_tbl)))
  
  # Verify that gate column contains numeric values
  expect_true(is.numeric(gate_tbl$gate))
  
  # Verify that the gate table contains data for the expected markers
  expect_true(all(example_data$marker %in% gate_tbl$chnl))
  
  # Verify that we have gate data for multiple samples
  expect_true(length(unique(gate_tbl$ind)) > 1)
  
  # Clean up
  if (dir.exists(path_project)) {
    unlink(path_project, recursive = TRUE)
  }
  if (dir.exists(example_data$path_gs)) {
    unlink(example_data$path_gs, recursive = TRUE)
  }
})
