library(testthat)

test_that("stimgate_gate_get returns gate table after stimgate_gate", {
  # Skip if we can't load the required package
  skip_if_not_installed("stimgate")

  # Load the stimgate package
  library(stimgate)

  # Get example data
  exampleData <- get_example_data()
  gs <- flowWorkspace::load_gs(exampleData$path_gs)
  pathProject <- file.path(
    dirname(exampleData$path_gs),
    "stimgate_gate_get_test"
  )

  # Run stimgate_gate first to create the project structure
  invisible(stimgate_gate(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  ))

  # Test that stimgate_gate_get function works
  gateTbl <- stimgate_gate_get(pathProject)

  # Verify the gate table has expected structure
  expect_true(is.data.frame(gateTbl))
  expect_true(nrow(gateTbl) > 0)

  # Verify expected columns exist
  expectedCols <- c("chnl", "marker", "batch", "ind", "gate")
  expect_true(all(expectedCols %in% colnames(gateTbl)))

  # Verify that gate column contains numeric values
  expect_true(is.numeric(gateTbl$gate))

  # Verify that the gate table contains data for the expected markers
  expect_true(all(exampleData$chnl %in% gateTbl$chnl))

  # Verify that we have gate data for multiple samples
  expect_true(length(unique(gateTbl$ind)) > 1)

  # Clean up
  if (dir.exists(pathProject)) {
    unlink(pathProject, recursive = TRUE)
  }
  if (dir.exists(exampleData$path_gs)) {
    unlink(exampleData$path_gs, recursive = TRUE)
  }
})
