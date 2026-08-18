test_that("getStimGates and getStimStats return gate table and stats after gateStim", {
  # Skip if we can't load the required package
  skip_if_not_installed("stimgate")

  # Get example data
  exampleData <- stimgate:::.getTestFixture()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(
    dirname(exampleData$pathGs),
    "stimgate_gate_get_test"
  )

  # Run gateStim first to create the project structure
  invisible(gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  ))

  # Test that getStimGates function works
  gateTbl <- getStimGates(pathProject)

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

  # Test that getStimStats function works
  statsTbl <- getStimStats(pathProject)

  # Verify the stats table has expected structure
  expect_true(is.data.frame(statsTbl))
  expect_true(nrow(statsTbl) > 0)

  # Verify error handling for getStimStats with invalid project path
  tmpEmptyDir <- file.path(tempdir(), "stimgate_empty_test_dir")
  dir.create(tmpEmptyDir, showWarnings = FALSE)
  expect_error(getStimStats(tmpEmptyDir), "No stats file found")
  unlink(tmpEmptyDir, recursive = TRUE)

  # Clean up
  if (dir.exists(pathProject)) {
    unlink(pathProject, recursive = TRUE)
  }
  if (dir.exists(exampleData$pathGs)) {
    unlink(exampleData$pathGs, recursive = TRUE)
  }
})
