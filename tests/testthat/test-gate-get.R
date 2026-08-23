test_that("getStimGates and getStimStats return gate table and stats after gateStim", {
  # Skip if we can't load the required package
  skip_if_not_installed("stimgate")

  # Get example data
  exampleData <- getExampleData()
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

  # Lightweight integration check: getStimStats returns a non-empty data frame
  statsTbl <- getStimStats(pathProject)
  expect_true(is.data.frame(statsTbl) && nrow(statsTbl) > 0)

  # Clean up
  if (dir.exists(pathProject)) {
    unlink(pathProject, recursive = TRUE)
  }
  if (dir.exists(exampleData$pathGs)) {
    unlink(exampleData$pathGs, recursive = TRUE)
  }
})
