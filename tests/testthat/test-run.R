test_that("stimgateGateRuns", {
  exampleData <- getExampleData(nCell = 1e3)
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "stimgate")
  # debugonce(.getCpUnsLocGetProb)
  # debugonce(stimgate_gate)
  # browser()
  debugonce(.getCpCluster)
  Sys.setenv("stimgateIntermediate" = "true")
  invisible(gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  ))
  # browser()
  expect_true(file.exists(file.path(pathProject, "gateStats.rds")))
})

test_that("stimgateGateRunsWithStimgateDebugEnvironmentVariable", {
  skip()
  # Run stimgateGate with stimgateDebug environment variable to ensure it works without error
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "stimgateDebug")

  # Set the environment variable for debug mode
  Sys.setenv("stimgateDebug" = "true")
  on.exit(Sys.unsetenv("stimgateDebug"))

  # Test that the function can be called with debug enabled without error
  # and that it returns the expected path
  resultPath <- stimgateGate(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  )

  # Verify the function returns the expected path
  expect_equal(resultPath, pathProject)

  # Verify basic output files still exist (same as non-debug version)
  expect_true(file.exists(file.path(pathProject, "gateStats.rds")))

  # Clean up debug test directory
  if (dir.exists(pathProject)) {
    unlink(pathProject, recursive = TRUE)
  }
})

test_that("stimgateGateRunsWithStimgateIntermediateEnvironmentVariable", {
  skip()
  # Run stimgateGate with stimgateIntermediate environment variable to ensure
  # intermediate data is saved correctly
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(
    dirname(exampleData$pathGs),
    "stimgateIntermediate"
  )

  # Set the environment variable for intermediate data saving
  Sys.setenv("stimgateIntermediate" = "true")
  on.exit(Sys.unsetenv("stimgateIntermediate"))

  # Test that the function can be called with intermediate saving enabled
  resultPath <- stimgateGate(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  )

  # Verify the function returns the expected path
  expect_equal(resultPath, pathProject)

  # Verify basic output files exist
  expect_true(file.exists(file.path(pathProject, "gateStats.rds")))

  # Verify intermediate data directory exists
  intermediateDir <- file.path(pathProject, "intermediateData")
  expect_true(dir.exists(intermediateDir))

  # Verify that intermediate files were created
  intermediateFiles <- list.files(intermediateDir, recursive = TRUE)
  expect_true(length(intermediateFiles) > 0)

  # Verify that files exist for the "init" stage
  initFiles <- list.files(
    file.path(intermediateDir, "init"),
    recursive = TRUE
  )
  expect_true(length(initFiles) > 0)

  # Clean up intermediate test directory
  if (dir.exists(pathProject)) {
    unlink(pathProject, recursive = TRUE)
  }
})
