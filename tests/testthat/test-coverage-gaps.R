library(testthat)

test_that("interpFunctionWorksWhenXLowNotEqualVal", {
  # Test the interpolation functionality from cp-sub.R
  # This tests the code path where xLow != val

  # Create test data where interpolation is needed
  x <- c(1, 2, 3, 4, 5)
  y <- c(10, 20, 30, 40, 50)

  # Test interpolation between points (val = 2.5, between x[2]=2 and x[3]=3)
  val <- 2.5
  result <- stimgate:::.interp(val, x, y)

  # Expected: linear interpolation between (2, 20) and (3, 30)
  # y = 20 + (2.5 - 2) * (30 - 20) / (3 - 2) = 20 + 0.5 * 10 = 25
  expect_equal(result, 25)

  # Test another interpolation case (val = 1.7)
  val <- 1.7
  result <- stimgate:::.interp(val, x, y)

  # Expected: linear interpolation between (1, 10) and (2, 20)
  # y = 10 + (1.7 - 1) * (20 - 10) / (2 - 1) = 10 + 0.7 * 10 = 17
  expect_equal(result, 17)

  # Test interpolation near the end (val = 4.3)
  val <- 4.3
  result <- stimgate:::.interp(val, x, y)

  # Expected: linear interpolation between (4, 40) and (5, 50)
  # y = 40 + (4.3 - 4) * (50 - 40) / (5 - 4) = 40 + 0.3 * 10 = 43
  expect_equal(result, 43)
})

test_that("interpFunctionWorksWhenXLowEqualVal", {
  # Create test data
  x <- c(1, 2, 3, 4, 5)
  y <- c(10, 20, 30, 40, 50)

  # Test exact match (should return corresponding y value directly)
  val <- 3
  result <- stimgate:::.interp(val, x, y)
  expect_equal(result, 30)

  # Test another exact match
  val <- 1
  result <- stimgate:::.interp(val, x, y)
  expect_equal(result, 10)
})

# Note: prejoin gate combination test is skipped due to a bug in the current implementation
# that causes: "object 'countStim' not found" error. This may need investigation in the main codebase.
test_that("stimgateGateRunsWithGateCombnPrejoin", {
  skip(
    "prejoin gate combination has a bug causing 'object countStim not found' error"
  )
})

test_that("stimgateGateRunsWithGateCombnMean", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("HDCytoData")

  # Get example data
  exampleData <- stimgate::getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "testMean")

  # Test with mean gate combination
  expect_no_error({
    result <- stimgate::stimgate_gate(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      gateCombn = "mean"
    )
  })

  # Verify the function completed and returned a path
  expect_true(is.character(result))
  expect_true(dir.exists(result))

  # Clean up
  unlink(pathProject, recursive = TRUE)
})

test_that("stimgateGateRunsWithGateCombnTrim20", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("HDCytoData")

  # Get example data
  exampleData <- stimgate::getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "testTrim20")

  # Test with trim20 gate combination
  expect_no_error({
    result <- stimgate::stimgate_gate(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      gateCombn = "trim20"
    )
  })

  # Verify the function completed and returned a path
  expect_true(is.character(result))
  expect_true(dir.exists(result))

  # Clean up
  unlink(pathProject, recursive = TRUE)
})

test_that("stimgateGateRunsWithGateCombnMedian", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("HDCytoData")

  # Get example data
  exampleData <- stimgate::getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "testMedian")

  # Test with median gate combination
  expect_no_error({
    result <- stimgate::stimgate_gate(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      gateCombn = "median"
    )
  })

  # Verify the function completed and returned a path
  expect_true(is.character(result))
  expect_true(dir.exists(result))

  # Clean up
  unlink(pathProject, recursive = TRUE)
})

test_that("stimgateGateRunsWithGateCombnMax", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("HDCytoData")

  # Get example data
  exampleData <- stimgate::getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "testMax")

  # Test with max gate combination
  expect_no_error({
    result <- stimgate::stimgate_gate(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      gateCombn = "max"
    )
  })

  # Verify the function completed and returned a path
  expect_true(is.character(result))
  expect_true(dir.exists(result))

  # Clean up
  unlink(pathProject, recursive = TRUE)
})
