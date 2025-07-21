library(testthat)

library(testthat)

test_that(".interp function works when x_low != val (interpolation case)", {
  # Test the interpolation functionality from cp-sub.R
  # This tests the code path where x_low != val (line 35-42 in cp-sub.R)
  
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

test_that(".interp function works when x_low == val (exact match case)", {
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

test_that("stimgate_gate runs with gate_combn = 'prejoin'", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("HDCytoData")
  
  # Get example data
  example_data <- stimgate::get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(tempdir(), "test_prejoin")
  
  # Test with prejoin gate combination
  expect_no_error({
    result <- stimgate::stimgate_gate(
      .data = gs,
      path_project = path_project,
      pop_gate = "root",
      batch_list = example_data$batch_list,
      marker = example_data$marker,
      gate_combn = "prejoin",
      debug = FALSE
    )
  })
  
  # Verify the function completed and returned a path
  expect_true(is.character(result))
  expect_true(dir.exists(result))
  
  # Clean up
  unlink(path_project, recursive = TRUE)
})

test_that("stimgate_gate runs with gate_combn = 'mean'", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("HDCytoData")
  
  # Get example data
  example_data <- stimgate::get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(tempdir(), "test_mean")
  
  # Test with mean gate combination
  expect_no_error({
    result <- stimgate::stimgate_gate(
      .data = gs,
      path_project = path_project,
      pop_gate = "root",
      batch_list = example_data$batch_list,
      marker = example_data$marker,
      gate_combn = "mean",
      debug = FALSE
    )
  })
  
  # Verify the function completed and returned a path
  expect_true(is.character(result))
  expect_true(dir.exists(result))
  
  # Clean up
  unlink(path_project, recursive = TRUE)
})

test_that("stimgate_gate runs with gate_combn = 'trim20'", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("HDCytoData")
  
  # Get example data
  example_data <- stimgate::get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(tempdir(), "test_trim20")
  
  # Test with trim20 gate combination
  expect_no_error({
    result <- stimgate::stimgate_gate(
      .data = gs,
      path_project = path_project,
      pop_gate = "root",
      batch_list = example_data$batch_list,
      marker = example_data$marker,
      gate_combn = "trim20",
      debug = FALSE
    )
  })
  
  # Verify the function completed and returned a path
  expect_true(is.character(result))
  expect_true(dir.exists(result))
  
  # Clean up
  unlink(path_project, recursive = TRUE)
})

test_that("stimgate_gate runs with gate_combn = 'median'", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("HDCytoData")
  
  # Get example data
  example_data <- stimgate::get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(tempdir(), "test_median")
  
  # Test with median gate combination
  expect_no_error({
    result <- stimgate::stimgate_gate(
      .data = gs,
      path_project = path_project,
      pop_gate = "root",
      batch_list = example_data$batch_list,
      marker = example_data$marker,
      gate_combn = "median",
      debug = FALSE
    )
  })
  
  # Verify the function completed and returned a path
  expect_true(is.character(result))
  expect_true(dir.exists(result))
  
  # Clean up
  unlink(path_project, recursive = TRUE)
})

test_that("stimgate_gate runs with gate_combn = 'max'", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("HDCytoData")
  
  # Get example data
  example_data <- stimgate::get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(tempdir(), "test_max")
  
  # Test with max gate combination
  expect_no_error({
    result <- stimgate::stimgate_gate(
      .data = gs,
      path_project = path_project,
      pop_gate = "root",
      batch_list = example_data$batch_list,
      marker = example_data$marker,
      gate_combn = "max",
      debug = FALSE
    )
  })
  
  # Verify the function completed and returned a path
  expect_true(is.character(result))
  expect_true(dir.exists(result))
  
  # Clean up
  unlink(path_project, recursive = TRUE)
})