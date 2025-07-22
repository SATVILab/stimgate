library(testthat)
library(flowCore)

# Source the function to test
source("/home/runner/work/stimgate/stimgate/R/UtilsCytoRSV-chnl_lab.R")

test_that("chnl_lab works with flowFrame objects", {
  # Load test data
  data(GvHD, package = "flowCore")
  ff <- GvHD[[1]]
  
  # Test the function
  result <- chnl_lab(ff)
  
  # Check that result is a character object (may have AsIs class)
  expect_type(result, "character")
  expect_true(!is.null(names(result)))
  
  # Check that names are channel names
  param_data <- flowCore::parameters(ff)@data
  expect_equal(names(result), as.character(param_data$name))
  
  # Check that values are marker descriptions (with NA handling)
  expected_values <- as.character(param_data$desc)
  for (i in seq_along(expected_values)) {
    if (is.na(expected_values[i])) {
      expected_values[i] <- as.character(param_data$name[i])
    }
  }
  expect_equal(as.character(result), expected_values)
})

test_that("chnl_lab works with flowSet objects", {
  # Load test data
  data(GvHD, package = "flowCore")
  fs <- GvHD[1:2]
  
  # Test the function
  result <- chnl_lab(fs)
  
  # Check that result is a character object (may have AsIs class)
  expect_type(result, "character")
  expect_true(!is.null(names(result)))
  
  # Check that it uses the first flowFrame's parameters
  param_data <- flowCore::parameters(fs[[1]])@data
  expect_equal(names(result), as.character(param_data$name))
  
  # Check that values are marker descriptions (with NA handling)
  expected_values <- as.character(param_data$desc)
  for (i in seq_along(expected_values)) {
    if (is.na(expected_values[i])) {
      expected_values[i] <- as.character(param_data$name[i])
    }
  }
  expect_equal(as.character(result), expected_values)
})

test_that("chnl_lab handles NA marker descriptions correctly", {
  # Load test data
  data(GvHD, package = "flowCore")
  ff <- GvHD[[1]]
  
  # Modify the flowFrame to have some NA descriptions
  param_data <- flowCore::parameters(ff)@data
  original_desc <- as.character(param_data$desc)
  
  # Set some descriptions to NA
  param_data$desc[1:2] <- NA
  flowCore::parameters(ff)@data <- param_data
  
  # Test the function
  result <- chnl_lab(ff)
  
  # Check that NA descriptions were replaced with channel names
  expect_equal(as.character(result[1]), as.character(param_data$name[1]))
  expect_equal(as.character(result[2]), as.character(param_data$name[2]))
  expect_equal(names(result)[1], as.character(param_data$name[1]))
  expect_equal(names(result)[2], as.character(param_data$name[2]))
  
  # Check that non-NA descriptions remain unchanged
  if (length(param_data$desc) > 2 && !is.na(original_desc[3])) {
    expect_equal(as.character(result[3]), original_desc[3])
  }
})

test_that("chnl_lab throws error for unsupported object types", {
  # Test with various unsupported objects
  expect_error(
    chnl_lab(data.frame(x = 1:5, y = 6:10)),
    "class of data not recognised"
  )
  
  expect_error(
    chnl_lab(matrix(1:10, nrow = 2)),
    "class of data not recognised"
  )
  
  expect_error(
    chnl_lab(list(a = 1, b = 2)),
    "class of data not recognised"
  )
  
  expect_error(
    chnl_lab("character_string"),
    "class of data not recognised"
  )
  
  expect_error(
    chnl_lab(123),
    "class of data not recognised"
  )
})

test_that("chnl_lab result structure is consistent", {
  # Load test data
  data(GvHD, package = "flowCore")
  ff <- GvHD[[1]]
  
  result <- chnl_lab(ff)
  
  # Check that every element has a name
  expect_true(all(nzchar(names(result))))
  expect_equal(length(result), length(names(result)))
  
  # Check that result length matches number of parameters
  param_data <- flowCore::parameters(ff)@data
  expect_equal(length(result), nrow(param_data))
  
  # Check that no values are NULL or missing
  expect_false(any(is.null(result)))
  expect_false(any(is.na(as.character(result))))
})

test_that("chnl_lab handles edge case with all NA descriptions", {
  # Load test data
  data(GvHD, package = "flowCore")
  ff <- GvHD[[1]]
  
  # Set all descriptions to NA
  param_data <- flowCore::parameters(ff)@data
  param_data$desc <- rep(NA, nrow(param_data))
  flowCore::parameters(ff)@data <- param_data
  
  # Test the function
  result <- chnl_lab(ff)
  
  # Check that all values are channel names (since all descriptions were NA)
  expect_equal(as.character(result), as.character(param_data$name))
  expect_equal(names(result), as.character(param_data$name))
})

test_that("chnl_lab handles edge case with no NA descriptions", {
  # Load test data
  data(GvHD, package = "flowCore")
  ff <- GvHD[[1]]
  
  # Ensure no descriptions are NA by replacing any existing NAs
  param_data <- flowCore::parameters(ff)@data
  original_desc <- as.character(param_data$desc)
  
  # Replace any NA with a dummy description
  for (i in seq_along(param_data$desc)) {
    if (is.na(param_data$desc[i])) {
      param_data$desc[i] <- paste0("Marker_", i)
    }
  }
  flowCore::parameters(ff)@data <- param_data
  
  # Test the function
  result <- chnl_lab(ff)
  
  # Check that all values are the actual descriptions (no channel names used)
  expect_equal(as.character(result), as.character(param_data$desc))
  expect_equal(names(result), as.character(param_data$name))
  expect_false(any(as.character(result) == names(result)))  # No fallback to channel names
})