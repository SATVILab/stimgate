library(testthat)
library(flowCore)

test_that("chnlLabWorksWithFlowFrameObjects", {
  # Load test data
  data(GvHD, package = "flowCore")
  ff <- GvHD[[1]]

  # Test the function
  result <- .chnlLab(ff)

  # Check that result is a character object (may have AsIs class)
  expect_type(result, "character")
  expect_true(!is.null(names(result)))

  # Check that names are channel names
  paramData <- flowCore::parameters(ff)@data
  expect_equal(names(result), as.character(paramData$name))

  # Check that values are marker descriptions (with NA handling)
  expectedValues <- as.character(paramData$desc)
  for (i in seq_along(expectedValues)) {
    if (is.na(expectedValues[i])) {
      expectedValues[i] <- as.character(paramData$name[i])
    }
  }
  expect_equal(as.character(result), expectedValues)
})

test_that("chnlLabWorksWithFlowSetObjects", {
  # Load test data
  data(GvHD, package = "flowCore")
  fs <- GvHD[1:2]

  # Test the function
  result <- .chnlLab(fs)

  # Check that result is a character object (may have AsIs class)
  expect_type(result, "character")
  expect_true(!is.null(names(result)))

  # Check that it uses the first flowFrame's parameters
  paramData <- flowCore::parameters(fs[[1]])@data
  expect_equal(names(result), as.character(paramData$name))

  # Check that values are marker descriptions (with NA handling)
  expectedValues <- as.character(paramData$desc)
  for (i in seq_along(expectedValues)) {
    if (is.na(expectedValues[i])) {
      expectedValues[i] <- as.character(paramData$name[i])
    }
  }
  expect_equal(as.character(result), expectedValues)
})

test_that("chnlLabHandlesNaMarkerDescriptionsCorrectly", {
  # Load test data
  data(GvHD, package = "flowCore")
  ff <- GvHD[[1]]

  # Modify the flowFrame to have some NA descriptions
  paramData <- flowCore::parameters(ff)@data
  originalDesc <- as.character(paramData$desc)

  # Set some descriptions to NA
  paramData$desc[1:2] <- NA
  flowCore::parameters(ff)@data <- paramData

  # Test the function
  result <- .chnlLab(ff)

  # Check that NA descriptions were replaced with channel names
  expect_equal(as.character(result[1]), as.character(paramData$name[1]))
  expect_equal(as.character(result[2]), as.character(paramData$name[2]))
  expect_equal(names(result)[1], as.character(paramData$name[1]))
  expect_equal(names(result)[2], as.character(paramData$name[2]))

  # Check that non-NA descriptions remain unchanged
  if (length(paramData$desc) > 2 && !is.na(originalDesc[3])) {
    expect_equal(as.character(result[3]), originalDesc[3])
  }
})

test_that("chnlLabThrowsErrorForUnsupportedObjectTypes", {
  # Test with various unsupported objects
  expect_error(
    .chnlLab(data.frame(x = 1:5, y = 6:10)),
    "classOfDataNotRecognised"
  )

  expect_error(
    .chnlLab(matrix(1:10, nrow = 2)),
    "classOfDataNotRecognised"
  )

  expect_error(
    .chnlLab(list(a = 1, b = 2)),
    "classOfDataNotRecognised"
  )

  expect_error(
    .chnlLab("characterString"),
    "classOfDataNotRecognised"
  )

  expect_error(
    .chnlLab(123),
    "classOfDataNotRecognised"
  )
})

test_that("chnlLabResultStructureIsConsistent", {
  # Load test data
  data(GvHD, package = "flowCore")
  ff <- GvHD[[1]]

  result <- .chnlLab(ff)

  # Check that every element has a name
  expect_true(all(nzchar(names(result))))
  expect_equal(length(result), length(names(result)))

  # Check that result length matches number of parameters
  paramData <- flowCore::parameters(ff)@data
  expect_equal(length(result), nrow(paramData))

  # Check that no values are NULL or missing
  expect_false(any(is.null(result)))
  expect_false(any(is.na(as.character(result))))
})

test_that("chnlLabHandlesEdgeCaseWithAllNaDescriptions", {
  # Load test data
  data(GvHD, package = "flowCore")
  ff <- GvHD[[1]]

  # Set all descriptions to NA
  paramData <- flowCore::parameters(ff)@data
  paramData$desc <- rep(NA, nrow(paramData))
  flowCore::parameters(ff)@data <- paramData

  # Test the function
  result <- .chnlLab(ff)

  # Check that all values are channel names (since all descriptions were NA)
  expect_equal(as.character(result), as.character(paramData$name))
  expect_equal(names(result), as.character(paramData$name))
})

test_that("chnlLabHandlesEdgeCaseWithNoNaDescriptions", {
  # Load test data
  data(GvHD, package = "flowCore")
  ff <- GvHD[[1]]

  # Ensure no descriptions are NA by replacing any existing NAs
  paramData <- flowCore::parameters(ff)@data
  originalDesc <- as.character(paramData$desc)

  # Replace any NA with a dummy description
  for (i in seq_along(paramData$desc)) {
    if (is.na(paramData$desc[i])) {
      paramData$desc[i] <- paste0("marker_", i)
    }
  }
  flowCore::parameters(ff)@data <- paramData

  # Test the function
  result <- .chnlLab(ff)

  # Check that all values are the actual descriptions (no channel names used)
  expect_equal(as.character(result), as.character(paramData$desc))
  expect_equal(names(result), as.character(paramData$name))
  expect_false(any(as.character(result) == names(result))) # No fallback to channel names
})
