makeChnlLabTestFrame <- function(desc = c("CD3", "CD4")) {
  exprs <- matrix(
    seq_len(8),
    ncol = 2,
    dimnames = list(NULL, c("FSC-A", "FL1-H"))
  )
  ff <- flowCore::flowFrame(exprs)
  paramData <- flowCore::parameters(ff)@data
  paramData$desc <- desc
  flowCore::parameters(ff)@data <- paramData
  ff
}

test_that("chnlLabWorksWithFlowFrameObjects", {
  ff <- makeChnlLabTestFrame()

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
  fs <- flowCore::flowSet(list(
    makeChnlLabTestFrame(),
    makeChnlLabTestFrame(c("CD8", "CD19"))
  ))

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

test_that("chnlLabWorksWithGatingSetObjects", {
  fs <- flowCore::flowSet(list(
    makeChnlLabTestFrame(),
    makeChnlLabTestFrame(c("CD8", "CD19"))
  ))
  gs <- flowWorkspace::GatingSet(fs)

  result_dot <- .chnlLab(gs)
  result_pub <- chnlLab(gs)

  expect_type(result_dot, "character")
  expect_type(result_pub, "character")
  expect_equal(result_dot, result_pub)
  expect_true(!is.null(names(result_dot)))
})

test_that("chnlLabHandlesNaMarkerDescriptionsCorrectly", {
  ff <- makeChnlLabTestFrame()

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
  ff <- makeChnlLabTestFrame()

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
  ff <- makeChnlLabTestFrame(rep(NA_character_, 2))

  # Test the function
  result <- .chnlLab(ff)

  # Check that all values are channel names (since all descriptions were NA)
  paramData <- flowCore::parameters(ff)@data
  expect_equal(as.character(result), as.character(paramData$name))
  expect_equal(names(result), as.character(paramData$name))
})

test_that("chnlLabHandlesEdgeCaseWithNoNaDescriptions", {
  ff <- makeChnlLabTestFrame()
  paramData <- flowCore::parameters(ff)@data

  # Test the function
  result <- .chnlLab(ff)

  # Check that all values are the actual descriptions (no channel names used)
  expect_equal(as.character(result), as.character(paramData$desc))
  expect_equal(names(result), as.character(paramData$name))
  expect_false(any(as.character(result) == names(result))) # No fallback to channel names
})
