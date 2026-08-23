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

test_that("stimgateGateRunsWithGateCombnPrejoin", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  batchList <- list(batch1 = c(1, 2, 4))
  stimInd <- as.character(batchList[[1]][-1])
  pathProject <- tempfile("testPrejoin")
  withr::defer(unlink(pathProject, recursive = TRUE))
  withr::local_envvar(STIMGATE_INTERMEDIATE = "all")

  result <- gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = batchList,
    marker = exampleData$marker,
    gateCombn = "prejoin"
  )

  expect_identical(result, pathProject)
  expect_identical(stimgateMetaReadBatchList(pathProject), batchList)

  gateTbl <- getStimGates(pathProject)
  stimGateTbl <- gateTbl[gateTbl$ind %in% stimInd, ]
  expect_gt(nrow(stimGateTbl), 0L)
  expect_true(all(grepl("prejoin", stimGateTbl$gateName, fixed = TRUE)))

  gateGroups <- split(
    stimGateTbl,
    interaction(stimGateTbl$gateName, stimGateTbl$chnl, drop = TRUE)
  )
  expect_true(all(vapply(gateGroups, function(x) {
    setequal(x$ind, stimInd) && length(unique(x$gate)) == 1L
  }, logical(1))))

  detailTbl <- getStimGatesDetailed(pathProject)
  stimDetailTbl <- detailTbl[
    detailTbl$detailLevel == "sample" & detailTbl$ind %in% stimInd,
  ]
  expect_gt(nrow(stimDetailTbl), 0L)
  detailGroups <- split(stimDetailTbl, stimDetailTbl$chnl)
  expect_true(all(vapply(detailGroups, function(x) {
    setequal(x$ind, stimInd) && length(unique(x$threshold)) == 1L
  }, logical(1))))
  expect_true(all(stimDetailTbl$locGenerated))
  expect_false(any(stimDetailTbl$locGeneratedDirect))
  expect_true(all(stimDetailTbl$locSource == "prejoin"))
  expect_true(all(
    stimDetailTbl$thresholdOrigin ==
      "prejoin_generated_from_joined_stim_conditions"
  ))
})

test_that("stimgateGateRunsWithGateCombnNo", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  # Get example data
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "testNo")

  # Test with no gate combination
  result <- expect_no_error({
    stimgate::gateStim(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      gateCombn = "no"
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

  # Get example data
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "testMedian")

  # Test with median gate combination
  result <- expect_no_error({
    stimgate::gateStim(
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

  # Get example data
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "testMax")

  # Test with max gate combination
  result <- expect_no_error({
    stimgate::gateStim(
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

test_that("completeChnlSettingsBiasUns defaults to 1/4 of bwFallback", {
  # Default case: biasUns is NULL, biasUnsFactor is 1, bwFallback is 0.4
  expect_equal(
    stimgate:::.completeChnlSettingsBiasUns(
      biasUns = NULL,
      biasUnsFactor = 1,
      bwMin = 0.1,
      bwMax = 1.0,
      bwFallback = 0.4
    ),
    0.1
  )

  # Scaled by biasUnsFactor: biasUnsFactor is 2, bwFallback is 0.4
  expect_equal(
    stimgate:::.completeChnlSettingsBiasUns(
      biasUns = NULL,
      biasUnsFactor = 2,
      bwMin = 0.1,
      bwMax = 1.0,
      bwFallback = 0.4
    ),
    0.2
  )

  # Explicit biasUns overrides bwFallback
  expect_equal(
    stimgate:::.completeChnlSettingsBiasUns(
      biasUns = 0.5,
      biasUnsFactor = 1,
      bwMin = 0.1,
      bwMax = 1.0,
      bwFallback = 0.4
    ),
    0.5
  )

  # bwFallback is NULL: falls back to 1/4 of mean(bwMin, bwMax)
  expect_equal(
    stimgate:::.completeChnlSettingsBiasUns(
      biasUns = NULL,
      biasUnsFactor = 1,
      bwMin = 0.2,
      bwMax = 0.6,
      bwFallback = NULL
    ),
    0.25 * 0.4
  )

  # bwFallback is NULL and no valid bw limits: returns 0
  expect_equal(
    stimgate:::.completeChnlSettingsBiasUns(
      biasUns = NULL,
      biasUnsFactor = 1,
      bwMin = -Inf,
      bwMax = Inf,
      bwFallback = NULL
    ),
    0
  )
})

test_that("gateStim defaults biasUns to 1/4 of bwFallback in metadata", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "testBiasUnsDefault")

  result <- gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  )

  chnlSettings <- stimgateMetaReadSettingsChnls(pathProject)
  for (chnlName in names(chnlSettings)) {
    bwFallback <- chnlSettings[[chnlName]]$bwFallback
    biasUns <- chnlSettings[[chnlName]]$biasUns
    expect_true(is.numeric(bwFallback) && bwFallback > 0)
    expect_equal(biasUns, 0.25 * bwFallback)
  }

  unlink(pathProject, recursive = TRUE)
})
