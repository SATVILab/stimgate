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
    control = stimControl(gateCombn = "prejoin")
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
      control = stimgate::stimControl(gateCombn = "no")
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
      control = stimgate::stimControl(gateCombn = "median")
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
      control = stimgate::stimControl(gateCombn = "max")
    )
  })

  # Verify the function completed and returned a path
  expect_true(is.character(result))
  expect_true(dir.exists(result))

  # Clean up
  unlink(pathProject, recursive = TRUE)
})

test_that("completeChnlSettingsBiasUns defaults to bwFallback", {
  # Default case: biasUns is NULL, biasUnsFactor is 1, bwFallback is 0.4
  expect_equal(
    stimgate:::.completeChnlSettingsBiasUns(
      biasUns = NULL,
      biasUnsFactor = 1,
      bwMin = 0.1,
      bwMax = 1.0,
      bwFallback = 0.4
    ),
    0.4
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
    0.8
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

  # bwFallback is NULL: falls back to mean(bwMin, bwMax)
  expect_equal(
    stimgate:::.completeChnlSettingsBiasUns(
      biasUns = NULL,
      biasUnsFactor = 1,
      bwMin = 0.2,
      bwMax = 0.6,
      bwFallback = NULL
    ),
    0.4
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

test_that("completeChnlSettingsBiasUns prefers the common bandwidth", {
  biasUns <- function(...) {
    stimgate:::.completeChnlSettingsBiasUns(
      biasUns = NULL, bwMin = 0.1, bwMax = 1, bwFallback = 0.4, ...
    )
  }
  expect_equal(biasUns(biasUnsFactor = 1, bwCommon = 0.25), 0.25)
  expect_equal(biasUns(biasUnsFactor = 2, bwCommon = 0.25), 0.5)
  expect_equal(biasUns(biasUnsFactor = 1, bwCommon = NULL), 0.4)

  common <- stimgate:::.completeChnlSettingsBwCommon
  # A fixed bandwidth, or a bandwidth shared by the whole channel.
  expect_equal(common(list(bw = 0.3, bwShared = 0.2)), 0.3)
  expect_equal(common(list(bwShared = 0.2)), 0.2)
  # Per-cluster, per-sample and missing bandwidths have no common value.
  expect_null(common(list(
    bwShared = 0.2,
    bwSharedTbl = data.frame(ind = c("1", "2"), bw = c(0.1, 0.3))
  )))
  expect_null(common(list()))
  expect_null(common(list(bwShared = NA_real_)))
})

test_that("gateStim defaults biasUns to the shared bandwidth in metadata", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- withr::local_tempdir("testBiasUnsDefault")

  gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  )

  chnlSettings <- stimgateMetaReadSettingsChnls(pathProject)
  for (chnlName in names(chnlSettings)) {
    bwShared <- chnlSettings[[chnlName]]$bwShared
    expect_true(is.numeric(bwShared) && bwShared > 0)
    expect_equal(chnlSettings[[chnlName]]$biasUns, bwShared)
  }

  # Per-sample bandwidths have no common value: the fallback is used.
  pathSample <- withr::local_tempdir("testBiasUnsSample")
  gateStim(
    .data = gs,
    pathProject = pathSample,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker,
    control = stimControl(bwScope = "sample")
  )
  chnlSettings <- stimgateMetaReadSettingsChnls(pathSample)
  for (chnlName in names(chnlSettings)) {
    bwFallback <- chnlSettings[[chnlName]]$bwFallback
    expect_true(is.numeric(bwFallback) && bwFallback > 0)
    expect_equal(chnlSettings[[chnlName]]$biasUns, bwFallback)
  }
})
