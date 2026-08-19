test_that("cytPos gates actually happen", {
  testData <- getTestData(
    scenario = "cytPos",
    dirCache = testthat::test_path("cache", "test_data", "cytPos"),
    clear = TRUE,
    nCell = 300L,
    nInd = 1
  )
  gs <- flowWorkspace::load_gs(testData$pathGs)
  pathProject <- file.path(dirname(testData$pathGs), "stimgate")
  gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = testData$batchList,
    marker = testData$marker
  )

  # Verify the gate table has expected structure and populated values
  gateTbl <- getStimGates(pathProject)
  expect_true(is.data.frame(gateTbl))
  expect_true(nrow(gateTbl) > 0)
  expect_true(all(c("chnl", "marker", "batch", "ind", "gate", "gateCyt") %in% colnames(gateTbl)))
  expect_true(is.numeric(gateTbl$gate))
  expect_true(all(is.finite(gateTbl$gate)))
  # cytPos scenario uses low separation; both markers should have finite gates
  expect_equal(nrow(gateTbl), length(testData$chnl))
})