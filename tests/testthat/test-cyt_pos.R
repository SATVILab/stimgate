library(testthat)

test_that("cytPos gates actually happen", {
  testData <- getTestData(
    scenario = "cytPos",
    dirCache = testthat::test_path("cache", "test_data", "default"),
    clear = TRUE,
    nInd = 1
  )
  gs <- flowWorkspace::load_gs(testData$pathGs)
  pathProject <- file.path(dirname(testData$pathGs), "stimgate")
  invisible(gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = testData$batchList,
    marker = testData$marker
  ))
  expect_true(file.exists(file.path(pathProject, "gateStats.rds")))
})