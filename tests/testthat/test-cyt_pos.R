library(testthat)

test_that("cytPos gates actually happen", {
  testData <- getTestData(
    scenario = "cytPos",
    dir_cache = testthat::test_path("cache", "test_data", "default"),
    clear = TRUE,
    n_ind = 1
  )
  gs <- flowWorkspace::load_gs(testData$path_gs)
  pathProject <- file.path(dirname(testData$path_gs), "stimgate")
  invisible(stimgate_gate(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = testData$batchList,
    marker = testData$marker
  ))
  expect_true(file.exists(file.path(pathProject, "gateStats.rds")))
})