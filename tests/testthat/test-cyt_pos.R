library(testthat)

test_that("cyt_pos gates actually happen", {
  testData <- getTestData(
    scenario = "cyt_pos",
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
    batchList = testData$batch_list,
    marker = testData$marker
  ))
  expect_true(file.exists(file.path(pathProject, "gate_stats.rds")))
})