test_that("stimgateGateRuns", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "stimgate")
  invisible(gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  ))
  expect_true(file.exists(file.path(pathProject, "gateStats.rds")))
})
