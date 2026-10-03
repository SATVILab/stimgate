test_that("getStimGatesDetailed populates and filters by population", {
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- tempfile("gate-pop-")
  withr::defer(unlink(pathProject, recursive = TRUE))
  withr::local_envvar(STIMGATE_INTERMEDIATE = "all")

  gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = list(batch1 = c(1, 2, 4)),
    marker = exampleData$marker[[1]]
  )

  detailAll <- getStimGatesDetailed(pathProject)
  expect_gt(nrow(detailAll), 0L)
  expect_false(any(is.na(detailAll$pop)))
  expect_true(all(detailAll$pop == "root"))

  detailRoot <- getStimGatesDetailed(pathProject, pop = "root")
  expect_identical(detailRoot, detailAll)

  detailNone <- getStimGatesDetailed(pathProject, pop = "otherPop")
  expect_equal(nrow(detailNone), 0L)
  expect_identical(names(detailNone), names(detailAll))
})
