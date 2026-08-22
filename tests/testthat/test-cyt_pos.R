test_that("cytPos gates are refined when the canonical example data exercises the branch", {
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "stimgate_cytpos_example")
  on.exit(unlink(pathProject, recursive = TRUE), add = TRUE)

  gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  )

  gateTbl <- getStimGates(pathProject)
  expect_true(is.data.frame(gateTbl))
  expect_true(nrow(gateTbl) > 0)
  expect_true(all(c("chnl", "marker", "batch", "ind", "gate", "gateCyt") %in% colnames(gateTbl)))
  expect_true(all(is.finite(gateTbl$gate)))
  expect_true(all(is.finite(gateTbl$gateCyt)))
  expect_true(all(gateTbl$gateCyt <= gateTbl$gate, na.rm = TRUE))
  expect_true(any(gateTbl$gateCyt < gateTbl$gate, na.rm = TRUE))
})
