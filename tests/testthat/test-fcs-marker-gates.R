test_that("control gates can be synthesized without marker labels", {
  example <- getExampleData()
  withr::defer(unlink(dirname(example$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(example$pathGs)
  gates <- tibble::tibble(
    chnl = example$chnl[[1]], batch = "batch_1", ind = "2",
    gate = 0.5, gateCyt = 0.25
  )
  actual <- .fcsWriteGetGateTbl(
    gates,
    chnl = example$chnl[[1]], pop = "root", .data = gs,
    indBatchList = list(c(1L, 2L)), gateUnsMethod = "min",
    pathProject = tempdir()
  )
  expect_identical(actual$ind, c("1", "2"))
  expect_identical(unname(actual$marker), rep(example$marker[[1]], 2L))
  expect_identical(actual$gate, c(0.5, 0.5))
})
