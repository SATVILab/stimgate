test_that("complete gate tables support a shared unstimulated sample", {
  example <- getExampleData()
  withr::defer(unlink(dirname(example$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(example$pathGs)
  gates <- tibble::tibble(
    chnl = example$chnl[[1]], marker = example$marker[[1]],
    batch = c("batch_1", "batch_1", "batch_2", "batch_2"),
    ind = c("1", "2", "1", "3"), gate = 0.5, gateCyt = 0.25
  )
  actual <- .fcsWriteGetGateTbl(
    gates, chnl = example$chnl[[1]], pop = "root", .data = gs,
    indBatchList = list(c(1L, 2L), c(1L, 3L)), gateUnsMethod = "min",
    gateTypeCytPos = "base", pathProject = tempdir()
  )
  actual$marker <- unname(actual$marker)
  expect_identical(actual, gates)
})
