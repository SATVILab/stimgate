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
    gates,
    chnl = example$chnl[[1]], pop = "root", .data = gs,
    indBatchList = list(c(1L, 2L), c(1L, 3L)), gateUnsMethod = "min",
    pathProject = tempdir()
  )
  actual$marker <- unname(actual$marker)
  expect_identical(actual, gates)
})

test_that("unstim gates match batches regardless of index digit width", {
  gates <- tibble::tibble(
    chnl = "c1", marker = "m1",
    batch = c("b1", "b1", "b2", "b2"),
    ind = c("9", "10", "12", "13"), gate = c(1, 2, 3, 4)
  )
  # stim 13 has no gate row; its batch must still be matched
  gates <- gates[gates$ind != "13", ]
  actual <- .fcsWriteGetGateTblAddUnsGetUnsImpl(
    gates,
    calc = min,
    indBatchList = list(b1 = c(8, 9, 10), b2 = c(11, 12, 13))
  )
  expect_identical(actual$ind, c("8", "11"))
  expect_identical(actual$gate, c(1, 3))
})
