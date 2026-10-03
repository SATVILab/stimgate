test_that("unstimulated gates can be synthesized without cytokine thresholds", {
  gates <- tibble::tibble(
    chnl = "BC1", marker = "MarkerF1", batch = "batch_1", ind = "2",
    gate = 0.5
  )
  actual <- .fcsWriteGetGateTblAddUnsGetUnsImpl(
    gates, calc = min, indBatchList = list(c(1L, 2L))
  )
  expect_identical(actual$ind, "1")
  expect_identical(actual$gate, 0.5)
  expect_false("gateCyt" %in% names(actual))
})
