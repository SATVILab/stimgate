test_that("ambiguous duplicate numeric thresholds are rejected", {
  gates <- tibble::tibble(
    chnl = "BC1", marker = "MarkerF1", batch = "batch_1", ind = "2",
    gate = c(111, 11), gateCyt = c(1, 11)
  )
  expect_error(
    .fcsWriteGetGateTblAddUnsGetUnsImpl(gates, min, list(c(1L, 2L))),
    "Gates are not the same for all duplicates"
  )
  expect_error(
    .fcsWriteGetGateTblAddUnsGetUnsImpl(gates[2:1, ], min, list(c(1L, 2L))),
    "Gates are not the same for all duplicates"
  )
  gates$gate <- c(11, 11)
  gates$gateCyt <- c(1, 1)
  actual <- .fcsWriteGetGateTblAddUnsGetUnsImpl(gates, min, list(c(1L, 2L)))
  expect_identical(actual$gate, 11)
  expect_identical(actual$gateCyt, 1)
})
