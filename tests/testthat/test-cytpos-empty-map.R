test_that("all unstimulated samples leave cytokine gates unavailable", {
  gates <- tibble::tibble(
    batch = "b", ind = 1L, chnl = "A", marker = "A",
    gateName = "gate", gate = 1
  )
  result <- .getCytPosGatesGateName(
    gateTblGn = gates, .data = NULL, indBatchList = list(b = 1L),
    chnlVec = "A", chnlLabVec = c(A = "A"), popGate = "root",
    bwMin = 0.1, calcCytPos = TRUE, stage = "cytPos",
    pathProject = tempdir()
  )
  expect_identical(result, dplyr::mutate(gates, gateCyt = NA_real_))
})
