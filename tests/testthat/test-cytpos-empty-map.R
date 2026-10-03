test_that("all unstimulated samples leave cytokine gates unavailable", {
  gates <- tibble::tibble(
    batch = "b", ind = 1L, chnl = "A", marker = "A",
    gateName = "gate", gate = 1
  )
  testthat::local_mocked_bindings(
    .getLabs = function(...) c(A = "A"),
    .getCytPosGatesGateTblGet = function(...) gates
  )
  result <- .gateCytPos(
    chnlSettings = list(list(chnlCut = "A", popGate = "root", bwMin = 0.1)),
    .data = list(NULL), indBatchList = list(b = 1L),
    stage = "cytPos", pathProject = tempdir()
  )
  expect_identical(result, dplyr::mutate(gates, gateCyt = NA_real_))
})
