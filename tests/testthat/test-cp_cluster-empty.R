test_that("clustering without stimulated samples retains its schema", {
  gate_tbl <- tibble::tibble(ind = "1", gate = 2, locSource = "unstim_summary")
  testthat::local_mocked_bindings(
    .getCpClusterLocExLookup = function(...) list(),
    .intSave = function(...) NULL
  )
  out <- .getCpCluster(
    .data = NULL, gateTbl = gate_tbl,
    chnlSettings = list(chnlCut = "expr"), stage = "init",
    pathProject = tempdir(), filterOtherCytPos = FALSE,
    calcCytPosGates = FALSE, indBatchList = list()
  )
  expect_equal(nrow(out), 0L)
  expect_identical(
    names(out),
    names(.getCpClusterLocSkipOut(
      .getCpClusterLocGateTblPrepare(gate_tbl), "test"
    ))
  )
  expect_identical(
    dplyr::select(out, ind, cpJoinTgOrig),
    tibble::tibble(ind = character(), cpJoinTgOrig = double())
  )
})
