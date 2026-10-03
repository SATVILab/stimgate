test_that("channel gating preserves the fixed temporary stimgate directory", {
  withr::local_envvar(STIMGATE_DEBUG = "false")
  path_project <- file.path(tempdir(), "stimgate")
  created <- dir.create(path_project, showWarnings = FALSE)
  sentinel <- tempfile("preserve-", tmpdir = path_project)
  saveRDS("existing project data", sentinel)
  withr::defer({
    unlink(sentinel)
    if (created) unlink(path_project)
  })
  gate_tbl <- tibble::tibble(ind = "2", gate = 1)
  testthat::local_mocked_bindings(
    .getLabs = function(...) NULL,
    .gateChnlPreAdjGatesGate = function(...) gate_tbl,
    .gateChnlGetAdjGatesAll = function(...) list(gateTbl = gate_tbl)
  )
  .gateChnl(
    .data = list(NULL), indBatchList = list(batch1 = c(1, 2)),
    chnlSettings = list(chnlCut = "BC1"), calcCytPosGates = FALSE,
    pathProject = path_project, stage = "init"
  )
  expect_true(file.exists(sentinel))
  if (file.exists(sentinel)) {
    expect_identical(readRDS(sentinel), "existing project data")
  }
})
