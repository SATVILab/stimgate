test_that("channel gating preserves the fixed temporary stimgate directory", {
  withr::local_envvar(STIMGATE_DEBUG = "false")
  pathProject <- file.path(tempdir(), "stimgate")
  dir.create(pathProject, showWarnings = FALSE)
  sentinel <- tempfile("preserve-", tmpdir = pathProject)
  saveRDS("existing project data", sentinel)
  withr::defer(unlink(sentinel))
  withr::defer(if (length(list.files(pathProject)) == 0L) unlink(pathProject))
  gateTbl <- tibble::tibble(ind = "2", gate = 1)
  testthat::local_mocked_bindings(
    .getLabs = function(...) NULL,
    .gateChnlPreAdjGatesGate = function(...) gateTbl,
    .gateChnlGetAdjGatesAll = function(...) list(gateTbl = gateTbl)
  )
  .gateChnl(
    .data = list(NULL), indBatchList = list(batch1 = c(1, 2)),
    chnlSettings = list(chnlCut = "BC1"), calcCytPosGates = FALSE,
    pathProject = pathProject, stage = "init"
  )
  expect_true(file.exists(sentinel))
  if (file.exists(sentinel)) {
    expect_identical(readRDS(sentinel), "existing project data")
  }
})
