test_that("cluster adjustment preserves order without statistics", {
  gate_tbl <- tibble::tibble(
    gateName = "base_each", gateType = "base", gateCombn = "each",
    batch = "batch1", ind = c("2", "4"), gate = c(1, 2), gateUse = "gate"
  )
  testthat::local_mocked_bindings(
    .getStats = function(...) stop("statistics must not be calculated"),
    .getCpCluster = function(...) {
      tibble::tibble(ind = c("4", "2"), cpJoinTgOrig = c(3, 4))
    }
  )
  result <- .gateChnlGetAdjGatesAll(
    gate_tbl,
    .data = NULL, pathProject = tempdir(), stage = "init",
    indBatchList = list(batch1 = c(1, 2, 4)),
    chnlSettings = list(clusterGates = TRUE), calcCytPosGates = FALSE
  )$gateTbl
  expected <- dplyr::bind_rows(
    dplyr::select(gate_tbl, -gateUse),
    tibble::tibble(
      gateName = "base_eachClust", gateType = "base", gateCombn = "eachClust",
      batch = "batch1", ind = c("4", "2"), gate = c(3, 4)
    )
  )
  expect_identical(result, expected)
})

test_that("legacy detailed column mappings preserve existing values", {
  object <- tibble::tibble(
    threshold = NA_real_, locFinalThreshold = 2,
    locFinalThresholdOrigin = "direct", locFinalNCellStim = 10L,
    locFinalNCellUns = 12L, locFinalPropStim = 0.2,
    locFinalPropUns = 0.1, locFinalPropBs = 0.1
  )
  expected <- object
  expected$detailLevel <- "cluster_final"
  expected$thresholdOrigin <- "direct"
  expected$nCellStim <- 10L
  expected$nCellUns <- 12L
  expected$propStim <- 0.2
  expected$propUns <- 0.1
  expected$propBs <- 0.1
  expect_identical(
    .gateGetDetailedNormaliseObject(object, "locDetailClusterFinal"), expected
  )
  expect_identical(
    .gateGetDetailedNormaliseObject(object, "locDetailSample"), object
  )
})

test_that("empty detailed output is saved; missing directories are empty", {
  path_project <- tempfile("gate-empty-")
  expect_identical(.gateGetPop(path_project), character(0))
  expect_identical(.gateGetChnl(path_project, "root"), character(0))
  dir.create(path_project)
  withr::defer(unlink(path_project, recursive = TRUE))
  result <- getStimGatesDetailed(path_project, save = TRUE)
  expect_identical(result, tibble::tibble())
  expect_identical(
    readRDS(file.path(path_project, "gatesDetailed.rds")), result
  )
})

test_that("detailed channel fallback preserves existing channels", {
  path_project <- tempfile("gate-coalesce-")
  path_dir <- file.path(
    path_project, "intermediateData", "init", "BC1", "ind", "2"
  )
  dir.create(path_dir, recursive = TRUE)
  withr::defer(unlink(path_project, recursive = TRUE))
  saveRDS(
    tibble::tibble(chnl = c(NA_character_, "BC2"), threshold = c(1, 2)),
    file.path(path_dir, "locDetailSample.rds")
  )
  detail <- getStimGatesDetailed(path_project, save = TRUE)
  expect_identical(detail$chnl, c("BC1", "BC2"))
  expect_identical(detail$pop, rep(NA_character_, 2))
  expect_identical(
    readRDS(file.path(path_project, "gatesDetailed.rds")), detail
  )
})
