test_that("final clustered gates preserve donor, transferred and failed provenance", {
  gates <- tibble::tibble(
    gateUse = "gate", gateName = "loc_min", gateType = "loc_",
    gateCombn = "min", batch = "a", ind = c("2", "3", "4"),
    gate = c(1, 100, 100), locGenerated = c(TRUE, FALSE, FALSE),
    locGeneratedDirect = c(TRUE, FALSE, FALSE),
    locSource = c("direct", "high_value", "high_value"),
    locReason = c(NA_character_, "no_threshold", "no_threshold")
  )
  clustered <- tibble::tibble(
    ind = gates$ind, cpJoinTgOrig = c(1, 2, 100),
    locGenerated = c(TRUE, TRUE, FALSE),
    locGeneratedDirect = c(TRUE, FALSE, FALSE),
    locSource = c("direct", "cluster_q60", "high_value"),
    locReason = c(NA_character_, "replaced_by_cluster_direct_threshold_q60", "no_threshold")
  )
  testthat::local_mocked_bindings(
    .getCpCluster = function(...) clustered,
    .package = "stimgate"
  )
  result <- stimgate:::.gateChnlGetAdjGatesAll(
    gateTbl = gates, .data = NULL, pathProject = tempdir(), stage = "init",
    indBatchList = list(a = 1:4), chnlSettings = list(clusterGates = TRUE),
    calcCytPosGates = TRUE
  )$gateTbl |>
    dplyr::filter(.data$gateName == "loc_minClust")
  expect_identical(result$gate, clustered$cpJoinTgOrig)
  expect_identical(result$locGenerated, clustered$locGenerated)
  expect_identical(result$locGeneratedDirect, clustered$locGeneratedDirect)
  expect_identical(result$locSource, clustered$locSource)
  expect_identical(result$locReason, clustered$locReason)
})

test_that("cytokine refinement retains provenance through final gate persistence", {
  project <- tempfile("cluster-provenance-")
  withr::defer(unlink(project, recursive = TRUE))
  initial <- tibble::tibble(
    gateName = "loc_minClust", batch = "a", ind = c("2", "3"),
    gate = c(2, 100), locGenerated = c(TRUE, FALSE),
    locGeneratedDirect = FALSE, locSource = c("cluster_q60", "high_value"),
    locReason = c("transferred", "no_threshold")
  )
  path <- stimgate:::.gatesGetPathAll(project, "root", "X", TRUE)
  dir.create(dirname(path), recursive = TRUE)
  saveRDS(initial, path)
  testthat::local_mocked_bindings(
    .getLabs = function(...) c(X = "IFNg"),
    .getCytPosGatesInd = function(ind, ...) {
      if (ind == 1L) return(NULL)
      tibble::tibble(batch = "a", ind = as.character(ind), chnl = "X",
                     marker = "IFNg", gateCyt = 1)
    },
    .package = "stimgate"
  )
  refined <- stimgate:::.gateCytPos(
    chnlSettings = list(list(chnlCut = "X", popGate = "root", bwMin = 1)),
    indBatchList = list(a = 1:3), .data = list(NULL),
    calcCytPos = TRUE, stage = "cytPos", pathProject = project
  )
  stimgate:::.getStatsGateTblSave(
    gateTbl = refined, pathProject = project, popGate = "root",
    chnlLab = c(X = "IFNg"), chnl = "X", save = TRUE
  )
  final <- readRDS(stimgate:::.gatesGetPathAll(project, "root", "X", FALSE))
  expect_identical(final$locGenerated, initial$locGenerated)
  expect_identical(final$locGeneratedDirect, initial$locGeneratedDirect)
  expect_identical(final$locSource, initial$locSource)
  expect_identical(final$locReason, initial$locReason)
  expect_equal(final$gateCyt, c(1, 1))
})
