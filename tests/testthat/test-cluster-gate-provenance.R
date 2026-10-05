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
