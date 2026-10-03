test_that("detailed gates include current cluster thresholds", {
  pathProject <- tempfile("gate-cluster-")
  pathDir <- file.path(pathProject, "intermediateData", "init", "BC1", "ind", "all")
  dir.create(pathDir, recursive = TRUE)
  withr::defer(unlink(pathProject, recursive = TRUE))
  saveRDS(
    tibble::tibble(ind = "2", cpJoinTgOrig = 1.5, locClusterAction = "direct_retained"),
    file.path(pathDir, "locClusterQuantileTbl.rds")
  )
  detail <- getStimGatesDetailed(pathProject)
  expect_identical(detail$detailObject, "locClusterQuantileTbl")
  expect_identical(detail$detailLevel, "cluster_final")
  expect_identical(detail$threshold, 1.5)
  expect_identical(detail$locClusterAction, "direct_retained")
})
