test_that("detailed gates include current cluster thresholds", {
  path_project <- tempfile("gate-cluster-")
  path_dir <- file.path(
    path_project, "intermediateData", "init", "BC1", "ind", "all"
  )
  dir.create(path_dir, recursive = TRUE)
  withr::defer(unlink(path_project, recursive = TRUE))
  saveRDS(
    tibble::tibble(
      ind = "2", cpJoinTgOrig = 1.5, locClusterAction = "direct_retained"
    ),
    file.path(path_dir, "locClusterQuantileTbl.rds")
  )
  detail <- getStimGatesDetailed(path_project)
  expect_identical(detail$detailObject, "locClusterQuantileTbl")
  expect_identical(detail$detailLevel, "cluster_final")
  expect_identical(detail$threshold, 1.5)
  expect_identical(detail$locClusterAction, "direct_retained")
})
