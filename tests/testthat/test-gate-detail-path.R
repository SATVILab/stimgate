test_that("diagnostic paths treat project regex characters literally", {
  path_project <- tempfile("gate[one]+(two)")
  path_int <- file.path(path_project, "intermediateData")
  path_dir <- file.path(path_int, "init", "BC1", "ind", "2")
  dir.create(path_dir, recursive = TRUE)
  withr::defer(unlink(path_project, recursive = TRUE))
  saveRDS(
    tibble::tibble(threshold = 1), file.path(path_dir, "locDetailSample.rds")
  )
  detail <- getStimGatesDetailed(path_project)
  expect_identical(detail$detailPathStage, "init")
  expect_identical(detail$chnl, "BC1")
  expect_identical(detail$detailPathInd, "2")
})
