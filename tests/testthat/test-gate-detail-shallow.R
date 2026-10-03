test_that("diagnostics without stage or channel directories remain readable", {
  path_project <- tempfile("gate-shallow-")
  path_int <- file.path(path_project, "intermediateData")
  dir.create(path_int, recursive = TRUE)
  withr::defer(unlink(path_project, recursive = TRUE))
  saveRDS(
    tibble::tibble(threshold = 1), file.path(path_int, "locDetailSample.rds")
  )
  detail <- getStimGatesDetailed(path_project)
  expect_identical(detail$threshold, 1)
  expect_identical(detail$detailPathStage, NA_character_)
  expect_identical(detail$chnl, NA_character_)
  expect_identical(detail$detailPathInd, NA_character_)
})
