test_that("unknown gate markers receive a meaningful validation error", {
  path_project <- tempfile("gate-marker-")
  dir.create(file.path(path_project, "metaData"), recursive = TRUE)
  dir.create(file.path(path_project, "gates", "poproot"), recursive = TRUE)
  withr::defer(unlink(path_project, recursive = TRUE))
  saveRDS(c(BC1 = "IL2"), file.path(path_project, "metaData", "chnlLab.rds"))
  expect_error(
    getStimGates(path_project, marker = "unknown"), "Unknown marker: unknown"
  )
})
