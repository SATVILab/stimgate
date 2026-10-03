test_that("unknown gate markers receive a meaningful validation error", {
  pathProject <- tempfile("gate-marker-")
  dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
  dir.create(file.path(pathProject, "gates", "poproot"), recursive = TRUE)
  withr::defer(unlink(pathProject, recursive = TRUE))
  saveRDS(c(BC1 = "IL2"), file.path(pathProject, "metaData", "chnlLab.rds"))
  expect_error(getStimGates(pathProject, marker = "unknown"), "Unknown marker: unknown")
})
