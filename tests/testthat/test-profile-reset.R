test_that("profile initialization without reset preserves existing records", {
  withr::local_envvar(c(STIMGATE_DEBUG = "true"))
  .profileStateReset()
  pathProject <- tempfile("profile-reset-")
  withr::defer({
    .profileStateReset()
    unlink(pathProject, recursive = TRUE)
  })
  rawDir <- file.path(pathProject, "profile", "raw")
  dir.create(rawDir, recursive = TRUE)
  pathRecord <- file.path(rawDir, "existing.rds")
  saveRDS(list(existing = TRUE), pathRecord)

  expect_true(.profileInit(pathProject, reset = FALSE))
  expect_true(file.exists(pathRecord))
})
