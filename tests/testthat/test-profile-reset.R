test_that("profile initialization without reset preserves existing records", {
  withr::local_envvar(c(STIMGATE_DEBUG = "true"))
  .profileStateReset()
  path_project <- tempfile("profile-reset-")
  withr::defer({
    .profileStateReset()
    unlink(path_project, recursive = TRUE)
  })
  raw_dir <- file.path(path_project, "profile", "raw")
  dir.create(raw_dir, recursive = TRUE)
  path_record <- file.path(raw_dir, "existing.rds")
  saveRDS(list(existing = TRUE), path_record)

  expect_true(.profileInit(path_project, reset = FALSE))
  expect_true(file.exists(path_record))
})
