test_that("failed debug file creation does not initialize debug state", {
  withr::local_envvar(c(STIMGATE_DEBUG = "true"))
  .debugStateReset()
  pathProject <- tempfile("debug-create-failure-")
  withr::defer({
    .debugStateReset()
    unlink(pathProject, recursive = TRUE)
  })
  writeLines("blocks directory creation", pathProject)

  expect_false(.debugInit(pathProject))
  expect_false(.debugState$initialized)
  expect_null(.debugState$file)
})
