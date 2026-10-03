test_that("failed debug file creation does not initialize debug state", {
  withr::local_envvar(c(STIMGATE_DEBUG = "true"))
  .debugStateReset()
  path_project <- tempfile("debug-create-failure-")
  withr::defer({
    .debugStateReset()
    unlink(path_project, recursive = TRUE)
  })
  writeLines("blocks directory creation", path_project)

  expect_false(.debugInit(path_project))
  expect_false(.debugState$initialized)
  expect_null(.debugState$file)
})
