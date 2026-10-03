test_that("fresh worker attachment preserves active sample context", {
  withr::local_envvar(c(STIMGATE_DEBUG = "true"))
  .profileStateReset()
  pathProject <- tempfile("profile-attach-")
  withr::defer({
    .profileStateReset()
    unlink(pathProject, recursive = TRUE)
  })

  timer <- .profileWithContext(
    .profileStart("sample", "initial_gating", pathProject = pathProject),
    marker = "IL2", channel = "IL2-A", batch = "donor",
    sample = "104", stage = "init"
  )
  expect_identical(timer$marker, "IL2")
  expect_identical(timer$channel, "IL2-A")
  expect_identical(timer$batch, "donor")
  expect_identical(timer$sample, "104")
  expect_identical(timer$stage, "init")
  expect_identical(.profileState$context, .profileContextDefault())
})
