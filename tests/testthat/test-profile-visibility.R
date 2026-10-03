test_that("statistics profiling preserves its debug-dependent visibility", {
  testthat::local_mocked_bindings(
    .profileOriginalGateStats = function(...) invisible(42)
  )
  path_project <- tempfile("profile-stats-visibility-")
  withr::defer({
    .profileStateReset()
    unlink(path_project, recursive = TRUE)
  })
  for (enabled in c("true", "false")) {
    withr::local_envvar(c(STIMGATE_DEBUG = enabled))
    .profileStateReset()
    result <- withVisible(.gateStats(
      .data = NULL, calcCytPosGates = FALSE, chnlSettings = NULL,
      indBatchList = NULL, pathProject = path_project
    ))
    expect_identical(result$value, 42)
    expect_identical(result$visible, identical(enabled, "true"))
  }
})
