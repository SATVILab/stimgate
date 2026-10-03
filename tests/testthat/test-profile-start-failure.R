test_that("timer-start errors leave the gating expression and RNG intact", {
  withr::local_envvar(c(STIMGATE_DEBUG = "true"))
  testthat::local_mocked_bindings(
    .profileStart = function(...) stop("timer start failure")
  )
  withr::local_seed(1)
  expected <- runif(1)
  set.seed(1)
  expect_identical(
    .profileTime(runif(1), level = "major", major = "initial_gating"),
    expected
  )
  expect_false(withVisible(
    .profileTime(invisible(42), level = "major", major = "initial_gating")
  )$visible)
})
