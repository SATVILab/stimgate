library(testthat)

test_that("pkgdown reference includes exported gating stats function", {
  pkgdown_config <- readLines(testthat::test_path("../../_pkgdown.yml"))
  expect_true(any(grepl("getStimStats", pkgdown_config, fixed = TRUE)))
})
