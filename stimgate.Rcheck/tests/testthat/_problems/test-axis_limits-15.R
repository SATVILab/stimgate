# Extracted from test-axis_limits.R:15

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "stimgate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
p <- readRDS(testthat::test_path("p_axis_limits.rds"))
p_adj <- axis_limits(
    p = p,
    limits_expand = list(-1e4)
  )
expect_identical(
    length(p_adj$layers),
    2L
  )
