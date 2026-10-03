test_that("cluster bandwidth respects custom normalization settings", {
  withr::local_preserve_seed()
  x <- exp(seq(-2, 2, length.out = 100))
  settings <- list(
    bwMtd = "nrd0Norm", bwAdj = 1, normMtd = "boxcox",
    normExtraFrac = 0, normDensityN = 64L, normLambda = c(0, 1)
  )
  set.seed(1)
  expected <- do.call(.bwCalcOne, c(list(x = x), settings))
  set.seed(1)
  actual <- .getCpClusterLocBwOne(x, settings)
  expect_equal(as.numeric(actual), as.numeric(expected))
})
