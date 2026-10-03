test_that("cluster grids reject sizes that cannot be represented as integers", {
  expect_error(.getCpClusterControlUpdate(list(nGrid = Inf)))
  expect_error(.getCpClusterControlUpdate(list(nGrid = .Machine$integer.max + 1)))
  control <- .getCpClusterControlUpdate(list(nGrid = 32))
  expect_identical(control$nGrid, 32L)
  expect_length(.getCpClusterLocDensityGrid(0, 1, control$nGrid), 32L)
})
