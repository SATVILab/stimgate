test_that("axis_limits works", {
  p <- readRDS(testthat::test_path("p_axis_limits.rds"))

  # tests
  # -----------------

  # one element of length 1, no name
  p_adj <- axis_limits(
    p = p,
    limits_expand = list(-1e4)
  )
  expect_identical(
    length(p_adj$layers),
    2L
  )
  expect_identical(
    p_adj$layers[[2]]$data,
    data.frame(
      x = c(-1e4, -1e4),
      y = c(-1e4, -1e4)
    )
  )

  # two elements, no name
  expect_error(
    axis_limits(
      p = p,
      limits_expand = list(1e4, -5e2)
    )
  )

  # one element, no name
  p_adj <- axis_limits(
    p = p,
    limits_expand = list(c(1e4, -5e2))
  )
  expect_identical(
    p_adj$layers[[2]]$data,
    data.frame(
      x = c(-5e2, 1e4),
      y = c(-5e2, 1e4)
    )
  )

  # one element, one name
  p_adj <- axis_limits(
    p = p,
    limits_expand = list(x = c(1e4, -5e2))
  )
  expect_identical(
    p_adj$layers[[2]]$data,
    data.frame(
      x = c(-5e2, 1e4)
    )
  )
  p_adj <- axis_limits(
    p = p,
    limits_expand = list(y = c(1e4, -5e2))
  )
  expect_identical(
    p_adj$layers[[2]]$data,
    data.frame(
      y = c(-5e2, 1e4)
    )
  )

  # two elements, both named
  p_adj <- axis_limits(
    p = p,
    limits_expand = list(
      y = c(1e4, -5e2),
      x = c(-1e4, 2e4)
    )
  )
  expect_identical(
    p_adj$layers[[2]]$data,
    data.frame(
      y = c(-5e2, 1e4),
      x = c(-1e4, 2e4)
    )
  )

  # axis range equal
  # --------------------

  # just axis range equal
  p_adj <- axis_limits(
    p = p,
    limits_equal = TRUE
  )

  expect_identical(
    p_adj$layers[[2]]$data[, 1],
    p_adj$layers[[2]]$data[, 2]
  )

  # with limits_expand
  # just axis range equal
  p_adj <- axis_limits(
    p = p,
    limits_equal = TRUE,
    limits_expand = list(
      y = c(1000, 200),
      x = c(-1e4, 500)
    )
  )
  expect_identical(
    p_adj$layers[[2]]$data[1, ] |>
      as.numeric(),
    c(-1e4, -1e4)
  )
  expect_identical(
    p_adj$layers[[2]]$data[2, ] |>
      as.numeric() |>
      round(),
    c(9222, 9222)
  )

  # just y-axis
  p_adj <- axis_limits(
    p = p,
    limits_equal = TRUE,
    limits_expand = list(y = c(1e4, 200))
  )

  expect_identical(
    p_adj$layers[[2]]$data[1, ] |>
      as.numeric(),
    c(1, 1)
  )
  expect_identical(
    p_adj$layers[[2]]$data[2, ] |>
      as.numeric() |>
      round(),
    c(1e4, 9222)
  )

  # just x-axis
  p_adj <- axis_limits(
    p = p,
    limits_equal = TRUE,
    limits_expand = list(x = c(1e4, 200))
  )

  expect_identical(
    p_adj$layers[[2]]$data[1, ] |>
      as.numeric(),
    c(1, 1)
  )
  expect_identical(
    p_adj$layers[[2]]$data[2, ] |>
      as.numeric() |>
      round(),
    c(1e4, 9222)
  )
})
