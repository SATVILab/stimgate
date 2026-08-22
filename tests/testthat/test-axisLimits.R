test_that("axisLimits works", {
  p <- ggplot2::ggplot(
    data.frame(x = c(1, 9222), y = c(1, 9222)),
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_point()

  # tests
  # -----------------

  # one element of length 1, no name
  pAdj <- axisLimits(
    p = p,
    limitsExpand = list(-1e4)
  )
  expect_equal(
    length(pAdj$layers),
    2L
  )
  dat <- ggplot2::layer_data(pAdj, 2)
  expect_equal(
    dat[c("x", "y")],
    data.frame(
      x = c(-1e4, -1e4),
      y = c(-1e4, -1e4)
    )
  )

  # two elements, no name
  expect_error(
    axisLimits(
      p = p,
      limitsExpand = list(1e4, -5e2)
    )
  )

  # one element, no name
  pAdj <- axisLimits(
    p = p,
    limitsExpand = list(c(1e4, -5e2))
  )
  dat <- ggplot2::layer_data(pAdj, 2)
  expect_equal(
    dat[c("x", "y")],
    data.frame(
      x = c(-5e2, 1e4),
      y = c(-5e2, 1e4)
    )
  )

  # one element, one name
  pAdj <- axisLimits(
    p = p,
    limitsExpand = list(x = c(1e4, -5e2))
  )
  dat <- ggplot2::layer_data(pAdj, 2)
  expect_equal(
    dat["x"],
    data.frame(
      x = c(-5e2, 1e4)
    )
  )
  pAdj <- axisLimits(
    p = p,
    limitsExpand = list(y = c(1e4, -5e2))
  )
  dat <- ggplot2::layer_data(pAdj, 2)
  expect_equal(
    dat["y"],
    data.frame(
      y = c(-5e2, 1e4)
    )
  )

  # two elements, both named
  pAdj <- axisLimits(
    p = p,
    limitsExpand = list(
      y = c(1e4, -5e2),
      x = c(-1e4, 2e4)
    )
  )
  dat <- ggplot2::layer_data(pAdj, 2)
  expect_equal(
    dat[c("y", "x")],
    data.frame(
      y = c(-5e2, 1e4),
      x = c(-1e4, 2e4)
    )
  )

  # axis range equal
  # --------------------

  # just axis range equal
  pAdj <- axisLimits(
    p = p,
    limitsEqual = TRUE
  )
  dat <- ggplot2::layer_data(pAdj, 2)
  expect_equal(
    dat$x,
    dat$y
  )

  # with limitsExpand
  # just axis range equal
  pAdj <- axisLimits(
    p = p,
    limitsEqual = TRUE,
    limitsExpand = list(
      y = c(1000, 200),
      x = c(-1e4, 500)
    )
  )
  dat <- ggplot2::layer_data(pAdj, 2)
  expect_equal(
    as.numeric(dat[1, c("x", "y")]),
    c(-1e4, -1e4)
  )
  expect_equal(
    round(as.numeric(dat[2, c("x", "y")])),
    c(9222, 9222)
  )

  # just y-axis
  pAdj <- axisLimits(
    p = p,
    limitsEqual = TRUE,
    limitsExpand = list(y = c(1e4, 200))
  )
  dat <- ggplot2::layer_data(pAdj, 2)
  expect_equal(
    as.numeric(dat[1, c("y", "x")]),
    c(1, 1)
  )
  expect_equal(
    round(as.numeric(dat[2, c("y", "x")])),
    c(1e4, 9222)
  )

  # just x-axis
  pAdj <- axisLimits(
    p = p,
    limitsEqual = TRUE,
    limitsExpand = list(x = c(1e4, 200))
  )
  dat <- ggplot2::layer_data(pAdj, 2)
  expect_equal(
    as.numeric(dat[1, c("x", "y")]),
    c(1, 1)
  )
  expect_equal(
    round(as.numeric(dat[2, c("x", "y")])),
    c(1e4, 9222)
  )
})
