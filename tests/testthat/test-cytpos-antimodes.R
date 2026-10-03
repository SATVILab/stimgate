test_that("cytokine antimodes preserve flat troughs and exclude boundaries", {
  expect_identical(.getCytPosTautStringAntimodes(NULL), numeric(0L))
  for (y in list(rep(1, 5), c(1, 1, 2, 2, 2), c(2, 2, 1, 1, 1))) {
    expect_identical(
      .getCytPosTautStringAntimodes(list(x = seq_along(y), y = y)),
      numeric(0L)
    )
  }
  expect_identical(
    .getCytPosTautStringAntimodes(list(x = 1:5, y = c(3, 1, 1, 1, 3))),
    3
  )
  expect_identical(
    .getCytPosTautStringAntimodes(list(x = 1:6, y = c(3, 1, NA, 1, 1, 3))),
    3.5
  )
})
