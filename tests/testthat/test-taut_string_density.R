test_that(".tautStringPmden returns a list with y of length n-1", {
  set.seed(42)
  x <- sort(rnorm(200))
  result <- .tautStringPmden(x)
  expect_type(result, "list")
  expect_named(result, "y")
  expect_length(result$y, length(x) - 1L)
  expect_true(all(is.finite(result$y)))
  expect_true(all(result$y >= 0))
})

test_that(".tautStringPmden returns non-negative finite values for unimodal data", {
  # Fixture: set.seed(1), sort(rnorm(100, 5, 1)) → nmax=1, all finite, all >= 0
  set.seed(1)
  x <- sort(rnorm(100, mean = 5, sd = 1))
  result <- .tautStringPmden(x)
  expect_length(result$y, length(x) - 1L)
  expect_true(all(result$y >= 0))
  expect_true(all(is.finite(result$y)))
})

test_that(".tautStringPmden has higher density in mode than antimode for bimodal data", {
  # Fixture: set.seed(7), two clusters at 0 and 5 (sd=0.5), n=100
  # density[49]=0.163882, density[50]=0.003675 — sharp drop at the gap
  set.seed(7)
  x <- sort(c(rnorm(50, mean = 0, sd = 0.5), rnorm(50, mean = 5, sd = 0.5)))
  result <- .tautStringPmden(x)

  # The gap between the two clusters falls between obs 50 and 51 (0-based: 49 and 50)
  # density index 49 is in the left mode, index 50 is the gap
  expect_true(result$y[49] > result$y[50])
  # Both should be non-negative and finite
  expect_true(all(result$y >= 0))
  expect_true(all(is.finite(result$y)))
})

test_that(".tautStringPmden handles short input (n < 50) gracefully", {
  # n < 50: cpPmden Kuiper-bound table starts at n=50; below this return zeros
  expect_type(.tautStringPmden(numeric(0L)), "list")
  expect_length(.tautStringPmden(numeric(0L))$y, 0L)

  expect_type(.tautStringPmden(1.0), "list")
  expect_length(.tautStringPmden(1.0)$y, 0L)

  res2 <- .tautStringPmden(c(1.0, 2.0))
  expect_length(res2$y, 1L)
  expect_equal(res2$y, 0.0)

  # n=10: all zeros
  res10 <- .tautStringPmden(sort(rnorm(10)))
  expect_true(all(res10$y == 0))
  expect_length(res10$y, 9L)
})

test_that(".tautStringPmden bimodal density has a clear dip between clusters", {
  # Two-cluster fixture with deterministic, well-separated groups (n=120)
  x <- c(seq(-1, 1, length.out = 60), seq(4, 6, length.out = 60))
  result <- .tautStringPmden(x)
  # The gap falls between obs 60 and 61 (gap density at index 60)
  # Left mode (indices 25-35): higher density
  # Right mode (indices 80-90): higher density
  # Gap (index 60): very low density
  gap_density <- result$y[60]
  left_density <- mean(result$y[25:35])
  right_density <- mean(result$y[80:90])
  expect_true(left_density > gap_density)
  expect_true(right_density > gap_density)
  expect_true(all(is.finite(result$y)))
  expect_true(all(result$y >= 0))
})

test_that(".getCytPosTautStringDensity returns NULL for degenerate input", {
  expect_null(.getCytPosTautStringDensity(numeric(0L)))
  expect_null(.getCytPosTautStringDensity(rep(1.0, 3L)))
  expect_null(.getCytPosTautStringDensity(c(1, 2, NA, NA, NA)))
})

test_that(".getCytPosTautStringAntimodes returns empty for unimodal density", {
  # Unimodal piecewise-constant density: no internal antimodes
  density <- list(
    x = 1:10 + 0.5,
    y = c(1, 2, 3, 4, 5, 5, 4, 3, 2, 1)
  )
  antimodes <- .getCytPosTautStringAntimodes(density)
  expect_length(antimodes, 0L)
})

test_that(".getCytPosTautStringAntimodes detects a clear antimode", {
  # Two modes with a clear trough between them
  density <- list(
    x = seq(0.5, 9.5, by = 1),
    y = c(5, 8, 5, 2, 1, 2, 5, 8, 5, 2)
  )
  antimodes <- .getCytPosTautStringAntimodes(density)
  # Should find the trough near x=4
  expect_gte(length(antimodes), 1L)
  expect_true(any(antimodes > 3 & antimodes < 5))
})
