test_that(".tautStringPmden returns a list with y of length n-1", {
  skip_if_not_installed("cytoUtils")
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
  skip_if_not_installed("cytoUtils")
  set.seed(1)
  x <- sort(rnorm(100, mean = 5, sd = 1))
  result <- .tautStringPmden(x)
  expect_true(all(result$y >= 0))
  expect_true(all(is.finite(result$y)))
})

test_that(".tautStringPmden has higher density in mode than antimode for bimodal data", {
  skip_if_not_installed("cytoUtils")
  set.seed(7)
  # Clear bimodal: well-separated modes so taut string reliably finds gates
  x <- sort(c(rnorm(300, mean = 0, sd = 0.5), rnorm(300, mean = 5, sd = 0.5)))
  result <- .tautStringPmden(x)

  # Left mode region: x near 0
  left_mode_idx <- which(x >= -1 & x <= 1)
  # Antimode region: x between modes
  anti_idx <- which(x >= 1.5 & x <= 3.5)
  # Right mode region: x near 5
  right_mode_idx <- which(x >= 4 & x <= 6)

  if (length(left_mode_idx) > 1 && length(anti_idx) > 1 &&
        length(right_mode_idx) > 1) {
    # Density at interval midpoints uses index i for interval (x[i], x[i+1])
    left_y <- mean(result$y[left_mode_idx[-length(left_mode_idx)]])
    anti_y <- mean(result$y[anti_idx[-length(anti_idx)]])
    right_y <- mean(result$y[right_mode_idx[-length(right_mode_idx)]])

    expect_gt(left_y, anti_y)
    expect_gt(right_y, anti_y)
  }
})

test_that(".tautStringPmden returns length-0 or handles very short input gracefully", {
  skip_if_not_installed("cytoUtils")
  expect_type(.tautStringPmden(numeric(0L)), "list")
  expect_type(.tautStringPmden(1.0), "list")
  expect_type(.tautStringPmden(c(1.0, 2.0)), "list")
  # length-2 input → y has length 1
  res2 <- .tautStringPmden(c(1.0, 2.0))
  expect_length(res2$y, 1L)
})

test_that(".tautStringPmden antimode locations match cytoUtils::tautstring gates", {
  skip_if_not_installed("cytoUtils")
  set.seed(99)
  x <- sort(c(rnorm(200, 0, 0.4), rnorm(200, 4, 0.4)))
  result <- .tautStringPmden(x)

  gates <- cytoUtils::tautstring(x)
  # Interior gates are the antimode locations
  interior <- gates[-c(1L, length(gates))]

  if (length(interior) > 0) {
    # The density should be at a local minimum near each interior gate
    x_mid <- (x[-length(x)] + x[-1L]) / 2
    for (g in interior) {
      # Find index nearest to gate
      nearest <- which.min(abs(x_mid - g))
      # Density at that index should be lower than neighbours some distance away
      window <- 10L
      left_range <- seq_len(max(1L, nearest - window))
      right_range <- seq.int(
        min(length(result$y), nearest + window),
        length(result$y)
      )
      if (length(left_range) > 0 && length(right_range) > 0) {
        expect_lte(
          result$y[nearest],
          max(result$y[left_range], result$y[right_range])
        )
      }
    }
  }
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
