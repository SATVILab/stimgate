library(testthat)

test_that(".interp function works with interpolation (x_low != val)", {
  # Create test data where interpolation is needed
  x <- c(1, 3, 5, 7, 9)
  y <- c(10, 30, 50, 70, 90)
  
  # Test interpolation when val is between two x values
  # val = 4 should interpolate between x=3 (y=30) and x=5 (y=50)
  result <- stimgate:::.interp(val = 4, x = x, y = y)
  expected <- 30 + (4 - 3) * (50 - 30) / (5 - 3)  # Should be 40
  expect_equal(result, expected)
  
  # Another test case
  # val = 6 should interpolate between x=5 (y=50) and x=7 (y=70)
  result2 <- stimgate:::.interp(val = 6, x = x, y = y)
  expected2 <- 50 + (6 - 5) * (70 - 50) / (7 - 5)  # Should be 60
  expect_equal(result2, expected2)
  
  # Test with non-integer values
  x3 <- c(0.5, 1.5, 2.5, 3.5)
  y3 <- c(5, 15, 25, 35)
  # val = 2.0 should interpolate between x=1.5 (y=15) and x=2.5 (y=25)
  result3 <- stimgate:::.interp(val = 2.0, x = x3, y = y3)
  expected3 <- 15 + (2.0 - 1.5) * (25 - 15) / (2.5 - 1.5)  # Should be 20
  expect_equal(result3, expected3)
})

test_that(".interp function works with exact match (x_low == val)", {
  # Create test data where exact match occurs
  x <- c(1, 3, 5, 7, 9)
  y <- c(10, 30, 50, 70, 90)
  
  # Test exact match - should return corresponding y value
  result <- stimgate:::.interp(val = 5, x = x, y = y)
  expect_equal(result, 50)
  
  # Test another exact match
  result2 <- stimgate:::.interp(val = 1, x = x, y = y)
  expect_equal(result2, 10)
})