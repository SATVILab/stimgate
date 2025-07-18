test_that("axis_limits works with basic functionality", {
  # Setup test data
  data("cars", package = "datasets")
  library(ggplot2)
  p <- ggplot(cars, aes(speed, dist)) + geom_point()
  
  # Test no modifications (should return original plot)
  result <- axis_limits(p)
  expect_identical(result, p)
  
  # Test with NULL expand and FALSE equal (should return original plot)
  result <- axis_limits(p, limits_expand = NULL, limits_equal = FALSE)
  expect_identical(result, p)
})

test_that("axis_limits works with limits_equal = TRUE", {
  # Setup test data
  data("cars", package = "datasets")
  library(ggplot2)
  p <- ggplot(cars, aes(speed, dist)) + geom_point()
  
  # Test limits_equal = TRUE
  result <- axis_limits(p, limits_equal = TRUE)
  expect_s3_class(result, c("gg", "ggplot"))
  expect_true(length(result$layers) >= length(p$layers))
  
  # The function should add expand_limits layer
  layer_types <- sapply(result$layers, function(x) class(x$geom)[1])
  expect_true(any(grepl("Geom", layer_types)))
})

test_that("axis_limits works with limits_expand unnamed list", {
  # Setup test data
  data("cars", package = "datasets")
  library(ggplot2)
  p <- ggplot(cars, aes(speed, dist)) + geom_point()
  
  # Test with unnamed list (applies to both axes)
  result <- axis_limits(p, limits_expand = list(200))
  expect_s3_class(result, c("gg", "ggplot"))
  
  # Test with multiple values in unnamed list
  result <- axis_limits(p, limits_expand = list(c(-50, 200)))
  expect_s3_class(result, c("gg", "ggplot"))
})

test_that("axis_limits works with limits_expand named list", {
  # Setup test data
  data("cars", package = "datasets")
  library(ggplot2)
  p <- ggplot(cars, aes(speed, dist)) + geom_point()
  
  # Test with x only
  result <- axis_limits(p, limits_expand = list(x = 75))
  expect_s3_class(result, c("gg", "ggplot"))
  
  # Test with y only
  result <- axis_limits(p, limits_expand = list(y = 200))
  expect_s3_class(result, c("gg", "ggplot"))
  
  # Test with both x and y
  result <- axis_limits(p, limits_expand = list(x = c(-10, 75), y = c(-50, 200)))
  expect_s3_class(result, c("gg", "ggplot"))
})

test_that("axis_limits works with combined limits_equal and limits_expand", {
  # Setup test data
  data("cars", package = "datasets")
  library(ggplot2)
  p <- ggplot(cars, aes(speed, dist)) + geom_point()
  
  # Test combined functionality
  result <- axis_limits(
    p, 
    limits_expand = list(x = c(-10, 75), y = c(-50, 200)),
    limits_equal = TRUE
  )
  expect_s3_class(result, c("gg", "ggplot"))
  
  # Test with unnamed expand and equal limits
  result <- axis_limits(
    p,
    limits_expand = list(300),
    limits_equal = TRUE
  )
  expect_s3_class(result, c("gg", "ggplot"))
})

test_that("axis_limits validates input correctly", {
  # Setup test data
  data("cars", package = "datasets")
  library(ggplot2)
  p <- ggplot(cars, aes(speed, dist)) + geom_point()
  
  # Test invalid ggplot object (only throws error when function tries to process it)
  expect_error(
    axis_limits("not_a_plot", limits_equal = TRUE),
    "p must be of class c\\('gg', 'ggplot'\\)"
  )
  
  # Test invalid limits_equal
  expect_error(
    axis_limits(p, limits_equal = "not_logical"),
    "limits_equal must be logical"
  )
  # Test invalid limits_equal - vector of logicals
  expect_error(
    axis_limits(p, limits_equal = c(TRUE, FALSE)),
    "'length = 2' in coercion to 'logical\\(1\\)'"
  )
  
  # Test invalid limits_expand - not a list
  expect_error(
    axis_limits(p, limits_expand = "not_a_list"),
    "limits_expand must be a list \\(if not NULL\\)"
  )
  
  # Test invalid limits_expand - unnamed list of length 2
  expect_error(
    axis_limits(p, limits_expand = list(1, 2)),
    "limits_expand must be named if of length 2"
  )
  
  # Test invalid limits_expand - too many elements
  expect_error(
    axis_limits(p, limits_expand = list(x = 1, y = 2, z = 3)),
    "limits_expand must have length 1 or 2 \\(if not NULL\\)"
  )
  
  # Test invalid limits_expand - wrong names
  expect_error(
    axis_limits(p, limits_expand = list(a = 1)),
    "limits_expand must have names of 'x' and/or 'y' \\(if named\\)"
  )
  
  # Test invalid limits_expand - non-numeric values
  expect_error(
    axis_limits(p, limits_expand = list(x = "not_numeric")),
    "input to limits_expand must be numeric \\(if limits_expand not NULL\\)"
  )
})

test_that("axis_limits handles edge cases", {
  library(ggplot2)
  
  # Test with single point data
  single_data <- data.frame(x = 5, y = 10)
  p_single <- ggplot(single_data, aes(x, y)) + geom_point()
  
  result <- axis_limits(p_single, limits_equal = TRUE)
  expect_s3_class(result, c("gg", "ggplot"))
  
  result <- axis_limits(p_single, limits_expand = list(20))
  expect_s3_class(result, c("gg", "ggplot"))
  
  # Test with identical x and y ranges
  identical_data <- data.frame(x = 1:5, y = 1:5)
  p_identical <- ggplot(identical_data, aes(x, y)) + geom_point()
  
  result <- axis_limits(p_identical, limits_equal = TRUE)
  expect_s3_class(result, c("gg", "ggplot"))
  
  # Test with negative values
  neg_data <- data.frame(x = -5:-1, y = -10:-6)
  p_neg <- ggplot(neg_data, aes(x, y)) + geom_point()
  
  result <- axis_limits(p_neg, limits_expand = list(x = 0, y = 0))
  expect_s3_class(result, c("gg", "ggplot"))
})

test_that("axis_limits preserves ggplot structure", {
  # Setup test data
  data("cars", package = "datasets")
  library(ggplot2)
  p <- ggplot(cars, aes(speed, dist)) + 
    geom_point() + 
    ggtitle("Test Plot") +
    theme_minimal()
  
  # Test that modifications preserve original structure
  result <- axis_limits(p, limits_equal = TRUE)
  
  # Should still be a ggplot
  expect_s3_class(result, c("gg", "ggplot"))
  
  # Should preserve original data
  expect_identical(result$data, p$data)
  
  # Should preserve original mapping
  expect_identical(result$mapping, p$mapping)
  
  # Should have at least the original layers plus any added
  expect_true(length(result$layers) >= length(p$layers))
})

test_that("axis_limits handles various expand limit combinations", {
  data("cars", package = "datasets")
  library(ggplot2)
  p <- ggplot(cars, aes(speed, dist)) + geom_point()
  
  # Test single value expansion
  result1 <- axis_limits(p, limits_expand = list(x = 0))
  expect_s3_class(result1, c("gg", "ggplot"))
  
  # Test vector expansion
  result2 <- axis_limits(p, limits_expand = list(y = c(-10, 150)))
  expect_s3_class(result2, c("gg", "ggplot"))
  
  # Test mixed single and vector
  result3 <- axis_limits(p, limits_expand = list(x = 0, y = c(-10, 150)))
  expect_s3_class(result3, c("gg", "ggplot"))
  
  # Test that values are properly sorted internally
  result4 <- axis_limits(p, limits_expand = list(x = c(100, -10))) # reversed order
  expect_s3_class(result4, c("gg", "ggplot"))
})