local({
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "scripts", "r", "coexpression-gates.R"), env)

  testthat::test_that("positivity follows the lowered gate only above the raised cut", {
    ex <- data.frame(A = c(3, 3, 1.5, 0), B = c(2, 0.5, 2, 5), C = c(0, 0, 0, 0))
    gate <- c(A = 2, B = 3, C = NA)
    low <- data.frame(a = "A", b = "B", cut = 1, condCut = 2.5, lowered = TRUE)
    pos <- env$.coexPositive(ex, gate, low)
    testthat::expect_equal(unname(pos[, "A"]), c(TRUE, TRUE, FALSE, FALSE))
    testthat::expect_equal(unname(pos[, "B"]), c(TRUE, FALSE, FALSE, TRUE))
    # A marker without a finite gate is negative throughout.
    testthat::expect_false(any(pos[, "C"]))
    testthat::expect_equal(env$.coexPositive(ex, gate)[, "B"], c(FALSE, FALSE, FALSE, TRUE))
  })

  testthat::test_that("combination counts cover every combination once", {
    pos <- cbind(A = c(TRUE, TRUE, FALSE), B = c(TRUE, FALSE, FALSE))
    n <- env$.coexCombnCounts(pos)
    testthat::expect_equal(n[["A+B+"]], 1L)
    testthat::expect_equal(n[["A+B-"]], 1L)
    testthat::expect_equal(n[["A-B-"]], 1L)
    testthat::expect_equal(sum(n), 3L)
    testthat::expect_length(n, 4L)
  })

  testthat::test_that("every ordered pair is evaluated and only co-expressed pairs move", {
    withr::local_seed(5)
    neg <- function(n) data.frame(A = abs(rnorm(n, 0.4, 0.3)), B = abs(rnorm(n, 0.6, 0.4)),
      C = abs(rnorm(n, 0.5, 0.3)))
    # A and B co-expressed by responders; C responds alone.
    resp <- data.frame(A = rnorm(200, 3.5, 0.4), B = runif(200, 2, 5), C = abs(rnorm(200, 0.5, 0.3)))
    c_only <- data.frame(A = abs(rnorm(200, 0.4, 0.3)), B = abs(rnorm(200, 0.6, 0.4)), C = rnorm(200, 4, 0.3))
    dat <- list(stim = rbind(neg(20000), resp, c_only), uns = neg(20000),
      gate = c(A = 2, B = 4.5, C = 2.5))
    low <- env$.coexLowerGates(dat)
    testthat::expect_equal(nrow(low), 6L)
    testthat::expect_true(low$lowered[low$a == "A" & low$b == "B"])
    testthat::expect_false(any(low$lowered[low$a == "C" | low$b == "C"]))
  })
})
