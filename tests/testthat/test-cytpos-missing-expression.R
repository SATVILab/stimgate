test_that("missing expression does not erase known other-marker positivity", {
  ex <- data.frame(A = c(0, NA, 2), B = c(2, 2, 0), C = c(NA, 0, NA))
  gates <- tibble::tibble(chnl = c("A", "B", "C"), gate = 1)
  base <- .getCytPosBasePos(ex, gates)
  expect_identical(base$nPos, c(1L, 1L, 1L))
  expect_identical(base$nPos - as.integer(base$pos$A) > 0L, c(TRUE, TRUE, FALSE))
  expect_identical(base$nPos - as.integer(base$pos$B) > 0L, c(FALSE, FALSE, TRUE))
})
