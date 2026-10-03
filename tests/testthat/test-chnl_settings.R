test_that("automatic settings use the same evenly spaced batches", {
  withr::local_preserve_seed()
  expect_identical(.completeChnlSettingsBatchInd(as.list(1:3)), c(1, 2, 3))
  expect_identical(
    .completeChnlSettingsBatchInd(as.list(1:9)),
    c(1, 3, 5, 7, 9)
  )
  expect_length(.completeChnlSettingsBatchInd(list()), 0L)
  set.seed(1)
  a <- .completeChnlSettingsBatchInd(as.list(1:20))
  set.seed(2)
  expect_identical(.completeChnlSettingsBatchInd(as.list(1:20)), a)
})
