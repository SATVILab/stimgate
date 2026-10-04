test_that("getBatchList positions unstimulated sample index first", {
  fnTblInfo <- data.frame(
    PID = c("P1", "P1", "P2", "P2"),
    Stim = c("SEB", "UNS", "UNS", "PMA"),
    NCells = c(1000, 1000, 1000, 1000),
    stringsAsFactors = FALSE
  )

  batchList <- getBatchList(
    fnTblInfo = fnTblInfo,
    colGrp = "PID",
    colStim = "Stim",
    unsChr = "UNS",
    colNCell = "NCells",
    minCell = 100
  )

  expect_named(batchList, c("P1", "P2"))
  # Batch P1: UNS is row 2, SEB is row 1 -> expected c(2, 1)
  expect_equal(batchList$P1, c(2, 1))
  # Batch P2: UNS is row 3, PMA is row 4 -> expected c(3, 4)
  expect_equal(batchList$P2, c(3, 4))
})

test_that("getExampleData positions unstimulated sample index first in batchList", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  expect_type(exampleData$batchList, "list")
  expect_equal(length(exampleData$batchList), 2)
  # For index 1: idxUnstim = 1, idxStim = 2 -> c(1, 2)
  expect_equal(exampleData$batchList[[1]], c(1, 2))
  # For index 2: idxUnstim = 3, idxStim = 4 -> c(3, 4)
  expect_equal(exampleData$batchList[[2]], c(3, 4))
})

test_that("getBatchList errors when a group has several unstimulated samples", {
  fnTblInfo <- data.frame(
    PID = c("P1", "P1", "P1"),
    Stim = c("UNS", "SEB", "UNS"),
    NCells = c(1000, 1000, 1000)
  )
  expect_error(
    getBatchList(
      fnTblInfo = fnTblInfo, colGrp = "PID", colStim = "Stim",
      unsChr = "UNS", colNCell = "NCells", minCell = 100
    ),
    "2 unstimulated samples"
  )
})

test_that("batchList validation enforces the unstim-first convention", {
  # a shared unstim that is first in every batch is allowed
  expect_true(.verifyBatchList(list(c(1, 2), c(1, 3))))
  expect_true(.verifyBatchList(list(c("u", "s1", "s2"))))
  expect_error(.verifyBatchList(list(c(1, 2), 3)), "at least two")
  expect_error(.verifyBatchList(list(c(1, NA))), "at least two")
  expect_error(.verifyBatchList(list(c(1, 2, 2))), "more than once")
  expect_error(.verifyBatchList(list(c(1, 1, 2))), "both unstimulated")
  expect_error(
    .verifyBatchList(list(c(1, 2), c(2, 3))),
    "both unstimulated"
  )
  expect_error(
    .verifyBatchList(list(c(1, 2), c(3, 2))),
    "more than once"
  )
})

test_that(".getIndUns returns the first sample of the containing batch", {
  indBatchList <- list(a = c(1, 2, 4), b = c(1, 3), c = c(5, 6))
  expect_identical(.getIndUns(4, indBatchList), 1)
  expect_identical(.getIndUns(1, indBatchList), 1)
  expect_identical(.getIndUns(6, indBatchList), 5)
  expect_identical(.getIndUns("6", indBatchList), 5)
})
