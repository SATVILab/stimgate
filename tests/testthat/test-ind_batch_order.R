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
  expect_type(exampleData$batchList, "list")
  expect_equal(length(exampleData$batchList), 2)
  # For index 1: idxUnstim = 1, idxStim = 2 -> c(1, 2)
  expect_equal(exampleData$batchList[[1]], c(1, 2))
  # For index 2: idxUnstim = 3, idxStim = 4 -> c(3, 4)
  expect_equal(exampleData$batchList[[2]], c(3, 4))
})
