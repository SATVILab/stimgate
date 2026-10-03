test_that("invalid FCS gates preserve existing output", {
  output <- tempfile("fcs-existing-")
  dir.create(output)
  withr::defer(unlink(output, recursive = TRUE))
  sentinel <- file.path(output, "previous.txt")
  writeLines("previous export", sentinel)

  expect_error(writeStimFCS(
    pathProject = tempfile("missing-project-"),
    .data = "invalid", indBatchList = list(c(1, 2)),
    pathDirSave = output, chnl = "missing"
  ))
  expect_true(file.exists(sentinel))
  expect_identical(readLines(sentinel), "previous export")
})
