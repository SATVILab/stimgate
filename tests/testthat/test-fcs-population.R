test_that("FCS export selects events from the requested population", {
  example <- getExampleData()
  withr::defer(unlink(dirname(example$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(example$pathGs)[1]
  chnl <- example$chnl
  subset_gate <- flowCore::rectangleGate(
    .gate = stats::setNames(list(c(0.5, Inf)), chnl[[1]]),
    filterId = "selected"
  )
  flowWorkspace::gs_pop_add(gs, subset_gate, parent = "root")
  flowWorkspace::recompute(gs)
  expected <- flowCore::exprs(
    flowWorkspace::gh_pop_get_data(gs[[1]], "selected")
  )
  expect_gt(nrow(expected), 0)
  expect_lt(nrow(expected), nrow(flowCore::exprs(
    flowWorkspace::gh_pop_get_data(gs[[1]], "root")
  )))
  gates <- data.frame(
    chnl = chnl, marker = example$marker,
    batch = "batch_1", ind = "1", gate = -Inf
  )
  output <- tempfile("fcs-population-")
  withr::defer(unlink(output, recursive = TRUE))
  writeStimFCS(tempdir(), gs,
    pop = "selected", indBatchList = list(1L),
    pathDirSave = output, chnl = chnl, gateTbl = gates,
    gateTypeCytPos = "base"
  )
  files <- list.files(output, full.names = TRUE)
  expect_length(files, 1L)
  actual <- flowCore::exprs(
    flowCore::read.FCS(files[[1]], transformation = FALSE)
  )
  expect_identical(dim(actual), dim(expected))
  expect_identical(colnames(actual), colnames(expected))
  expect_equal(as.numeric(actual), as.numeric(expected))

  project <- tempfile("fcs-pop-project-")
  withr::defer(unlink(project, recursive = TRUE))
  for (channel in chnl) {
    path <- .gatesGetPathAll(project, "selected", channel, FALSE)
    dir.create(dirname(path), recursive = TRUE)
    saveRDS(gates[gates$chnl == channel, ], path)
  }
  writeStimFCS(project, gs,
    indBatchList = list(1L), pathDirSave = output,
    chnl = chnl, gateTypeCytPos = "base"
  )
  actual <- flowCore::exprs(
    flowCore::read.FCS(files[[1]], transformation = FALSE)
  )
  expect_identical(dim(actual), dim(expected))
  expect_identical(colnames(actual), colnames(expected))
  expect_equal(as.numeric(actual), as.numeric(expected))
})
