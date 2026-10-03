test_that("inferred FCS channels match explicit channels with exclusions", {
  example <- getExampleData()
  withr::defer(unlink(dirname(example$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(example$pathGs)[1]
  chnl <- example$chnl
  gates <- data.frame(chnl = chnl, marker = example$marker,
                     batch = "batch_1", ind = "1", gate = 0.5)
  output <- tempfile("fcs-channels-")
  withr::defer(unlink(output, recursive = TRUE))
  writeStimFCS(tempdir(), gs, indBatchList = list(1L), pathDirSave = output,
               chnl = chnl, gateTbl = gates, gateTypeCytPos = "base",
               combnExc = list(chnl[[1]]))
  file <- list.files(output, full.names = TRUE)[[1]]
  expected <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE))
  writeStimFCS(tempdir(), gs, indBatchList = list(1L), pathDirSave = output,
               gateTbl = gates, gateTypeCytPos = "base",
               combnExc = list(chnl[[1]]))
  actual <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE))
  expect_identical(actual, expected)
})

test_that("inferred FCS channels support a single event", {
  example <- getExampleData()
  withr::defer(unlink(dirname(example$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(example$pathGs)[1]
  fr <- flowWorkspace::gh_pop_get_data(gs[[1]])
  if (inherits(fr, "cytoframe")) {
    fr <- flowWorkspace::cytoframe_to_flowFrame(fr)
  }
  fr <- fr[1, ]
  gs <- flowWorkspace::GatingSet(flowCore::flowSet(fr))
  gates <- data.frame(chnl = example$chnl, marker = example$marker,
                     batch = "batch_1", ind = "1", gate = -Inf)
  output <- tempfile("fcs-one-event-")
  withr::defer(unlink(output, recursive = TRUE))
  writeStimFCS(tempdir(), gs, indBatchList = list(1L), pathDirSave = output,
               gateTbl = gates, gateTypeCytPos = "base")
  files <- list.files(output, full.names = TRUE)
  expect_length(files, 1L)
  expect_identical(nrow(flowCore::exprs(flowCore::read.FCS(files[[1]]))), 1L)
})
