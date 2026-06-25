library(testthat)
exampleData <- getExampleData()
gs <- flowWorkspace::load_gs(exampleData$path_gs)
pathProject <- file.path(dirname(exampleData$path_gs), "stimgate")

# First run gating to create necessary gate data
invisible(gateStim(
  .data = gs,
  pathProject = pathProject,
  popGate = "root",
  batchList = exampleData$batchList,
  chnl = exampleData$chnl
))

test_that("stimgate_plot function exists", {
  # Just test that the function exists and is callable
  expect_true(exists("stimgate_plot", envir = asNamespace("stimgate")))
  expect_true(is.function(stimgate_plot))
})

test_that("stimgate_plot runs", {
  # debugonce(.plotGateBv)
  p <- stimgate_plot(
    ind = exampleData$batchList[[1]], # indices in `gs` to plot
    .data = gs, # GatingSet
    pathProject = pathProject,
    chnl = exampleData$chnl,
    grid = TRUE
  )
  expect_true(inherits(p, "ggplot"))
})

test_that("stimgate_plot returns NULL when pList is empty", {
  # Test with empty ind list to generate empty pList
  result <- stimgate_plot(
    ind = list(),
    .data = gs,
    pathProject = pathProject,
    chnl = exampleData$chnl,
    grid = TRUE
  )
  expect_null(result)
})

test_that(".plotGateBv returns NULL for single marker", {
  # Test single marker scenario
  singleMarker <- exampleData$chnl[1]
  result <- stimgate:::.plotGateBv(
    chnl = singleMarker,
    marker = NULL,
    ind = exampleData$batchList[[1]],
    indLab = NULL,
    .data = gs,
    axisLab = NULL,
    pathProject = pathProject,
    excMin = TRUE,
    limitsExpand = NULL,
    limitsEqual = FALSE,
    showGate = TRUE,
    minCell = 10
  )
  expect_null(result)
})

test_that("plot functions handle minCell threshold correctly", {
  # Test with very high minCell to trigger early return
  # debugonce(.getExNew)
  result <- stimgate:::.plotGateBv(
    chnl = exampleData$chnl,
    pop = "root",
    marker = NULL,
    ind = exampleData$batchList[[1]],
    indLab = NULL,
    .data = gs,
    axisLab = NULL,
    pathProject = pathProject,
    excMin = TRUE,
    limitsExpand = NULL,
    limitsEqual = FALSE,
    showGate = TRUE,
    chnlGate = NULL,
    markerGate = NULL,
    bias = FALSE,
    combnExc = NULL,
    gateTypeCytPos = "cyt",
    gateTypeSinglePos = "single",
    mult = FALSE,
    gateUnsMethod = "min",
    minCell = 999999 # Very high threshold
  )
  # Should return a list with NULLs filtered out, or NULL
  expect_true(is.null(result) || (is.list(result) && length(result) == 0))
})

test_that(".plotGetLab handles various valLab configurations", {
  # Test with NULL valLab
  result1 <- stimgate:::.plotGetLab(
    val = c("A", "B"),
    valLab = NULL,
    i = NULL
  )
  expect_equal(result1, c("A", "B"))

  # Test with named valLab
  result2 <- stimgate:::.plotGetLab(
    val = c("A", "B"),
    valLab = c("A" = "Label A", "B" = "Label B"),
    i = NULL
  )
  expect_equal(result2, c("Label A", "Label B"))
  expect_null(names(result2))

  # Test with unnamed valLab and i NULL
  result3 <- stimgate:::.plotGetLab(
    val = c("A", "B"),
    valLab = c("Label A", "Label B"),
    i = NULL
  )
  expect_equal(result3, c("Label A", "Label B"))
  expect_null(names(result3))

  # Test with unnamed valLab and i non-null
  result4 <- stimgate:::.plotGetLab(
    val = c("A", "B"),
    valLab = c("Label A", "Label B"),
    i = 1
  )
  expect_equal(result4, "Label A")
  expect_null(names(result4))
})

test_that(".plotGateUv returns NULL when all markers return NULL", {
  # Test with empty ind to generate NULL results
  result <- stimgate:::.plotGateUv(
    ind = list(),
    indLab = NULL,
    .data = gs,
    chnl = exampleData$chnl,
    marker = NULL,
    excMin = TRUE,
    axisLab = NULL,
    showGate = TRUE,
    pathProject = pathProject,
    minCell = 10
  )
  expect_null(result)
})

test_that(".plotGateUvMarker returns NULL when plotTbl is NULL", {
  # Test with empty ind list to generate NULL plotTbl
  result <- stimgate:::.plotGateUvMarker(
    chnl = exampleData$chnl[1],
    marker = NULL,
    pop = "root",
    ind = list(),
    .data = gs,
    excMin = TRUE,
    indLab = NULL,
    axisLab = NULL,
    showGate = TRUE,
    pathProject = pathProject,
    minCell = 10
  )
  expect_null(result)
})

test_that(".plotGateUvMarkerGetPlotTbl returns NULL for insufficient cells", {
  result <- stimgate:::.plotGateUvMarkerGetPlotTbl(
    ind = exampleData$batchList[[1]],
    .data = gs,
    chnl = exampleData$chnl[1],
    marker = NULL,
    pop = "root",
    excMin = TRUE,
    indLab = NULL,
    pathProject = pathProject,
    bias = FALSE,
    combnExc = NULL,
    chnlGate = NULL,
    markerGate = NULL,
    gateTypeCytPos = "cyt",
    gateTypeSinglePos = "single",
    mult = FALSE,
    gateUnsMethod = "min",
    minCell = 999999 # Very high threshold
  )
  expect_null(result)
})

test_that(".plotGateUvMarkerPlotInit handles different condition branches", {
  # Create mock plotTbl for testing
  plotTbl <- data.frame(
    x = 1:10,
    y = 1:10,
    type = rep(c("raw", "adj"), 5),
    indLab = rep(c("Sample1", "Sample2"), 5)
  )

  # Test excMin = TRUE, multiple ind
  p1 <- stimgate:::.plotGateUvMarkerPlotInit(
    plotTbl = plotTbl,
    excMin = TRUE,
    ind = c(1, 2),
    indLab = c("Sample1", "Sample2")
  )
  expect_s3_class(p1, "ggplot")

  # Test excMin = TRUE, single ind
  p2 <- stimgate:::.plotGateUvMarkerPlotInit(
    plotTbl = plotTbl,
    excMin = TRUE,
    ind = c(1),
    indLab = c("Sample1")
  )
  expect_s3_class(p2, "ggplot")

  # Test excMin = FALSE, multiple ind
  p3 <- stimgate:::.plotGateUvMarkerPlotInit(
    plotTbl = plotTbl,
    excMin = FALSE,
    ind = c(1, 2),
    indLab = c("Sample1", "Sample2")
  )
  expect_s3_class(p3, "ggplot")

  # Test excMin = FALSE, single ind
  p4 <- stimgate:::.plotGateUvMarkerPlotInit(
    plotTbl = plotTbl,
    excMin = FALSE,
    ind = c(1),
    indLab = c("Sample1")
  )
  expect_s3_class(p4, "ggplot")
})

test_that(".plotGrid returns pList when plot = FALSE", {
  # Create mock plot list
  pList <- list(
    plot1 = ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = 1, y = 1)),
    plot2 = ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = 2, y = 2))
  )

  # Test with plot = FALSE
  result <- stimgate:::.plotGrid(
    plot = FALSE,
    pList = pList,
    nCol = 2
  )
  expect_identical(result, pList)

  # Test with plot = TRUE (should return combined plot)
  result2 <- stimgate:::.plotGrid(
    plot = TRUE,
    pList = pList,
    nCol = 2
  )
  expect_s3_class(result2, "ggplot")
})

# Additional comprehensive edge case tests
test_that("comprehensive edge case coverage for plot_gate functions", {
  # Test .plotGateUvMarkerGetPlotTblInd with insufficient cells
  expect_null(stimgate:::.plotGateUvMarkerGetPlotTblInd(
    ind = c(1), # Single index
    .data = gs,
    chnl = exampleData$chnl[1],
    marker = NULL,
    pop = "root",
    excMin = TRUE,
    pathProject = pathProject,
    bias = FALSE,
    combnExc = NULL,
    chnlGate = NULL,
    markerGate = NULL,
    gateTypeCytPos = "cyt",
    gateTypeSinglePos = "single",
    mult = FALSE,
    gateUnsMethod = "min",
    minCell = 999999 # Impossible threshold
  ))

  # Test .plotGetExTbl with single index
  exTbl <- stimgate:::.plotGetExTbl(
    ind = c(1),
    .data = gs,
    chnl = exampleData$chnl[1],
    marker = NULL,
    pop = "root",
    excMin = TRUE,
    pathProject = pathProject,
    bias = FALSE,
    combnExc = NULL,
    chnlGate = NULL,
    markerGate = NULL,
    gateTypeCytPos = "cyt",
    gateTypeSinglePos = "single",
    mult = FALSE,
    gateUnsMethod = "min"
  )
  expect_true(is.data.frame(exTbl))
  expect_true(nrow(exTbl) > 0)

  # Test .plotAddAxisTitle with single and multiple markers
  pBase <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = 1, y = 1))
  pSingle <- stimgate:::.plotAddAxisTitle(
    pBase,
    exampleData$chnl[1],
    NULL,
    NULL
  )
  expect_s3_class(pSingle, "ggplot")

  pDouble <- stimgate:::.plotAddAxisTitle(
    pBase,
    exampleData$chnl,
    NULL,
    NULL
  )
  expect_s3_class(pDouble, "ggplot")

  # Test .plotAddTitle
  pTitled <- stimgate:::.plotAddTitle(pBase, c(1, 2), 1, NULL)
  expect_s3_class(pTitled, "ggplot")

  # Test .plotAddGate when showGate = FALSE
  pNoGate <- stimgate:::.plotAddGate(
    pBase,
    gs,
    c(1),
    exampleData$chnl[1],
    pathProject,
    showGate = FALSE
  )
  expect_identical(pNoGate, pBase)
})

test_that("test plot_cyto import and dependencies", {
  # Test that plot_cyto function is accessible
  expect_true(exists("plot_cyto", envir = asNamespace("stimgate")))

  # Test hexbin namespace checking
  hexbinAvailable <- requireNamespace("hexbin", quietly = TRUE)
  expect_true(is.logical(hexbinAvailable))
})

# if (dir.exists(exampleData$path_gs)) {
#   unlink(exampleData$path_gs, recursive = TRUE)
# }
