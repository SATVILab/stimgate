# Share one expensive fixture within this file and clean it on exit.
local({
exampleData <- getExampleData()
withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
gs <- flowWorkspace::load_gs(exampleData$pathGs)
pathProject <- file.path(dirname(exampleData$pathGs), "stimgate")

# First run gating to create necessary gate data
invisible(gateStim(
  .data = gs,
  pathProject = pathProject,
  popGate = "root",
  batchList = exampleData$batchList,
  chnl = exampleData$chnl
))

test_that("plotStim runs", {
  p <- plotStim(
    ind = exampleData$batchList[[1]], # indices in `gs` to plot
    .data = gs, # GatingSet
    pathProject = pathProject,
    chnl = exampleData$chnl,
    grid = TRUE
  )
  expect_true(inherits(p, "ggplot"))
})

test_that("plotStim returns NULL when pList is empty", {
  # Test with empty ind list to generate empty pList
  result <- plotStim(
    ind = list(),
    .data = gs,
    pathProject = pathProject,
    chnl = exampleData$chnl,
    grid = TRUE
  )
  expect_null(result)
})

test_that("plotStim returns only univariate plots for a single channel", {
  pList <- plotStim(
    ind = exampleData$batchList[[1]],
    .data = gs,
    pathProject = pathProject,
    chnl = exampleData$chnl[1],
    grid = FALSE
  )
  expect_length(pList, 1L)
  expect_s3_class(pList[[1]]$layers[[1]]$geom, "GeomLine")
})

test_that("plotStim returns NULL when no sample has minCell cells", {
  result <- plotStim(
    ind = exampleData$batchList[[1]],
    .data = gs,
    pathProject = pathProject,
    chnl = exampleData$chnl,
    grid = FALSE,
    minCell = 999999
  )
  expect_null(result)
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

test_that(".plotGateUvMarker returns NULL when ind is empty", {
  result <- stimgate:::.plotGateUvMarker(
    ind = list(),
    indLab = NULL,
    marker = NULL,
    chnl = exampleData$chnl[1],
    pop = "root",
    excMin = TRUE,
    axisLab = NULL,
    showGate = TRUE,
    pathProject = pathProject,
    minCell = 10,
    exArgs = list()
  )
  expect_null(result)
})

test_that("univariate plots map alpha to excMin and colour to sample", {
  aesNames <- function(ind, excMin) {
    pList <- plotStim(
      ind = ind,
      .data = gs,
      pathProject = pathProject,
      chnl = exampleData$chnl[1],
      excMin = excMin,
      grid = FALSE
    )
    sort(names(pList[[1]]$mapping))
  }
  indVec <- exampleData$batchList[[1]]
  expect_identical(aesNames(indVec, TRUE), c("alpha", "colour", "x", "y"))
  expect_identical(aesNames(indVec[[2]], TRUE), c("alpha", "x", "y"))
  expect_identical(aesNames(indVec, FALSE), c("colour", "x", "y"))
  expect_identical(aesNames(indVec[[2]], FALSE), c("x", "y"))
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
test_that("plot helpers handle axis titles and disabled gates", {
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

  pNoGate <- stimgate:::.plotAddGate(
    pBase,
    ind = 1,
    marker = NULL,
    chnl = exampleData$chnl[1],
    pop = "root",
    pathProject = pathProject,
    showGate = FALSE
  )
  expect_identical(pNoGate, pBase)
})

test_that("plotStim keeps univariate plots when excMin = FALSE", {
  pList <- plotStim(
    ind = exampleData$batchList[[1]],
    .data = gs,
    pathProject = pathProject,
    chnl = exampleData$chnl[1],
    excMin = FALSE,
    grid = FALSE
  )
  expect_length(pList, 1L)
  expect_s3_class(pList[[1]], "ggplot")
  expect_no_warning(layerTbl <- ggplot2::layer_data(pList[[1]], 1L))
  expect_gt(nrow(layerTbl), 0L)
})

test_that("bivariate gate lines are drawn on the axis of their own channel", {
  pathProjectCopy <- withr::local_tempdir()
  file.copy(pathProject, pathProjectCopy, recursive = TRUE)
  pathProjectCopy <- file.path(pathProjectCopy, basename(pathProject))
  # x channel has no gate; only the y channel is gated
  unlink(
    file.path(
      pathProjectCopy, "gates", "poproot", paste0("chnl", exampleData$chnl[1])
    ),
    recursive = TRUE
  )
  pList <- plotStim(
    ind = exampleData$batchList[[1]][[2]],
    .data = gs,
    pathProject = pathProjectCopy,
    chnl = exampleData$chnl,
    grid = FALSE
  )
  geomVec <- vapply(
    pList[[1]]$layers,
    function(l) class(l$geom)[[1]],
    character(1)
  )
  expect_false("GeomVline" %in% geomVec)
  expect_true("GeomHline" %in% geomVec)
})
})
