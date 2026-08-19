test_that("cytPos gates actually happen", {
  set.seed(42L)

  probVecUns <- c(0.90, 0.04, 0.04, 0.02)
  probResponseVec <- list(c(-0.03, 0.01, 0.01, 0.02))
  meanExprMat <- matrix(
    c(0, 0, 0, 4, 4, 0, 4, 4),
    nrow = 4L,
    byrow = TRUE
  )

  simRes <- simCytExperiment(
    nSample = 1L,
    nMarker = 2L,
    nCondition = 2L,
    nCluster = 4L,
    nCellByCondition = 1e4,
    transformationFunc = function(x) x,
    mixtureType = "gaussianOnly",
    meanExprMat = meanExprMat,
    clusterLabelVec = c("negNeg", "negPos", "posNeg", "posPos"),
    probVecUns = probVecUns,
    probResponseVecByStimCondition = probResponseVec,
    clusterPerturbationSd = 0
  )

  ffList <- lapply(simRes$flowFrameList, function(ff) {
    flowCore::colnames(ff) <- c("BC1(La139)Dd", "BC2(Pr141)Dd")
    ff
  })
  gs <- flowWorkspace::GatingSet(flowCore::flowSet(ffList))
  pathProject <- file.path(tempdir(), "stimgate_cytpos_test")
  on.exit(unlink(pathProject, recursive = TRUE), add = TRUE)

  gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = list(c(1L, 2L)),
    marker = c("MarkerF1", "MarkerF2")
  )

  gateTbl <- getStimGates(pathProject)
  expect_true(is.data.frame(gateTbl))
  expect_true(nrow(gateTbl) > 0)
  expect_true(all(c("chnl", "marker", "batch", "ind", "gate", "gateCyt") %in% colnames(gateTbl)))
  expect_true(is.numeric(gateTbl$gate))
  expect_true(is.numeric(gateTbl$gateCyt))
  expect_true(all(is.finite(gateTbl$gate)))
  expect_true(all(is.finite(gateTbl$gateCyt)))
  expect_true(any(gateTbl$gateCyt != gateTbl$gate, na.rm = TRUE))
  expect_equal(nrow(gateTbl), 2L)
})