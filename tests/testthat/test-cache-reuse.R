test_that("gateStim reads GatingSet population exactly once per sample", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "read_once")
  withr::defer(unlink(pathProject, recursive = TRUE))

  nReads <- 0L
  origGetData <- flowWorkspace::gh_pop_get_data
  testthat::local_mocked_bindings(
    gh_pop_get_data = function(x, y, ...) {
      if (!missing(y) && !is.null(y)) {
        nReads <<- nReads + 1L
        origGetData(x, y = y, ...)
      } else {
        origGetData(x, ...)
      }
    },
    .package = "flowWorkspace"
  )

  invisible(gateStim(
    pathProject = pathProject,
    .data = gs,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker,
    parallel = FALSE
  ))

  expectedReads <- length(unique(unlist(exampleData$batchList)))
  expect_equal(nReads, expectedReads)
})

test_that("gateStim invalidates stale expression cache on rerun", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "stale_cache")
  withr::defer(unlink(pathProject, recursive = TRUE))

  set.seed(1)
  invisible(gateStim(
    pathProject = pathProject,
    .data = gs,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  ))
  gates1 <- getStimGates(pathProject)

  cachedFiles <- list.files(
    file.path(pathProject, "sampleData", "pop_root"),
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )
  expect_gt(length(cachedFiles), 0L)
  exCorrupt <- readRDS(cachedFiles[[1L]]) + 100
  saveRDS(exCorrupt, cachedFiles[[1L]])

  set.seed(1)
  invisible(gateStim(
    pathProject = pathProject,
    .data = gs,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  ))
  gates2 <- getStimGates(pathProject)

  expect_equal(gates2, gates1)
})

test_that("rerun with fewer markers removes gates for omitted markers", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "fewer_markers")
  withr::defer(unlink(pathProject, recursive = TRUE))

  markers2 <- exampleData$marker[1:2]
  invisible(gateStim(
    pathProject = pathProject,
    .data = gs,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = markers2
  ))
  gatesBefore <- getStimGates(pathProject)
  expect_equal(sort(unique(gatesBefore$marker)), sort(markers2))

  marker1 <- markers2[1L]
  invisible(gateStim(
    pathProject = pathProject,
    .data = gs,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = marker1
  ))
  gatesAfter <- getStimGates(pathProject)
  expect_equal(unique(gatesAfter$marker), marker1)
})

test_that("gateStim results are bit-for-bit deterministic across projects", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  path1 <- file.path(dirname(exampleData$pathGs), "det1")
  path2 <- file.path(dirname(exampleData$pathGs), "det2")
  withr::defer(unlink(c(path1, path2), recursive = TRUE))

  set.seed(1)
  invisible(gateStim(
    pathProject = path1,
    .data = gs,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  ))
  gates1 <- getStimGates(path1)

  set.seed(1)
  invisible(gateStim(
    pathProject = path2,
    .data = gs,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  ))
  gates2 <- getStimGates(path2)

  expect_identical(gates1, gates2)
})

test_that(".completeChnlSettingsBwShared memoises readBatch calls", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "bw_shared_memo")
  withr::defer(unlink(pathProject, recursive = TRUE))

  invisible(gateStim(
    pathProject = pathProject,
    .data = gs,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker[1L]
  ))
  settingsList <- stimgateMetaReadSettingsChnls(pathProject)
  chnlSettings <- settingsList[[1L]]
  chnlSettings$bwScope <- "cytokine"
  chnlSettings$bw <- NULL

  calls <- integer()
  origReadBatch <- .bwSharedReadBatch
  testthat::local_mocked_bindings(
    .bwSharedReadBatch = function(i, ...) {
      calls <<- c(calls, i)
      origReadBatch(i, ...)
    }
  )

  .completeChnlSettingsBwShared(
    chnlSettings = chnlSettings,
    indBatchList = exampleData$batchList,
    .data = gs,
    pathProject = pathProject
  )

  expect_equal(anyDuplicated(calls), 0L)
  expect_equal(length(calls), length(exampleData$batchList))
})

test_that(".getCpCluster accepts precomputed exLookup", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "cp_cluster_ex")
  withr::defer(unlink(pathProject, recursive = TRUE))

  invisible(gateStim(
    pathProject = pathProject,
    .data = gs,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker[1L]
  ))
  settingsList <- stimgateMetaReadSettingsChnls(pathProject)
  chnlSettings <- settingsList[[1L]]

  gatePath <- .gatesGetPathAll(
    pathProject = pathProject,
    pop = chnlSettings$popGate,
    chnlCut = chnlSettings$chnlCut,
    init = TRUE
  )
  gateTbl <- readRDS(gatePath)

  exLookup <- .getCpClusterLocExLookup(
    .data = gs,
    indBatchList = exampleData$batchList,
    chnlSettings = chnlSettings,
    filterOtherCytPos = FALSE,
    calcCytPosGates = FALSE,
    gateTbl = gateTbl,
    pathProject = pathProject
  )

  lookupCalls <- 0L
  testthat::local_mocked_bindings(
    .getCpClusterLocExLookup = function(...) {
      lookupCalls <<- lookupCalls + 1L
      list()
    }
  )

  res <- .getCpCluster(
    .data = gs,
    gateTbl = gateTbl,
    chnlSettings = chnlSettings,
    stage = "init",
    pathProject = pathProject,
    filterOtherCytPos = FALSE,
    calcCytPosGates = FALSE,
    indBatchList = exampleData$batchList,
    exLookup = exLookup
  )

  expect_equal(lookupCalls, 0L)
  expect_true(is.data.frame(res))
})

test_that(".completeChnlSettingsInd reads batch expression once for limits", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "complete_chnl_ex")
  withr::defer(unlink(pathProject, recursive = TRUE))

  invisible(gateStim(
    pathProject = pathProject,
    .data = gs,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker[1L]
  ))
  chnlLab <- stimgateMetaReadChnlLab(pathProject)
  chnlCurr <- names(chnlLab)[1L]

  readCount <- 0L
  origGetExList <- .getExList
  testthat::local_mocked_bindings(
    .getExList = function(...) {
      readCount <<- readCount + 1L
      origGetExList(...)
    }
  )

  commonSettings <- list(
    popGate = "root",
    biasUns = NULL,
    biasUnsFactor = 1,
    excMin = TRUE,
    cpMin = NULL,
    bwMin = "auto",
    bw = NULL,
    bwMax = "auto",
    bwFallback = "auto",
    bwMtd = "nrd0",
    bwAdj = 1,
    bwNcellMin = 100,
    bwNcellMax = 1000,
    bwCluster = FALSE,
    bwScope = "sample",
    bwAdaptive = FALSE,
    bwAdaptiveDensityN = NULL,
    bwAdaptivePadFrac = 0.15,
    bwAdaptiveCore = NULL,
    bwAdaptiveExtra = NULL,
    bwAdaptiveCrossover = NULL,
    bwAdaptiveTransitionWidth = 0,
    normPeakMinRel = 0.75,
    normExtraFrac = 0.2,
    normExtraMax = Inf,
    normLambda = seq(-2, 2, length.out = 81),
    normDensityN = 512L,
    normExcessBwMtd = "hpi3",
    normExcessNcell = 10000L,
    normAdaptiveNcell = 2500L,
    normMtd = "moments",
    minCell = 100,
    tolClust = NULL,
    locProbCol = "pred",
    locMinPeakProb = 0.05,
    locEnforceShapeThreshold = TRUE,
    locDipAlpha = 0.05,
    locAntimodeHeightFrac = 0.1,
    locAntimodeLowRel = 0.5,
    locAntimodeLowAbs = 0,
    locFlatDerivFrac = 0.01,
    locFlatHardDerivFrac = 0.05,
    locMarginalPurityRel = 0.5,
    locMarginalCellBinRatio = 0.5,
    locMarginalRefQuantile = 0.95,
    gateCombn = "min",
    maxPosProbX = Inf,
    gateQuant = c(0.25, 0.75)
  )

  .completeChnlSettingsInd(
    chnlSettingsCommon = commonSettings,
    chnlSettingsSpec = list(chnlCut = chnlCurr, marker = chnlLab[[chnlCurr]]),
    chnl = chnlCurr,
    .data = gs,
    indBatchList = exampleData$batchList,
    pathProject = pathProject
  )

  nBatchesExpected <- length(
    .completeChnlSettingsBatchInd(exampleData$batchList)
  )
  expect_equal(readCount, nBatchesExpected)
})
