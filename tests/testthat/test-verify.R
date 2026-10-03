# Behavioural safety net for the tuning-settings validation in R/control.R and
# the gateStim-level validation in R/verify.R: each case pins only whether a
# value errors, not the error wording.

test_that("stimControl validates every tuning setting eagerly", {
  expect_no_error(stimControl())

  # Each case: list(setting, invalid value, valid value).
  cases <- list(
    list("excMin", "yes", FALSE),
    list("biasUnsFactor", 0, 2),
    list("cpMin", "a", 0.1),
    list("calcCytPosGates", "yes", FALSE),
    list("minCell", 0, 1e2),
    list("bwAdj", 0, 2),
    list("bwCluster", -1, 0.5),
    list("bwScope", "batch", "cluster"),
    list("clusterGates", "yes", TRUE),
    list("bwMtd", "foo", "sj"),
    list("gateCombn", "foo", c("min", "max")),
    list("bwMin", "foo", 0.1),
    list("bwMin", "foo", "none"),
    list("bwMin", "foo", -Inf),
    list("bwMin", "foo", -1),
    list("bwMax", -Inf, Inf),
    list("bwMax", 0, 0.5),
    list("bwMax", "foo", "none"),
    list("bwFallback", "none", 0.1),
    list("bwAdaptive", "yes", TRUE),
    list("bwAdaptiveDensityN", 0, 256),
    list("bwAdaptivePadFrac", -1, 0),
    list("bwAdaptiveCore", 0, 1),
    list("bwAdaptiveExtra", 0, 1),
    list("bwAdaptiveCrossover", Inf, 0.5),
    list("bwAdaptiveTransitionWidth", -1, 0),
    list("normPeakMinRel", -1, 0.5),
    list("normExtraFrac", 2, 0.1),
    list("normExtraMax", 0, 10),
    list("normLambda", c(1, Inf), c(-1, 1)),
    list("normDensityN", 0, 256),
    list("normExcessBwMtd", "nrd0Norm", "sj"),
    list("normExcessNcell", -1, 100),
    list("normAdaptiveNcell", 0, 100),
    list("normMtd", "foo", "boxcox"),
    list("locProbCol", "foo", "probSmooth"),
    list("locMinPeakProb", 1.5, 0.3),
    list("locEnforceShapeThreshold", "yes", TRUE),
    list("locDipAlpha", -0.1, 0.3),
    list("locAntimodeHeightFrac", 1.5, 0.3),
    list("locAntimodeLowRel", 1.5, 0.3),
    list("locAntimodeLowAbs", 1.5, 0.3),
    list("locFlatDerivFrac", 1.5, 0.3),
    list("locFlatHardDerivFrac", 1.5, 0.3),
    list("locMarginalPurityRel", 1.5, 0.3),
    list("locMarginalRefQuantile", "a", 0.3),
    list("locMarginalCellBinRatio", 0, 3)
  )
  for (case in cases) {
    nm <- case[[1]]
    bad <- stats::setNames(list(case[[2]]), nm)
    good <- stats::setNames(list(case[[3]]), nm)
    expect_error(do.call(stimControl, bad), info = paste("invalid", nm))
    expect_no_error(do.call(stimControl, good))
  }

  # Cross-setting rules.
  expect_error(stimControl(bwMin = 1, bwMax = 0.5))
  expect_error(stimControl(bwAdaptive = TRUE, normMtd = "boxcox"))
  expect_no_error(stimControl(bwMin = 0.1, bwMax = 0.5))

  # Settings that must always be supplied.
  for (nm in c("excMin", "biasUnsFactor", "bwAdj", "bwMtd", "gateCombn")) {
    expect_error(
      do.call(stimControl, stats::setNames(list(NULL), nm)),
      info = paste("missing", nm)
    )
  }
})

test_that("verifyGateInputs rejects invalid gateStim-level arguments", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- withr::local_tempdir()

  nmsFill <- setdiff(
    names(formals(.verifyGateInputs)),
    c("pathProject", ".data", "batchList", "marker")
  )
  base <- c(
    lapply(formals(gateStim)[nmsFill], eval),
    list(
      pathProject = pathProject,
      .data = gs,
      batchList = exampleData$batchList,
      marker = exampleData$marker
    )
  )
  globalCall <- function(vals) {
    args <- base
    args[names(vals)] <- vals
    do.call(.verifyGateInputs, args)
  }

  invalid <- list(
    list(pathProject = ""),
    list(pathProject = c("a", "b")),
    list(popGate = c("root", "x")),
    list(.data = data.frame(x = 1)),
    list(batchList = list()),
    list(bw = -1),
    list(marker = "notAMarker"),
    list(marker = NULL, chnl = "notAChannel"),
    list(chnl = exampleData$chnl),
    list(marker = NULL),
    list(popGate = NULL)
  )
  for (vals in invalid) {
    expect_error(globalCall(vals), info = paste(names(vals), collapse = ","))
  }

  expect_no_error(globalCall(list(bw = 0.5)))
  expect_no_error(globalCall(list(marker = NULL, chnl = exampleData$chnl)))
})

test_that("markerControl rejects malformed structures", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(withr::local_tempdir(), "project")
  markerChnl <- exampleData$chnl[[1]]

  gateWith <- function(markerControl) {
    gateStim(
      pathProject = pathProject,
      .data = gs,
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      markerControl = markerControl
    )
  }

  expect_error(gateWith("notAList"))
  expect_error(gateWith(list(list(bwAdj = 2))))
  expect_error(gateWith(stats::setNames(list("notAList"), markerChnl)))
  expect_error(gateWith(stats::setNames(list(list(popGate = 1)), markerChnl)))
  expect_error(gateWith(stats::setNames(list(list(biasUns = "a")), markerChnl)))

  # Validation happens before any project directory is created.
  expect_false(dir.exists(pathProject))
})
