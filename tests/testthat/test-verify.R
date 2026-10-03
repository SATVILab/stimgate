# Behavioural safety net for R/verify.R: each case pins only whether a value
# errors, not the error wording.

test_that("verifyGlobalAndPerChannelAgreeOnSharedSettings", {
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
  chnlCall <- function(vals) {
    .verifyChnlSettings(
      chnlSettings = list(ch = vals),
      chnl = "ch",
      markerSettings = NULL,
      marker = NULL
    )
  }
  markerCall <- function(vals) {
    .verifyChnlSettings(
      chnlSettings = NULL,
      chnl = NULL,
      markerSettings = list(mk = vals),
      marker = "mk"
    )
  }

  passes <- function(f, vals) {
    tryCatch(
      {
        f(vals)
        "ok"
      },
      error = conditionMessage
    )
  }

  # Default arguments pass both validators.
  expect_no_error(globalCall(list()))
  expect_no_error(chnlCall(list()))

  # Each case: list(setting, invalid value, valid value).
  cases <- list(
    list("excMin", "yes", FALSE),
    list("biasUns", "a", 0.5),
    list("biasUnsFactor", 0, 2),
    list("cpMin", "a", 0.1),
    list("maxPosProbX", "a", 0.9),
    list("bw", -1, 0.5),
    list("bwAdj", 0, 2),
    list("bwCluster", -1, 0.5),
    list("tolClust", 0, 1e-5),
    list("bwMtd", "foo", "sj"),
    list("gateCombn", "foo", c("min", "max")),
    list("gateQuant", c(-0.1, 0.5), c(0.1, 0.9)),
    list("gateQuant", 0.5, c(0.1, 0.9)),
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
    expect_error(globalCall(bad), info = paste("global invalid", nm))
    expect_error(chnlCall(bad), info = paste("chnl invalid", nm))
    expect_identical(passes(globalCall, good), "ok", info = nm)
    expect_identical(passes(chnlCall, good), "ok", info = nm)
  }

  # Cross-setting rules shared by both paths.
  combos <- list(
    list(bwMin = 1, bwMax = 0.5),
    list(bwAdaptive = TRUE, normMtd = "boxcox")
  )
  for (vals in combos) {
    expect_error(globalCall(vals), info = names(vals)[[1]])
    expect_error(chnlCall(vals), info = names(vals)[[1]])
  }
  expect_no_error(globalCall(list(bwMin = 0.1, bwMax = 0.5)))
  expect_no_error(chnlCall(list(bwMin = 0.1, bwMax = 0.5)))

  # markerSettings routes through the same per-channel checks.
  expect_error(markerCall(list(bwAdj = 0)))
  expect_no_error(markerCall(list(bwAdj = 2)))
})

test_that("verifyGateInputsRejectsInvalidGlobalOnlyArguments", {
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
    list(calcCytPosGates = "yes"),
    list(minCell = 0),
    list(bwNcellMin = "a"),
    list(bwNcellMax = "a"),
    list(bwNcellMin = 100, bwNcellMax = 10),
    list(marker = "notAMarker"),
    list(marker = NULL, chnl = "notAChannel"),
    list(chnl = exampleData$chnl),
    list(marker = NULL),
    list(marker = NULL, chnl = exampleData$chnl, markerSettings = list()),
    list(chnlSettings = list()),
    list(popGate = NULL),
    list(excMin = NULL),
    list(biasUnsFactor = NULL),
    list(maxPosProbX = NULL),
    list(bwAdj = NULL),
    list(bwMtd = NULL),
    list(gateCombn = NULL),
    list(gateQuant = NULL)
  )
  for (vals in invalid) {
    expect_error(globalCall(vals), info = paste(names(vals), collapse = ","))
  }

  expect_no_error(globalCall(list(minCell = 10, bwNcellMax = 1e4)))
  # Globally, `bwMtd` is unused and so unchecked when `bw` is fixed.
  expect_no_error(globalCall(list(bw = 0.5, bwMtd = NULL)))
  expect_no_error(globalCall(list(bw = 0.5, bwMtd = "foo")))
  expect_no_error(globalCall(list(marker = NULL, chnl = exampleData$chnl)))
})

test_that("verifyChnlSettingsRejectsInvalidStructure", {
  chnlCall <- function(chnlSettings, chnl = "ch") {
    .verifyChnlSettings(
      chnlSettings = chnlSettings,
      chnl = chnl,
      markerSettings = NULL,
      marker = NULL
    )
  }

  expect_error(chnlCall(list(ch = list(notASetting = 1))))
  # Legacy setting stays global only.
  expect_error(chnlCall(list(ch = list(locEnforceShapeThreshold = TRUE))))
  expect_error(
    .verifyChnlSettings(
      chnlSettings = NULL,
      chnl = NULL,
      markerSettings = list(mk = list(notASetting = 1)),
      marker = "mk"
    )
  )
  expect_error(chnlCall("notAList"))
  expect_error(chnlCall(list(ch = "notAList")))
  expect_error(chnlCall(list(list(bwAdj = 2))))
  expect_error(chnlCall(list(ch = list(), ch = list())))
  expect_error(chnlCall(list(other = list())))
  expect_error(chnlCall(list(ch = list(popGate = 1))))

  expect_no_error(chnlCall(NULL))
  expect_no_error(chnlCall(list(ch = list(popGate = "root", bwAdj = 2))))
})

test_that("gateStimFailsFastOnInvalidSettings", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(withr::local_tempdir(), "project")

  expect_error(
    gateStim(
      pathProject = pathProject,
      .data = gs,
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      biasUnsFactor = 0
    ),
    "biasUnsFactor"
  )
  expect_false(dir.exists(pathProject))

  chnlSettings <- stats::setNames(
    list(list(bwAdj = 0)),
    exampleData$chnl[[1]]
  )
  expect_error(
    gateStim(
      pathProject = pathProject,
      .data = gs,
      batchList = exampleData$batchList,
      chnl = exampleData$chnl,
      chnlSettings = chnlSettings
    ),
    "bwAdj"
  )
})
