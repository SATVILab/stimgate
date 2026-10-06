# Contract tests for the gateStim() tuning controls: the stimControl() object
# holds the tuning settings, gateStim() keeps only its data/marker arguments
# plus biasUns and the fixed bw, and markerControl holds per-channel overrides
# keyed by marker label or channel name.

test_that("stimControl returns a documented stimControl object", {
  control <- stimControl()
  expect_s3_class(control, "stimControl")
  expect_type(control, "list")
  expect_setequal(names(control), names(formals(stimControl)))

  # Documented defaults.
  expect_identical(control$biasUnsFactor, 1)
  expect_true(control$excMin)
  expect_null(control$cpMin)
  expect_true(control$calcCytPosGates)
  expect_identical(control$minCell, 1e2)
  expect_identical(control$bwMtd, "nrd0")
  expect_identical(control$bwNcellMax, 1e4)
  # Tubes below bwNcellMin are upsampled to it and tubes above bwNcellMax
  # downsampled to it. bwNcellMin defaults to bwNcellMax, so every tube is
  # resampled to exactly bwNcellMax cells; a lower bwNcellMin leaves tubes
  # between the two unresampled.
  expect_identical(control$bwNcellMin, 1e4)
  expect_identical(stimControl(bwNcellMax = 5000)$bwNcellMin, 5000)
  expect_identical(stimControl(bwNcellMin = 100)$bwNcellMin, 100)
  expect_identical(control$bwScope, "cytokine")
  expect_identical(control$bwAdj, 1)
  expect_identical(control$bwFallback, "auto")
  expect_true(control$clusterGates)
  expect_identical(control$gateCombn, "min")
  expect_identical(control$locProbCol, "pred")
  expect_identical(control$locMinPeakProb, 0.25)
  expect_false(control$locEnforceShapeThreshold)
  expect_false(control$bwAdaptive)
  expect_identical(control$normMtd, "moments")

  # Supplied values are stored and the class is retained.
  override <- stimControl(
    clusterGates = FALSE,
    bwMtd = "nrd0",
    locProbCol = "probSmooth"
  )
  expect_s3_class(override, "stimControl")
  expect_false(override$clusterGates)
  expect_identical(override$bwMtd, "nrd0")
  expect_identical(override$locProbCol, "probSmooth")
})

test_that("stimControl rejects invalid settings eagerly", {
  # No project, data or marker is needed to surface these errors.
  expect_error(stimControl(bwMtd = "bad"))
  expect_error(stimControl(bwAdj = 0))
  expect_error(stimControl(bwScope = "batch"))
  expect_error(stimControl(gateCombn = "foo"))
  expect_error(stimControl(locProbCol = "foo"))
  expect_error(stimControl(locMinPeakProb = 1.5))
  expect_error(stimControl(normDensityN = 0))
  expect_error(stimControl(excMin = NULL))
  expect_error(stimControl(bwMtd = NULL))
  expect_error(stimControl(calcCytPosGates = "yes"))
  expect_error(stimControl(minCell = 0))
  expect_error(stimControl(locEnforceShapeThreshold = "yes"))
  expect_no_error(stimControl(bwMtd = "sj", locProbCol = "probSmooth"))
  expect_no_error(stimControl(locEnforceShapeThreshold = TRUE))
})

test_that("bw is a gateStim argument, not a stimControl setting", {
  # `bw` is the user-facing fixed bandwidth and lives on gateStim().
  expect_false("bw" %in% names(formals(stimControl)))
  expect_true("bw" %in% names(formals(gateStim)))
  expect_error(stimControl(bw = 0.1))
})

test_that("clusterGates is a logical switch", {
  expect_true(stimControl()$clusterGates)
  expect_true(stimControl(clusterGates = TRUE)$clusterGates)
  expect_false(stimControl(clusterGates = FALSE)$clusterGates)

  # The old numeric tolerance is not a valid value.
  expect_error(stimControl(clusterGates = 1e-7))
  expect_error(stimControl(clusterGates = "yes"))
  expect_error(stimControl(clusterGates = NULL))
  expect_error(stimControl(clusterGates = c(TRUE, FALSE)))
})

test_that("gateStim requires a stimControl control object", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(withr::local_tempdir(), "project")

  base <- list(
    pathProject = pathProject,
    .data = gs,
    batchList = exampleData$batchList,
    marker = exampleData$marker
  )
  # The message names the control class and points at stimControl().
  expect_error(
    do.call(gateStim, c(base, list(control = list()))),
    "stimControl"
  )
  expect_error(
    do.call(gateStim, c(base, list(control = "hpi1"))),
    "stimControl"
  )

  # The error is raised before any project directory is created.
  expect_false(dir.exists(pathProject))
})

test_that("markerControl accepts marker and channel names", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(withr::local_tempdir(), "project")

  # One key is a marker label, the other the channel name of another channel.
  # `bw` is a per-marker fixed-bandwidth override.
  markerControl <- stats::setNames(
    list(list(bw = 0.1), list(bwAdj = 2)),
    c(exampleData$marker[[1]], exampleData$chnl[[2]])
  )
  gateStim(
    pathProject = pathProject,
    .data = gs,
    batchList = exampleData$batchList,
    marker = exampleData$marker,
    control = stimControl(calcCytPosGates = FALSE, clusterGates = FALSE),
    markerControl = markerControl
  )

  expect_true(file.exists(file.path(pathProject, "gateStats.rds")))
  expect_gt(nrow(getStimGates(pathProject)), 0L)

  # The per-marker fixed bandwidth reaches the saved channel settings.
  settings <- stimgateMetaReadSettingsChnls(pathProject)
  expect_identical(settings[[exampleData$marker[[1]]]]$bw, 0.1)
})

test_that("markerControl rejects unknown names and settings", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(withr::local_tempdir(), "project")
  markerLabel <- exampleData$marker[[1]]
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

  # Unknown marker or channel name.
  expect_error(gateWith(list(notAChannel = list(bwAdj = 2))), "notAChannel")
  # Unknown setting name.
  expect_error(
    gateWith(stats::setNames(list(list(notASetting = 1)), markerChnl))
  )
  # locEnforceShapeThreshold and calcCytPosGates are global only.
  expect_error(gateWith(stats::setNames(
    list(list(locEnforceShapeThreshold = TRUE)), markerChnl
  )))
  expect_error(gateWith(stats::setNames(
    list(list(calcCytPosGates = FALSE)), markerChnl
  )))
  # Two keys that resolve to the same channel.
  expect_error(gateWith(stats::setNames(
    list(list(bwAdj = 2), list(bwAdj = 2)),
    c(markerLabel, markerChnl)
  )))
  # Invalid override values are rejected like the global settings.
  expect_error(gateWith(stats::setNames(list(list(bwAdj = 0)), markerChnl)))

  expect_false(dir.exists(pathProject))
})
