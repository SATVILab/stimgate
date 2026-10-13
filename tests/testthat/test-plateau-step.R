test_that("stimControl() validates the plateau-step settings", {
  ctrl <- stimControl()
  expect_false(ctrl$locPlateauStep)
  expect_identical(ctrl$locPlateauMinCells, 2)
  expect_identical(ctrl$locPlateauGapBw, 1)
  expect_identical(ctrl$locPlateauAlignMax, 1.1)
  expect_identical(ctrl$locPlateauAlignUpper, 1.25)
  expect_error(stimControl(locPlateauStep = NA), "locPlateauStep")
  expect_error(stimControl(locPlateauMinCells = 0), "locPlateauMinCells")
  expect_error(stimControl(locPlateauGapBw = 0), "locPlateauGapBw")
  expect_silent(stimControl(locPlateauGapBw = Inf))
  expect_error(stimControl(locPlateauAlignMax = 0.9), "locPlateauAlignMax")
  expect_error(
    stimControl(locPlateauAlignMax = 1.3, locPlateauAlignUpper = 1.2),
    "locPlateauAlignUpper"
  )
})

test_that("the plateau step moves through pure spans and stops at gaps", {
  withr::local_seed(1)
  neg <- stats::rnorm(1e5)
  # responders spread from 3 to 6: the band just below a gate at 4.5 is
  # mostly responders, so the gate may move down
  xStim <- c(neg, stats::runif(2000, 3, 6))
  xUns <- stats::rnorm(1e5)
  step <- .getCpUnsLocPlateauStepGate(
    xStim = xStim, xUns = xUns, gate = 4.5, edge = 3, bw = 0.2,
    purityRel = 0.5
  )
  expect_lt(step$gate, 4.5)
  expect_gte(step$gate, 3)

  # one stray stimulated cell in the second step below the gate, after an
  # empty step: the gap rule stops the scan at the empty step, and without it
  # the minimum cell count still refuses a single cell
  xStim2 <- c(stats::rnorm(1e5), 5.65, stats::runif(200, 7, 9))
  xUns2 <- c(stats::rnorm(1e5), stats::runif(40, 7, 9))
  gapStop <- .getCpUnsLocPlateauStepGate(
    xStim = xStim2, xUns = xUns2, gate = 6, edge = 4.5, bw = 0.2,
    purityRel = 0.5, minCells = 1, gapBw = 1
  )
  expect_identical(gapStop$gate, 6)
  expect_match(gapStop$reason, "empty_gap")
  twoCells <- .getCpUnsLocPlateauStepGate(
    xStim = xStim2, xUns = xUns2, gate = 6, edge = 4.5, bw = 0.2,
    purityRel = 0.5, minCells = 2, gapBw = Inf
  )
  expect_identical(twoCells$gate, 6)
  noGuard <- .getCpUnsLocPlateauStepGate(
    xStim = xStim2, xUns = xUns2, gate = 6, edge = 4.5, bw = 0.2,
    purityRel = 0.5, minCells = 1, gapBw = Inf
  )
  expect_lt(noGuard$gate, 6)

  # no edge below the gate: unchanged
  none <- .getCpUnsLocPlateauStepGate(
    xStim = xStim, xUns = xUns, gate = 4.5, edge = NA_real_, bw = 0.2,
    purityRel = 0.5
  )
  expect_identical(none$gate, 4.5)
  expect_identical(none$reason, "no_plateau_edge")
})

test_that("tail alignment switches the plateau step off for mismatched tails", {
  mk <- function(nStim, nUns, n = 1e5) {
    data.frame(
      ind = as.character(seq_along(nStim)), nStim = nStim, nUns = nUns,
      nStimTotal = n, nUnsTotal = n
    )
  }
  aligned <- .getCpUnsLocPlateauAlign(
    mk(c(200, 210, 190, 205, 198), c(200, 200, 200, 200, 200)),
    alignMax = 1.1, alignUpper = 1.25
  )
  expect_true(aligned$group$pass)
  expect_true(all(aligned$tbl$allowed))

  # stimulated tails heavier throughout, as with a batch shift
  shifted <- .getCpUnsLocPlateauAlign(
    mk(c(300, 310, 290, 305, 298), c(200, 200, 200, 200, 200)),
    alignMax = 1.1, alignUpper = 1.25
  )
  expect_false(shifted$group$pass)
  expect_false(any(shifted$tbl$allowed))

  # an aligned group keeps the step, except for one clearly heavier tube
  oneHigh <- .getCpUnsLocPlateauAlign(
    mk(c(200, 205, 195, 200, 202, 198, 201, 199, 400), rep(200, 9)),
    alignMax = 1.1, alignUpper = 1.25
  )
  expect_identical(oneHigh$tbl$allowed[9], FALSE)

  # too few cells in the band: such tubes follow the group
  sparse <- .getCpUnsLocPlateauAlign(
    mk(c(200, 205, 195, 5), c(200, 200, 200, 4)),
    alignMax = 1.1, alignUpper = 1.25
  )
  expect_false(sparse$tbl$usable[4])
  expect_identical(sparse$tbl$allowed[4], sparse$group$pass)
})

test_that("gateStim() records the plateau step only when requested", {
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  runGate <- function(control) {
    pathProject <- tempfile("plateau")
    withr::local_seed(1)
    suppressMessages(gateStim(
      pathProject = pathProject,
      .data = gs,
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      calcCytPosGates = FALSE,
      control = control
    ))
    gates <- getStimGates(pathProject)
    unlink(pathProject, recursive = TRUE)
    gates
  }
  gatesDefault <- runGate(stimControl())
  gatesOff <- runGate(stimControl(locPlateauStep = FALSE))
  gatesOn <- runGate(stimControl(locPlateauStep = TRUE))
  expect_false("locPlateauStep" %in% names(gatesDefault))
  expect_identical(gatesDefault$gate, gatesOff$gate)
  expect_true("locPlateauStep" %in% names(gatesOn))
  expect_type(gatesOn$locPlateauStep, "logical")
})
