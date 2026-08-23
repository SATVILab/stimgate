test_that("gateStim with gateCombn = 'min' applies minimum generated threshold across multi-condition batch", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProjectMin <- file.path(tempdir(), "test_gate_combn_min")

  on.exit(
    {
      if (dir.exists(pathProjectMin)) unlink(pathProjectMin, recursive = TRUE)
      if (dir.exists(exampleData$pathGs)) unlink(exampleData$pathGs, recursive = TRUE)
    },
    add = TRUE
  )

  # Multi-condition batch: sample 1 is unstim, samples 2 and 4 are two stimulated conditions
  batchListMulti <- list(batch1 = c(1L, 2L, 4L))

  # Run with gateCombn = "min" on the multi-condition batch
  invisible(gateStim(
    .data = gs,
    pathProject = pathProjectMin,
    popGate = "root",
    batchList = batchListMulti,
    marker = exampleData$marker,
    gateCombn = "min",
    tolClust = NULL,
    calcCytPosGates = FALSE
  ))

  gatesMin <- getStimGates(pathProjectMin)
  expect_true(is.data.frame(gatesMin) && nrow(gatesMin) > 0L)

  statsMin <- getStimStats(pathProjectMin)
  expect_true(is.data.frame(statsMin) && nrow(statsMin) > 0L)

  # Verify across markers that gateCombn = "min" applies identical combined threshold across the batch
  for (m in unique(as.character(gatesMin$marker))) {
    mGatesMin <- gatesMin[as.character(gatesMin$marker) == m, , drop = FALSE]
    stimGatesMin <- mGatesMin[as.character(mGatesMin$ind) %in% c("2", "4"), , drop = FALSE]

    # Both conditions in the batch receive identical combined gate
    expect_equal(unname(stimGatesMin$gate[1]), unname(stimGatesMin$gate[2]))
    expect_true(is.finite(stimGatesMin$gate[1]))
  }

  # Verify that statsMin has valid frequencies computed from the combined gate
  expect_true(all(!is.na(statsMin$propStim)))
  expect_true(all(!is.na(statsMin$propBs)))
  expect_true(all(!is.na(statsMin$freqBs)))
})

test_that("gateStim with gateCombn = 'min' combines generated thresholds on distinct responses", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  set.seed(123)
  nCells <- 1500L

  # Sample 1: Unstimulated control (negative peak at 1.0)
  m1_uns <- stats::rnorm(nCells, mean = 1.0, sd = 0.25)

  # Sample 2: Strong response (75% negative at 1.0, 25% positive at 4.5)
  m1_stim_strong <- c(
    stats::rnorm(nCells * 0.75, mean = 1.0, sd = 0.25),
    stats::rnorm(nCells * 0.25, mean = 4.5, sd = 0.3)
  )

  # Sample 3: Moderate response (80% negative at 1.0, 20% positive at 3.2)
  m1_stim_mod <- c(
    stats::rnorm(nCells * 0.8, mean = 1.0, sd = 0.25),
    stats::rnorm(nCells * 0.2, mean = 3.2, sd = 0.3)
  )

  mat1 <- matrix(m1_uns, ncol = 1, dimnames = list(NULL, "MarkerF1"))
  mat2 <- matrix(m1_stim_strong, ncol = 1, dimnames = list(NULL, "MarkerF1"))
  mat3 <- matrix(m1_stim_mod, ncol = 1, dimnames = list(NULL, "MarkerF1"))

  ff1 <- flowCore::flowFrame(mat1)
  ff2 <- flowCore::flowFrame(mat2)
  ff3 <- flowCore::flowFrame(mat3)

  fs <- flowCore::flowSet(list(ff1, ff2, ff3))
  flowCore::sampleNames(fs) <- c("sample1", "sample2", "sample3")
  gs <- flowWorkspace::GatingSet(fs)

  pathProjectNo <- file.path(tempdir(), "test_distinct_no")
  pathProjectMin <- file.path(tempdir(), "test_distinct_min")
  on.exit(
    {
      if (dir.exists(pathProjectNo)) unlink(pathProjectNo, recursive = TRUE)
      if (dir.exists(pathProjectMin)) unlink(pathProjectMin, recursive = TRUE)
    },
    add = TRUE
  )

  batchList <- list(batch1 = c(1L, 2L, 3L))

  # Gating with gateCombn = "no"
  invisible(gateStim(
    .data = gs,
    pathProject = pathProjectNo,
    popGate = "root",
    batchList = batchList,
    marker = "MarkerF1",
    gateCombn = "no",
    tolClust = NULL,
    calcCytPosGates = FALSE
  ))
  gatesNo <- getStimGates(pathProjectNo)
  gateStrong <- gatesNo$gate[as.character(gatesNo$ind) == "2"]
  gateMod <- gatesNo$gate[as.character(gatesNo$ind) == "3"]
  expect_true(gateStrong > gateMod)

  # Gating with gateCombn = "min"
  invisible(gateStim(
    .data = gs,
    pathProject = pathProjectMin,
    popGate = "root",
    batchList = batchList,
    marker = "MarkerF1",
    gateCombn = "min",
    tolClust = NULL,
    calcCytPosGates = FALSE
  ))
  gatesMin <- getStimGates(pathProjectMin)
  expect_equal(
    unname(gatesMin$gate[as.character(gatesMin$ind) == "2"]),
    unname(gateMod),
    tolerance = 1e-6
  )
  expect_equal(
    unname(gatesMin$gate[as.character(gatesMin$ind) == "3"]),
    unname(gateMod),
    tolerance = 1e-6
  )
})

test_that("gateStim with gateCombn = 'min' ignores fallback non-generated cutpoints in multi-condition batch", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  set.seed(42)
  nCells <- 1200L

  # Sample 1: Unstimulated control (negative peak at 1.0)
  m1_uns <- stats::rnorm(nCells, mean = 1.0, sd = 0.25)

  # Sample 2: Responsive condition (80% negative at 1.0, 20% positive at 4.0)
  m1_stim_resp <- c(
    stats::rnorm(nCells * 0.8, mean = 1.0, sd = 0.25),
    stats::rnorm(nCells * 0.2, mean = 4.0, sd = 0.3)
  )

  # Sample 3: Non-responsive condition (100% negative at 1.0, identical distribution to unstim)
  m1_stim_nonresp <- stats::rnorm(nCells, mean = 1.0, sd = 0.25)

  mat1 <- matrix(m1_uns, ncol = 1, dimnames = list(NULL, "MarkerF1"))
  mat2 <- matrix(m1_stim_resp, ncol = 1, dimnames = list(NULL, "MarkerF1"))
  mat3 <- matrix(m1_stim_nonresp, ncol = 1, dimnames = list(NULL, "MarkerF1"))

  ff1 <- flowCore::flowFrame(mat1)
  ff2 <- flowCore::flowFrame(mat2)
  ff3 <- flowCore::flowFrame(mat3)

  fs <- flowCore::flowSet(list(ff1, ff2, ff3))
  flowCore::sampleNames(fs) <- c("sample1", "sample2", "sample3")
  gs <- flowWorkspace::GatingSet(fs)

  pathProject <- file.path(tempdir(), "test_gate_combn_fallback")
  on.exit(
    {
      if (dir.exists(pathProject)) unlink(pathProject, recursive = TRUE)
    },
    add = TRUE
  )

  batchList <- list(batch1 = c(1L, 2L, 3L))

  invisible(gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = batchList,
    marker = "MarkerF1",
    gateCombn = "min",
    tolClust = NULL,
    calcCytPosGates = FALSE
  ))

  gates <- getStimGates(pathProject)
  expect_true(is.data.frame(gates) && nrow(gates) > 0L)

  gatesStim <- gates[as.character(gates$ind) %in% c("2", "3"), , drop = FALSE]

  # Sample 2 generated a real threshold in [1.8, 3.5]
  gateSample2 <- gatesStim$gate[as.character(gatesStim$ind) == "2"]
  expect_true(gateSample2 > 1.8 && gateSample2 < 3.5)

  # Sample 3 received Sample 2's generated threshold via min combination (not a fallback cutpoint)
  expect_equal(
    unname(gatesStim$gate[as.character(gatesStim$ind) == "3"]),
    unname(gateSample2),
    tolerance = 1e-6
  )

  # Stats reflect the combined gate
  stats <- getStimStats(pathProject)
  expect_true(is.data.frame(stats) && nrow(stats) > 0L)
  expect_true(all(!is.na(stats$propStim)))
  expect_true(all(!is.na(stats$propBs)))
})
