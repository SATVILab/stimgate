# Optional shifted-peak rule (`stimControl(locShiftedPeakRef = TRUE)`): when
# most stimulated cells respond, the stimulated main peak is the responders,
# so the search for responding cells starts above the unstimulated peak.

.shiftedPeakSim <- function(seed, fracResp, respMean = 3, n = 5000L) {
  withr::with_seed(seed, {
    mk <- function(f) {
      nResp <- round(n * f)
      matrix(
        c(stats::rnorm(n - nResp, 0, 0.3), stats::rnorm(nResp, respMean, 0.4)),
        ncol = 1L,
        dimnames = list(NULL, "cyt")
      )
    }
    list(
      uns1 = mk(0.002), stim1 = mk(fracResp),
      uns2 = mk(0.002), stim2 = mk(fracResp)
    )
  })
}

.shiftedPeakGate <- function(data, ...) {
  pathProject <- withr::local_tempdir("shifted_peak")
  ctrl <- stimControl(clusterGates = FALSE, gateCombn = "no", ...)
  withr::with_seed(1L, suppressMessages(gateStim(
    .data = data, pathProject = pathProject,
    batchList = list(c(1, 2), c(3, 4)), chnl = "cyt", control = ctrl
  )))
  gates <- getStimGates(pathProject)
  gates[gates$gateName == "loc_no", , drop = FALSE]
}

.shiftedPeakFreq <- function(data, gate) {
  mean(data$stim1 > gate) - mean(data$uns1 > gate)
}

test_that("stimControl() validates the shifted-peak settings", {
  expect_false(stimControl()$locShiftedPeakRef)
  expect_identical(stimControl()$locShiftedPeakBwMult, 2)
  expect_error(stimControl(locShiftedPeakRef = NA), "locShiftedPeakRef")
  expect_error(stimControl(locShiftedPeakRef = "yes"), "locShiftedPeakRef")
  expect_error(stimControl(locShiftedPeakBwMult = 0), "locShiftedPeakBwMult")
  expect_error(
    stimControl(locShiftedPeakBwMult = c(1, 2)), "locShiftedPeakBwMult"
  )
})

test_that("the rule fires only beyond the bandwidth multiple", {
  settings <- list(mult = 2, bw = 0.1, bwSource = "shared")
  fired <- stimgate:::.getCpUnsLocShiftedPeakRule(settings, 0.25, 0, 0.5)
  notFired <- stimgate:::.getCpUnsLocShiftedPeakRule(settings, 0.15, 0, 0.5)
  expect_true(fired$applied)
  expect_false(notFired$applied)
  expect_null(stimgate:::.getCpUnsLocShiftedPeakRule(NULL, 3, 0, 0.5))
  # Without an unstimulated left window there is nothing to start from.
  expect_false(
    stimgate:::.getCpUnsLocShiftedPeakRule(settings, 3, 0, NA_real_)$applied
  )
  # Adaptive bandwidths are read at the unstimulated peak.
  adaptive <- list(
    mult = 2, bw = NA_real_, bwSource = "adaptive",
    bwGrid = list(x = c(-1, 1), y = c(0.1, 0.3))
  )
  rule <- stimgate:::.getCpUnsLocShiftedPeakRule(adaptive, 0.5, 0, 0.5)
  expect_equal(rule$bw, 0.2)
  expect_true(rule$applied)
})

test_that("the reference bandwidth is the unscaled shared bandwidth", {
  exStim <- data.frame(val = 1:3)
  attr(exStim, "ind") <- 2L
  exUns <- data.frame(val = 1:3)
  attr(exUns, "ind") <- 1L
  get <- function(chnlSettings, densityBw = 0.5) {
    stimgate:::.getCpUnsLocShiftedPeakSettings(
      c(chnlSettings, locShiftedPeakRef = TRUE), exStim, exUns, densityBw
    )
  }
  expect_null(stimgate:::.getCpUnsLocShiftedPeakSettings(
    list(), exStim, exUns, 0.5
  ))
  shared <- get(list(bwShared = 0.1, sampleScale = 1.5))
  expect_identical(shared$bw, 0.1)
  expect_identical(shared$bwSource, "shared")
  fixed <- get(list(bw = 0.3, bwShared = 0.1))
  expect_identical(fixed$bw, 0.3)
  expect_identical(fixed$bwSource, "fixed")
  sample <- get(list(), densityBw = 0.4)
  expect_identical(sample$bw, 0.4)
  expect_identical(sample$bwSource, "sample")
  adaptive <- get(list(), densityBw = list(
    adaptive = TRUE, grid = 1:2, sharedGrid = c(0.1, 0.2)
  ))
  expect_identical(adaptive$bwSource, "adaptive")
  expect_identical(adaptive$bwGrid$y, c(0.1, 0.2))
})

test_that("gates are unchanged when the rule is off", {
  data <- .shiftedPeakSim(11L, fracResp = 0.8)
  default <- .shiftedPeakGate(data)
  off <- .shiftedPeakGate(data, locShiftedPeakRef = FALSE)
  expect_identical(default, off)
  expect_false("locShiftedPeakRef" %in% names(default))
})

test_that("most stimulated cells shifted right gives a much lower gate", {
  data <- .shiftedPeakSim(11L, fracResp = 0.8)
  default <- .shiftedPeakGate(data)
  shifted <- .shiftedPeakGate(data, locShiftedPeakRef = TRUE)
  expect_true(all(shifted$locShiftedPeakRef))
  expect_true(all(shifted$gate < default$gate - 1))
  # Default misses most responders (about 9% of 80% since the search start is
  # measured from the main negative peak's half-width); the rule recovers
  # about 80%.
  expect_lt(.shiftedPeakFreq(data, default$gate[[1]]), 0.15)
  expect_equal(.shiftedPeakFreq(data, shifted$gate[[1]]), 0.8, tolerance = 0.05)
})

test_that("a stimulated peak within two bandwidths leaves gates unchanged", {
  data <- .shiftedPeakSim(11L, fracResp = 0.05)
  default <- .shiftedPeakGate(data)
  shifted <- .shiftedPeakGate(data, locShiftedPeakRef = TRUE)
  expect_false(any(shifted$locShiftedPeakRef))
  expect_identical(
    shifted[, setdiff(names(shifted), "locShiftedPeakRef")],
    default
  )
})
