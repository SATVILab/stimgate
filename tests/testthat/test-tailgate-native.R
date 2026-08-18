.sourceLegacyTailgate <- function() {
  candidate_dirs <- c(
    file.path(testthat::test_path(), "..", ".."),
    file.path(testthat::test_path(), ".."),
    getwd(),
    file.path(getwd(), "..")
  )

  found <- FALSE
  for (dir in candidate_dirs) {
    p1 <- file.path(dir, "scripts", "r", "openCyto-find_peaks_and_valleys.R")
    p2 <- file.path(dir, "scripts", "r", "cytoUtils-cytokine_cutpoint.R")
    if (file.exists(p1) && file.exists(p2)) {
      source(p1, local = parent.frame())
      source(p2, local = parent.frame())
      found <- TRUE
      break
    }
  }

  if (!found) {
    .loadLegacyTailgateFallback(env = parent.frame())
  }
}

test_that("native tailgate preserves the relative-derivative shoulder rule", {
  .sourceLegacyTailgate()
  x <- seq(0, 8, length.out = 500)
  y <- stats::dnorm(x, mean = 2.5, sd = 0.8)

  tail <- stimgate:::.getStimGateTailgate(
    density = data.frame(x = x, y = y),
    peakX = 2.5,
    fraction = 1 / 200
  )

  expect_true(is.finite(tail$lowerBoundX))
  expect_gt(tail$lowerBoundX, 2.5)
  expect_lt(tail$lowerBoundX, max(x))
  expect_equal(tail$info$reason, "identified_stimulated_density_shoulder_lower_bound")
})

test_that("native tailgate returns NA when the descending shoulder never reaches the target", {
  x <- seq(0, 1, length.out = 300)
  y <- seq(0.1, 1, length.out = 300)

  tail <- stimgate:::.getStimGateTailgate(
    density = data.frame(x = x, y = y),
    peakX = max(x),
    fraction = 1 / 200
  )

  expect_false(is.finite(tail$lowerBoundX))
  expect_equal(tail$info$reason, "no_negative_derivative_right_of_stimulated_peak")
})

test_that("native tailgate has a deterministic multimodal fixture and differs from legacy tol gating", {
  .sourceLegacyTailgate()
  set.seed(1)
  x <- c(rnorm(200, mean = 0, sd = 0.6), rnorm(200, mean = 4, sd = 0.8))
  density <- stats::density(x, n = 512)

  native <- stimgate:::.getStimGateTailgate(
    density = data.frame(x = density$x, y = density$y),
    fraction = 1 / 200
  )
  wrapper <- stimgate:::.getCpUnsLocMarginalDensityLowerBound(
    density = data.frame(x = density$x, y = density$y),
    peakX = density$x[which.max(density$y)],
    fraction = 1 / 200
  )
  legacy <- .cytokineCutpoint(
    x = x,
    numPeaks = 1,
    refPeak = 1,
    tol = 1e-2,
    side = "right",
    strict = FALSE,
    plot = FALSE
  )

  expect_true(is.finite(native$lowerBoundX))
  expect_true(is.finite(legacy))
  expect_equal(wrapper$lowerBoundX, native$lowerBoundX, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(native$lowerBoundX, legacy, tolerance = 1e-8)))
  expect_true(native$lowerBoundX > min(x))
  expect_true(native$lowerBoundX < max(x))
})

test_that("legacy .getCpTg delegates to the native tailgate helper", {
  exList <- list(
    uns = data.frame(marker = c(1, 2, 3, 4, 5)),
    stim = data.frame(marker = c(1, 3, 5, 7, 9))
  )
  attr(exList[[1]], "chnlCut") <- "marker"
  attr(exList[[2]], "chnlCut") <- "marker"
  attr(exList[[1]], "batch") <- "batch1"
  attr(exList[[2]], "batch") <- "batch1"

  chnlSettings <- list(
    chnlCut = "marker",
    gateCombn = "no",
    excMin = FALSE,
    minCell = 5,
    cpMin = 0,
    bw = 1
  )

  local_mocked_bindings(
    .getCpTailgate = function(density, peakX = NULL, fraction = 1 / 200) {
      expect_equal(fraction, 1 / 200)
      list(
        lowerBoundX = 6,
        info = list(
          reason = "identified_stimulated_density_shoulder_lower_bound"
        )
      )
    }
  )

  out <- stimgate:::.getCpTg(
    exList = exList,
    chnlSettings = chnlSettings,
    tgType = "tg",
    stage = "init",
    pathProject = tempdir()
  )

  expect_equal(out[["no"]][["stim"]], 6)
})
