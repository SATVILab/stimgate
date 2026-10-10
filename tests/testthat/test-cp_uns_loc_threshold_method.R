# locThresholdMethod: "region" gates at the lower boundary of the filtered
# region (xSum); "match" keeps the probability-sum matching method; "cap" keeps
# the boundary unless its frequency exceeds the estimate by more than a factor.

pkg_ns <- asNamespace("stimgate")

# Ten stimulated and ten unstimulated cells. Filtering keeps the region from
# 2.2 upwards; fitted probabilities of 0.4 give a probability-sum estimate of
# 0.28, which matching reaches only at the cell at 3.8.
.locMethodFixture <- function() {
  ex <- function(x, ind) {
    structure(data.frame(IFNg = x), chnlCut = "IFNg", ind = ind)
  }
  dataMod <- structure(
    data.frame(
      IFNg = c(2.2, 2.3, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6),
      probSmooth = 0.4,
      pred = 0.4
    ),
    chnlCut = "IFNg"
  )
  list(
    dataMod = dataMod,
    stim = ex(c(0.5, 1, 1.5, 2.3, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6), 2L),
    uns = ex(seq(0.2, 2, by = 0.2), 1L)
  )
}

.locMethodGetCp <- function(fx, chnlSettings, xSum = 2.2, pathProject) {
  testthat::local_mocked_bindings(
    .getCpUnsLocFilterAfterSmoothing = function(dataMod, ...) {
      list(dataMod = dataMod, cp = NULL, info = list(final = list(xSum = xSum)))
    },
    .package = "stimgate"
  )
  pkg_ns$.getCpUnsLocGetCp(
    dataMod = fx$dataMod,
    exTblStimOrig = fx$stim,
    exTblStimNoMin = fx$stim,
    exTblUnsOrig = fx$uns,
    exTblUnsBias = fx$uns,
    bias = 0,
    cpMin = 0,
    stage = "init",
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )
}

test_that("stimControl defaults to the cap method and validates it", {
  expect_identical(stimControl()$locThresholdMethod, "cap")
  expect_identical(
    stimControl(locThresholdMethod = "match")$locThresholdMethod, "match"
  )
  expect_error(stimControl(locThresholdMethod = "sum"), "'region', 'match' or 'cap'")
  expect_error(stimControl(locThresholdMethod = c("region", "match")))
  expect_error(stimControl(locThresholdMethod = TRUE))
  expect_error(stimControl(locThresholdMethod = NULL), "locThresholdMethod")
  expect_error(stimControl(locThresholdMethod = NA), "locThresholdMethod")
  # Per-marker overrides use the same validation.
  expect_error(
    pkg_ns$.verifyChnlSettingsChnl("IFNg", list(locThresholdMethod = "x")),
    "'region', 'match' or 'cap'"
  )
})

test_that("missing per-marker and old control values resolve by inheritance", {
  resolve <- pkg_ns$.completeChnlSettingsLocThresholdMethod
  expect_identical(resolve("match", "region"), "match")
  expect_identical(resolve(NULL, "match"), "match")
  expect_identical(resolve(NA, "match"), "match")
  # A control object saved before the option existed has no value.
  expect_identical(resolve(NULL, NULL), "region")
  expect_identical(resolve(NA_character_, NA), "region")

  getMethod <- pkg_ns$.getCpUnsLocThresholdMethod
  expect_identical(getMethod(list()), "region")
  expect_identical(getMethod(list(locThresholdMethod = NA)), "region")
  expect_identical(getMethod(list(locThresholdMethod = "match")), "match")
})

test_that("region returns the boundary when matching moves the gate up", {
  fx <- .locMethodFixture()
  pathProject <- withr::local_tempdir()

  match <- .locMethodGetCp(
    fx, list(locThresholdMethod = "match"),
    pathProject = pathProject
  )
  # Matching selects the cell at 3.8 and places the gate halfway to the cell
  # below it (3.4), well above the region boundary.
  expect_equal(attr(match, "cpSelected"), 3.8)
  expect_equal(match$cp, 3.6)
  expect_identical(match$locReason, "local_fdr_threshold_selected")

  region <- .locMethodGetCp(fx, list(), pathProject = pathProject)
  expect_identical(region$cp, 2.2)
  expect_null(attr(region, "cpSelected"))
  expect_true(region$locGenerated)
  expect_true(region$locGeneratedDirect)
  expect_identical(region$locSource, "direct")
  expect_identical(region$locReason, "local_fdr_region_boundary_selected")
  expect_identical(attr(region, "locThresholdMethod"), "region")
  expect_identical(attr(region, "locRegionX"), 2.2)
})

test_that("explicit match uses the shared discounted estimate and selection", {
  fx <- .locMethodFixture()
  pathProject <- withr::local_tempdir()
  viaMethod <- .locMethodGetCp(
    fx, list(locThresholdMethod = "match"),
    pathProject = pathProject
  )
  # Verify routing through the shared estimator, not parity with legacy gates:
  # the bounded proportional discount intentionally changes matching too.
  dataThreshold <- pkg_ns$.getCpUnsLocGetCpDataThreshold(
    dataMod = fx$dataMod,
    exTblStimOrig = fx$stim,
    exTblStimNoMin = fx$stim,
    exTblUnsOrig = fx$uns,
    pathProject = pathProject,
    stage = "init"
  )
  shared <- pkg_ns$.getCpUnsLocGetCpActual(
    dataThreshold = dataThreshold,
    exTblStimNoMin = fx$stim,
    exTblUnsBias = fx$uns,
    cpMin = 0,
    stage = "init",
    exTblStimOrig = fx$stim,
    exTblUnsOrig = fx$uns,
    densityBw = NULL
  )
  for (nm in names(shared)) {
    expect_identical(viaMethod[[nm]], shared[[nm]])
  }
  expect_identical(attr(viaMethod, "cpSelected"), attr(shared, "cpSelected"))
})

test_that("region diagnostics count at the applied gate", {
  fx <- .locMethodFixture()
  withr::local_envvar(STIMGATE_INTERMEDIATE = "all")
  pathProject <- withr::local_tempdir()
  for (method in c("region", "match")) {
    cpObj <- .locMethodGetCp(
      fx, list(locThresholdMethod = method),
      pathProject = pathProject
    )
    detail <- readRDS(file.path(
      pathProject, "intermediateData", "init", "IFNg", "ind", "2",
      "locDetailCondition.rds"
    ))
    gate <- cpObj$cp
    expect_identical(detail$locThresholdMethod, method)
    expect_identical(detail$locRegionX, 2.2)
    expect_identical(detail$threshold, gate)
    expect_equal(detail$propStim, mean(fx$stim$IFNg > gate))
    expect_equal(detail$propUns, mean(fx$uns$IFNg > gate))
    expect_equal(detail$propBs, detail$propStim - detail$propUns)
    # The probability-sum estimate is reported separately, not forced to
    # agree with the region gate's frequency.
    expect_equal(detail$propBsEst, 0.28)
  }
  expect_equal(detail$propBsDiff, 0.3 - 0.28) # match: 3/10 cells above 3.6
  region <- .locMethodGetCp(fx, list(), pathProject = pathProject)
  detail <- readRDS(file.path(
    pathProject, "intermediateData", "init", "IFNg", "ind", "2",
    "locDetailCondition.rds"
  ))
  expect_equal(detail$propBs, 0.7)
  expect_equal(detail$propBsDiff, 0.7 - 0.28)
})

test_that("region keeps the no-response and failure fallbacks", {
  fx <- .locMethodFixture()
  pathProject <- withr::local_tempdir()

  # A non-finite boundary falls back rather than gating at NA.
  noBoundary <- .locMethodGetCp(
    fx, list(),
    xSum = NA_real_, pathProject = pathProject
  )
  expect_false(noBoundary$locGenerated)
  expect_false(noBoundary$locGeneratedDirect)
  expect_identical(noBoundary$locSource, "not_calculated")
  expect_identical(noBoundary$locReason, "No finite local-FDR region boundary")
  expect_true(is.finite(noBoundary$cp))
  expect_gt(noBoundary$cp, max(fx$stim$IFNg))

  # No responding cells: both methods give the same labelled fallback.
  empty <- structure(
    data.frame(IFNg = numeric(0), propBsDiff = numeric(0)),
    chnlCut = "IFNg"
  )
  args <- list(
    dataThreshold = empty, exTblStimNoMin = fx$stim, exTblUnsBias = fx$uns,
    cpMin = 0, stage = "init"
  )
  region <- do.call(pkg_ns$.getCpUnsLocGetCpRegion, c(args, regionX = 2.2))
  match <- do.call(pkg_ns$.getCpUnsLocGetCpActual, args)
  expect_identical(region, match)
  expect_false(region$locGenerated)
  expect_identical(region$locReason, "Too few responding cells")

  # Filtering that already returned a fallback cutpoint is left alone.
  testthat::local_mocked_bindings(
    .getCpUnsLocFilterAfterSmoothing = function(dataMod, ...) {
      list(
        dataMod = dataMod[0, , drop = FALSE], cp = 99,
        info = list(reason = "max_response_probability_below_minimum")
      )
    },
    .package = "stimgate"
  )
  filtered <- pkg_ns$.getCpUnsLocGetCp(
    dataMod = fx$dataMod, exTblStimOrig = fx$stim, exTblStimNoMin = fx$stim,
    exTblUnsOrig = fx$uns, exTblUnsBias = fx$uns, bias = 0, cpMin = 0,
    stage = "init", pathProject = pathProject, chnlSettings = list()
  )
  expect_identical(filtered$cp, 99)
  expect_false(filtered$locGenerated)
  expect_identical(filtered$locReason, "max_response_probability_below_minimum")
})

test_that("the shape-enforced route exposes its filtering boundary", {
  boundary <- pkg_ns$.getCpUnsLocShapeRegionBoundary
  dataMod <- structure(data.frame(IFNg = c(3, 4, 5)), chnlCut = "IFNg")
  # The kept region starts at the largest applied cut.
  out <- boundary(
    dataMod,
    shapeLowerBoundX = 1,
    globalInfo = list(applied = TRUE, thresholdX = 2),
    marginalInfo = list(finalStartX = 2.5)
  )
  expect_identical(out$xSum, 2.5)
  expect_identical(out$xSumSource, "marginal")
  # A global threshold that removed nothing is not a cut.
  out <- boundary(
    dataMod,
    shapeLowerBoundX = NA_real_,
    globalInfo = list(applied = FALSE, thresholdX = 2.9),
    marginalInfo = list()
  )
  expect_identical(out$xSum, 3)
  expect_identical(out$xSumSource, "lowest_kept_value")
})

test_that("gateStim applies the region gate on both filtering routes", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  withr::local_envvar(STIMGATE_INTERMEDIATE = "all")
  marker <- exampleData$marker[[1]]
  chnl <- exampleData$chnl[[1]]

  gateWith <- function(shape, method = NULL) {
    pathProject <- file.path(
      withr::local_tempdir(.local_envir = parent.frame(2)), "p"
    )
    control <- if (is.null(method)) {
      stimControl(
        clusterGates = FALSE, calcCytPosGates = FALSE,
        locEnforceShapeThreshold = shape
      )
    } else {
      stimControl(
        clusterGates = FALSE, calcCytPosGates = FALSE,
        locEnforceShapeThreshold = shape, locThresholdMethod = method
      )
    }
    withr::with_seed(1, gateStim(
      pathProject, gs, exampleData$batchList,
      marker = marker, control = control
    ))
  }

  for (shape in c(FALSE, TRUE)) {
    pathRegion <- gateWith(shape, "region")
    pathMatch <- gateWith(shape, "match")

    expect_identical(
      stimgateMetaReadSettingsChnls(pathRegion)[[marker]]$locThresholdMethod,
      "region"
    )
    expect_identical(
      stimgateMetaReadSettingsChnls(pathMatch)[[marker]]$locThresholdMethod,
      "match"
    )

    detail <- getStimGatesDetailed(pathRegion) |>
      dplyr::filter(.data$detailLevel == "condition", .data$stage == "init")
    direct <- detail$locGeneratedDirect
    expect_true(any(direct))
    expect_identical(unique(detail$locThresholdMethod), "region")
    expect_true(all(is.finite(detail$locRegionX[direct])))
    expect_identical(detail$threshold[direct], detail$locRegionX[direct])
    expect_true(all(
      detail$locReason[direct] == "local_fdr_region_boundary_selected"
    ))

    # Both methods share the filtering; matching may only raise the gate.
    detailMatch <- getStimGatesDetailed(pathMatch) |>
      dplyr::filter(.data$detailLevel == "condition", .data$stage == "init")
    detailMatch <- detailMatch[match(detail$ind, detailMatch$ind), ]
    expect_identical(detailMatch$locRegionX, detail$locRegionX)
    expect_true(all(detailMatch$threshold[direct] >= detail$threshold[direct]))

    # Public counts and frequencies describe the applied gate.
    gates <- getStimGates(pathRegion)
    stats <- getStimStats(pathRegion) |>
      dplyr::filter(grepl("~\\+~$", .data$cytCombn))
    for (i in seq_len(nrow(gates))) {
      ind <- gates$ind[[i]]
      # The unstimulated sample is first in the stimulated sample's batch.
      indUns <- Find(
        function(b) as.integer(ind) %in% b, exampleData$batchList
      )[[1]]
      expr <- function(j) {
        flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[j]], "root"))[, chnl]
      }
      row <- stats[stats$ind == ind, ]
      gate <- gates$gate[[i]]
      expect_identical(row$countStim, sum(expr(as.integer(ind)) > gate))
      expect_identical(row$countUns, sum(expr(indUns) > gate))
    }
  }
})

# The fixture's probability-sum estimate is 0.28 and the frequency above the
# region boundary (2.2) is 0.7. Candidate cells from 2.3 to 4.6 have
# frequencies (cells at or above them) of 0.7, 0.6, ..., 0.1.
test_that("cap keeps the region boundary when its frequency is within the cap", {
  fx <- .locMethodFixture()
  pathProject <- withr::local_tempdir()
  for (cap in c(3, Inf)) {
    cpObj <- .locMethodGetCp(
      fx, list(locThresholdMethod = "cap", locThresholdCap = cap),
      pathProject = pathProject
    )
    expect_identical(cpObj$cp, 2.2)
    expect_null(attr(cpObj, "cpSelected"))
    expect_identical(cpObj$locReason, "local_fdr_cap_region_boundary_selected")
    expect_identical(attr(cpObj, "locThresholdMethod"), "cap")
    expect_false(attr(cpObj, "locCapExceededAbove"))
  }
})

test_that("cap moves to the lowest candidate within the cap", {
  fx <- .locMethodFixture()
  pathProject <- withr::local_tempdir()
  # Limit 0.56: the first candidate at or below it is the cell at 3.0, and the
  # gate is placed halfway to the cell below it (2.6).
  cpObj <- .locMethodGetCp(
    fx, list(locThresholdMethod = "cap", locThresholdCap = 2),
    pathProject = pathProject
  )
  expect_equal(attr(cpObj, "cpSelected"), 3.0)
  expect_equal(cpObj$cp, 2.8)
  expect_identical(cpObj$locReason, "local_fdr_cap_threshold_selected")
  expect_true(cpObj$locGenerated)
  expect_false(attr(cpObj, "locCapExceededAbove"))

  # The default cap (1.3, limit 0.364) reaches the cell matching selects.
  cpDefault <- .locMethodGetCp(
    fx, list(locThresholdMethod = "cap"),
    pathProject = pathProject
  )
  match <- .locMethodGetCp(
    fx, list(locThresholdMethod = "match"),
    pathProject = pathProject
  )
  expect_equal(attr(cpDefault, "cpSelected"), attr(match, "cpSelected"))
  expect_equal(cpDefault$cp, match$cp)
})

test_that("cap falls back to matching when no candidate is within the cap", {
  fx <- .locMethodFixture()
  pathProject <- withr::local_tempdir()
  dataThreshold <- pkg_ns$.getCpUnsLocGetCpDataThreshold(
    dataMod = fx$dataMod,
    exTblStimOrig = fx$stim,
    exTblStimNoMin = fx$stim,
    exTblUnsOrig = fx$uns,
    pathProject = pathProject,
    stage = "init"
  )
  # An estimate of 0.05 gives a limit (0.065) below every candidate.
  dataThreshold$propBsDiff <- dataThreshold$propBs - 0.05
  cpObj <- pkg_ns$.getCpUnsLocGetCpCap(
    dataThreshold = dataThreshold, regionX = 2.2, cap = 1.3,
    exTblStimNoMin = fx$stim, exTblUnsBias = fx$uns, cpMin = 0,
    stage = "init", exTblStimOrig = fx$stim, exTblUnsOrig = fx$uns
  )
  expect_equal(attr(cpObj, "cpSelected"), 4.6)
  expect_equal(cpObj$cp, 4.4)
  expect_identical(cpObj$locReason, "local_fdr_cap_match_fallback")
})

test_that("cap diagnostics count at the applied gate", {
  fx <- .locMethodFixture()
  withr::local_envvar(STIMGATE_INTERMEDIATE = "all")
  pathProject <- withr::local_tempdir()
  for (cap in c(2, 3)) {
    cpObj <- .locMethodGetCp(
      fx, list(locThresholdMethod = "cap", locThresholdCap = cap),
      pathProject = pathProject
    )
    detail <- readRDS(file.path(
      pathProject, "intermediateData", "init", "IFNg", "ind", "2",
      "locDetailCondition.rds"
    ))
    gate <- cpObj$cp
    expect_identical(detail$locThresholdMethod, "cap")
    expect_identical(detail$threshold, gate)
    expect_equal(detail$propBs, mean(fx$stim$IFNg > gate) - mean(fx$uns$IFNg > gate))
    expect_equal(detail$propBsEst, 0.28)
    expect_false(detail$locCapExceededAbove)
  }
})

test_that("stimControl validates locThresholdCap", {
  expect_identical(stimControl()$locThresholdCap, 1.3)
  expect_identical(stimControl(locThresholdMethod = "cap")$locThresholdMethod, "cap")
  expect_no_error(stimControl(locThresholdCap = Inf))
  expect_error(stimControl(locThresholdCap = 0.9), "at least 1")
  expect_error(stimControl(locThresholdCap = "1.3"), "at least 1")
  expect_error(stimControl(locThresholdCap = c(1.2, 1.5)), "at least 1")
})

.locLeadingRunFixture <- function(x, pred, ratio, bw) {
  fx <- .locMethodFixture()
  fx$dataMod <- structure(
    data.frame(IFNg = x, pred = pred, probSmooth = pred * ratio),
    chnlCut = "IFNg", locDensityBw = bw
  )
  fx
}

.locLeadingRunThreshold <- function(fx) {
  pkg_ns$.getCpUnsLocGetCpDataThreshold(
    dataMod = fx$dataMod, exTblStimOrig = fx$stim,
    exTblStimNoMin = fx$stim, exTblUnsOrig = fx$uns,
    pathProject = withr::local_tempdir(), stage = "init",
    densityBw = attr(fx$dataMod, "locDensityBw")
  )
}

test_that("tiny disagreement over a long leading run retains candidates", {
  fx <- .locLeadingRunFixture(0:100, 0.99, 0.997, 2)
  fx$stim <- structure(data.frame(IFNg = 0:100), chnlCut = "IFNg", ind = 2L)
  out <- .locLeadingRunThreshold(fx)
  # The minimum is still excluded, but the old unlimited leading-run drop
  # would have removed every remaining row. Only x = 1 is discounted here.
  expect_equal(out$IFNg, 1:100)
  plain <- sum(fx$dataMod$pred[-1]) / nrow(fx$stim)
  estimate <- pkg_ns$.getCpUnsLocProbBsEst(out)
  expect_equal(estimate, (99 - 0.99 * 0.012) / nrow(fx$stim))
  expect_lt(abs(estimate - plain), 0.0002)
})

test_that("zero-weight leading rows are removed only within half a bandwidth", {
  fx <- .locLeadingRunFixture(0:6, 0.8, 0.5, 4)
  out <- .locLeadingRunThreshold(fx)
  expect_equal(out$IFNg, 3:6)
  expect_equal(out$weight, rep(1, 4))
  expect_equal(pkg_ns$.getCpUnsLocProbBsEst(out), 4 * 0.8 / nrow(fx$stim))

  # Adaptive shared bandwidth at x0 is 4, not the bandwidth farther right.
  attr(fx$dataMod, "locDensityBw") <- list(
    grid = c(0, 6), sharedGrid = c(4, 12)
  )
  expect_equal(.locLeadingRunThreshold(fx)$IFNg, 3:6)
})

test_that("partial leading weights are linear and end with the initial run", {
  fx <- .locLeadingRunFixture(0:4, 0.8, c(0.5, 0.875, 1, 0.5, 0.5), 20)
  out <- .locLeadingRunThreshold(fx)
  expect_equal(out$IFNg, 1:4)
  expect_equal(out$weight, c(0.5, 1, 1, 1))
  expect_equal(pkg_ns$.getCpUnsLocProbBsEst(out), (0.4 + 3 * 0.8) / nrow(fx$stim))
})

test_that("cap can select a gate inside the former leading run", {
  fx <- .locMethodFixture()
  fx$dataMod$probSmooth <- c(rep(0.4 * 0.997, 7), 0.4)
  attr(fx$dataMod, "locDensityBw") <- 1
  pathProject <- withr::local_tempdir()
  for (method in c("cap", "region", "match")) {
    out <- .locMethodGetCp(
      fx, list(locThresholdMethod = method, locThresholdCap = 2),
      pathProject = pathProject
    )
    # Every method uses the same discounted estimate, including match.
    expect_equal(out$propBsEst, (2 * 0.4 * 0.988 + 5 * 0.4) / 10)
    if (method == "cap") {
      # The old drop left only 4.6 as a candidate. The cap can now select 3.0.
      expect_equal(attr(out, "cpSelected"), 3.0)
      expect_equal(out$cp, 2.8)
      expect_identical(out$locReason, "local_fdr_cap_threshold_selected")
    }
  }
})

test_that("unavailable bandwidth keeps full weights and candidates", {
  bandwidths <- list(
    NULL, NA_real_, NaN, Inf, -Inf, 0, -1,
    list(grid = 0:1, sharedGrid = c(NA_real_, -1))
  )
  for (bw in bandwidths) {
    fx <- .locLeadingRunFixture(0:4, 0.8, 0.5, bw)
    out <- .locLeadingRunThreshold(fx)
    expect_equal(out$IFNg, 1:4)
    expect_equal(out$weight, rep(1, 4))
    expect_equal(pkg_ns$.getCpUnsLocProbBsEst(out), 4 * 0.8 / nrow(fx$stim))
  }
})

test_that("minimum exclusion retains its single-row exception", {
  count <- pkg_ns$.getCpUnsLocGetCpDataThresholdCount
  fx <- .locLeadingRunFixture(c(0, 0, 1), 0.8, 1, 2)
  expect_equal(count(fx$dataMod, 2)$IFNg, 1)
  fx <- .locLeadingRunFixture(0, 0.8, 0.997, 2)
  out <- count(fx$dataMod, 2)
  expect_equal(out$IFNg, 0)
  expect_equal(out$weight, 0.988)

  # Non-positive fits and non-finite ratios in the discount window count zero.
  fx <- .locLeadingRunFixture(0:3, c(0.8, 0, 0.8, 0.8), 1, 20)
  fx$dataMod$probSmooth <- c(0.4, -1, -Inf, 0.8)
  expect_equal(count(fx$dataMod, 20)$IFNg, 3)
})


test_that("the discount window starts after margin exclusion and sorts rows", {
  fx <- .locLeadingRunFixture(0:6, 0.8, 0.5, 2)
  attr(fx$dataMod, "minProbXPos") <- 2
  fx$dataMod <- fx$dataMod[c(7, 3, 1, 6, 4, 2, 5), , drop = FALSE]
  out <- .locLeadingRunThreshold(fx)
  # x0 is 2 after margin exclusion, so zero-weight candidates stop at 3.
  expect_equal(out$IFNg, 4:6)
  expect_equal(pkg_ns$.getCpUnsLocProbBsEst(out), 3 * 0.8 / nrow(fx$stim))
})
