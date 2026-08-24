root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
script_gate <- file.path(root_dir, "scripts", "r", "acs_cytof-gate.R")
script_methods <- file.path(root_dir, "scripts", "r", "acs_cytof-methods.R")
script_manual <- file.path(root_dir, "scripts", "r", "acs_cytof-manual.R")

.load_acs_method_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_gate, local = env)
  source(script_methods, local = env)
  source(script_manual, local = env)
  env
}

test_that("ACS comparator settings pin F-beta defaults and Tailgate auto tuning", {
  env <- .load_acs_method_env()

  fbeta <- env$.acsCytofComparatorSettings("fbeta")
  expect_equal(fbeta$params$beta, 0.8)
  expect_equal(fbeta$params$theta, 2)
  expect_equal(fbeta$params$width, 10L)
  expect_null(fbeta$params$numBins)

  tailgate <- env$.acsCytofComparatorSettings("tailgate")
  expect_identical(tailgate$params$tailgateX, "stim")
  expect_null(tailgate$params$bandwidth)
  expect_true(tailgate$params$autoTol)
  expect_identical(tailgate$params$derivativeMethod, "firstDeriv")
})

test_that("thresholded cells are saved as a complete combination table", {
  env <- .load_acs_method_env()
  channels <- c("A", "B")
  x_stim <- matrix(
    c(0, 0, 2, 0, 0, 2, 2, 2),
    ncol = 2L,
    byrow = TRUE
  )
  x_uns <- matrix(c(0, 0, 2, 2), ncol = 2L, byrow = TRUE)

  out <- env$.acsCytofCombinationCounts(
    xStim = x_stim,
    xUns = x_uns,
    thresholds = c(A = 1, B = 1),
    channels = channels
  )

  expect_equal(nrow(out), 4L)
  expect_equal(sum(out$countStim), nrow(x_stim))
  expect_equal(sum(out$countUns), nrow(x_uns))
  expect_equal(out$countStim, rep(1L, 4L))
  expect_equal(out$countUns, c(1L, 0L, 0L, 1L))
})

test_that("comparator cache validation includes settings and sample count", {
  env <- .load_acs_method_env()
  settings <- env$.acsCytofComparatorSettings("fbeta")
  cache <- list(
    settings = settings,
    nSample = 10L,
    stats = tibble::tibble(
      method = "fbeta",
      pop = "cd4",
      ind = "2",
      cytCombn = "A~+~",
      countStim = 1L,
      nCellStim = 2L,
      countUns = 0L,
      nCellUns = 2L
    ),
    thresholds = tibble::tibble(
      method = "fbeta",
      pop = "cd4",
      ind = "2",
      chnl = "A",
      threshold = 1,
      thresholdOrigin = "calculated",
      thresholdFallbackUsed = FALSE
    )
  )

  expect_true(env$.acsCytofCacheIsCurrent(cache, settings, nSample = 10L))
  expect_false(env$.acsCytofCacheIsCurrent(cache, settings, nSample = 20L))

  changed_settings <- settings
  changed_settings$params$beta <- 1
  expect_false(env$.acsCytofCacheIsCurrent(cache, changed_settings, nSample = 10L))
})

test_that("ACS comparison methods are deliberately sequential and cache each method", {
  env <- .load_acs_method_env()
  runner_body <- paste(
    deparse(body(env$.acsCytofRunComparisonMethods)),
    collapse = "\n"
  )

  expect_true(grepl("lapply(methodVec", runner_body, fixed = TRUE))
  expect_false(grepl("future_map", runner_body, fixed = TRUE))
  expect_true(grepl(".acsCytofWriteComparatorCache", runner_body, fixed = TRUE))
  expect_true(grepl("Using cached", runner_body, fixed = TRUE))
})

test_that("manual formatting uses the shared cytometry combination utilities", {
  env <- .load_acs_method_env()
  formatter_body <- paste(
    deparse(body(env$.acsCytofStatsSingleMarkers)),
    collapse = "\n"
  )
  converter_body <- paste(
    deparse(body(env$.acsCytofCombinationToStandard)),
    collapse = "\n"
  )

  expect_true(grepl(
    "UtilsCytoRSV::sum_over_markers",
    formatter_body,
    fixed = TRUE
  ))
  expect_true(grepl(
    "UtilsCompassSV::convert_cyt_combn_format",
    converter_body,
    fixed = TRUE
  ))
})

test_that("combination counts collapse to one positive row per cytokine", {
  skip_if_not_installed("UtilsCytoRSV")
  skip_if_not_installed("UtilsCompassSV")
  env <- .load_acs_method_env()
  channels <- names(env$.acsCytofChannelMap())
  combinations <- env$.acsCytofCombinationLevels(channels)
  stats_tbl <- tibble::tibble(
    gateName = "loc_min",
    ind = "2",
    cytCombn = combinations,
    countStim = 1L,
    nCellStim = length(combinations),
    countUns = 0L,
    nCellUns = length(combinations)
  )
  sample_map <- tibble::tibble(
    popCode = "cd4",
    ind = "2",
    SampleID = "000001_D0",
    stim = "mtb"
  )

  out <- env$.acsCytofStatsSingleMarkers(
    statsTbl = stats_tbl,
    method = "stimgate",
    popCode = "cd4",
    popLabel = "CD4 T cells",
    sampleMap = sample_map
  )

  expect_equal(nrow(out), length(channels))
  expect_equal(out$countStim, rep(32L, length(channels)))
  expect_equal(as.character(out$cyt), unname(env$.acsCytofChannelMap()))
  expect_equal(out$cytCombn, paste0(out$cyt, "+"))
})
