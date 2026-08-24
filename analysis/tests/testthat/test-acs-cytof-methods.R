root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
script_gate <- file.path(root_dir, "scripts", "r", "acs_cytof-gate.R")
script_methods <- file.path(root_dir, "scripts", "r", "acs_cytof-methods.R")
script_manual <- file.path(root_dir, "scripts", "r", "acs_cytof-manual.R")
script_compare <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")
script_fbeta <- file.path(root_dir, "scripts", "python", "fbeta.py")
qmd_path <- file.path(root_dir, "analysis", "9-real-compare-acs-cytof.qmd")

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
  expect_equal(fbeta$cacheVersion, 2L)

  tailgate <- env$.acsCytofComparatorSettings("tailgate")
  expect_identical(tailgate$params$tailgateX, "stim")
  expect_null(tailgate$params$bandwidth)
  expect_true(tailgate$params$autoTol)
  expect_identical(tailgate$params$derivativeMethod, "firstDeriv")
  expect_equal(tailgate$cacheVersion, 1L)
})

test_that("F-beta histogram support includes the stimulated response tail", {
  skip_if_not_installed("reticulate")

  env <- new.env(parent = getNamespace("stimgate"))
  source(script_compare, local = env)

  x_uns <- seq(0, 1, length.out = 400L)
  x_stim <- c(seq(0, 1, length.out = 350L), rep(4, 50L))
  out <- env$.simCompareFbetaThreshold(
    xUns = x_uns,
    xStim = x_stim,
    pathFbeta = script_fbeta,
    width = 2L,
    numBins = 20L
  )

  expect_gt(max(as.numeric(out$fbeta$pdfx)), max(x_uns))
  expect_gt(
    max(as.numeric(out$fbeta$pdfx)),
    max(x_stim) - diff(range(c(x_uns, x_stim))) / 20
  )
})

test_that("F-beta does not retain Python objects in a global R cache", {
  skip_if_not_installed("reticulate")

  env <- new.env(parent = getNamespace("stimgate"))
  source(script_compare, local = env)

  expect_false(exists(
    ".simCompareCacheEnv",
    envir = env,
    inherits = FALSE
  ))

  fbeta_env <- env$.simCompareFbetaEnvironment(
    pathFbeta = script_fbeta
  )
  out <- env$.simCompareFbetaThreshold(
    xUns = seq(0, 1, length.out = 400L),
    xStim = c(seq(0, 1, length.out = 350L), rep(4, 50L)),
    pathFbeta = script_fbeta,
    fbetaEnv = fbeta_env,
    width = 2L,
    numBins = 20L
  )

  expect_true(is.finite(out$threshold))
  expect_false(exists(
    ".simCompareCacheEnv",
    envir = env,
    inherits = FALSE
  ))
})

test_that("comparator execution errors are not converted into fallback gates", {
  env <- .load_acs_method_env()
  env$.simCompareFbetaThreshold <- function(...) {
    stop("reticulate worker failure")
  }

  expect_error(
    env$.acsCytofThresholdOne(
      method = "fbeta",
      xUns = c(0, 1),
      xStim = c(0, 2),
      settings = env$.acsCytofComparatorSettings("fbeta"),
      fbetaEnv = new.env(parent = emptyenv())
    ),
    "reticulate worker failure"
  )
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

test_that("requested comparator runs remove both previous result files", {
  env <- .load_acs_method_env()
  path_dir <- tempfile("acs-comparator-results-")
  dir.create(path_dir, recursive = TRUE)
  withr::defer(unlink(path_dir, recursive = TRUE))

  paths <- list(
    fbeta = file.path(path_dir, "fbeta", "result.rds"),
    tailgate = file.path(path_dir, "tailgate", "result.rds"),
    gs = file.path(path_dir, "gs")
  )
  dir.create(dirname(paths$fbeta), recursive = TRUE)
  dir.create(dirname(paths$tailgate), recursive = TRUE)
  saveRDS("old fbeta", paths$fbeta)
  saveRDS("old tailgate", paths$tailgate)
  writeLines("keep", paths$gs)

  removed <- env$.acsCytofRemoveComparatorResults(paths)

  expect_setequal(removed, c(paths$fbeta, paths$tailgate))
  expect_false(file.exists(paths$fbeta))
  expect_false(file.exists(paths$tailgate))
  expect_true(file.exists(paths$gs))
})

test_that("ACS methods stay sequential and always replace previous results", {
  env <- .load_acs_method_env()
  runner_body <- paste(
    deparse(body(env$.acsCytofRunComparisonMethods)),
    collapse = "\n"
  )

  expect_true(grepl("lapply(methodVec", runner_body, fixed = TRUE))
  expect_false(grepl("future_map", runner_body, fixed = TRUE))
  expect_true(grepl(".acsCytofWriteComparatorCache", runner_body, fixed = TRUE))
  expect_true(grepl(
    ".acsCytofRemoveComparatorResults",
    runner_body,
    fixed = TRUE
  ))
  expect_false(grepl("Using cached", runner_body, fixed = TRUE))
  expect_false("overwrite" %in% names(formals(
    env$.acsCytofRunComparisonMethods
  )))
})

test_that("ACS comparator populations run in parallel", {
  qmd_lines <- readLines(qmd_path, warn = FALSE)
  expect_false(any(grepl(
    "overwrite_comparator_cache",
    qmd_lines,
    fixed = TRUE
  )))
  chunk_start <- grep(
    "#| label: run-comparison-methods-in-parallel",
    qmd_lines,
    fixed = TRUE
  )
  expect_length(chunk_start, 1L)

  chunk_end_offset <- which(
    qmd_lines[(chunk_start + 1L):length(qmd_lines)] == "```"
  )[1L]
  expect_false(is.na(chunk_end_offset))

  chunk_lines <- qmd_lines[
    chunk_start:(chunk_start + chunk_end_offset)
  ]
  chunk_body <- paste(chunk_lines, collapse = "\n")

  expect_true(grepl("furrr::future_map", chunk_body, fixed = TRUE))
  expect_true(grepl(
    ".acsCytofRunComparisonMethodsSafe",
    chunk_body,
    fixed = TRUE
  ))
  expect_true(grepl("future::multisession", chunk_body, fixed = TRUE))
  expect_true(grepl("pathFbeta = path_fbeta", chunk_body, fixed = TRUE))
  expect_true(grepl(
    "finally = future::plan(old_comparator_plan)",
    chunk_body,
    fixed = TRUE
  ))
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
