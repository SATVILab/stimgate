root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
script_runtime <- file.path(root_dir, "scripts", "r", "analysis-runtime.R")
script_helper <- file.path(root_dir, "scripts", "r", "acs_cytof-helper.R")
script_gate <- file.path(root_dir, "scripts", "r", "acs_cytof-gate.R")
script_methods <- file.path(root_dir, "scripts", "r", "acs_cytof-methods.R")
script_manual <- file.path(root_dir, "scripts", "r", "acs_cytof-manual.R")
script_compare <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")
script_fbeta <- file.path(root_dir, "scripts", "python", "fbeta.py")
qmd_path <- file.path(root_dir, "analysis", "9-real-compare-acs-cytof.qmd")

.load_acs_method_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_runtime, local = env)
  source(script_helper, local = env)
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

test_that("Tailgate execution errors are not converted into fallback gates", {
  env <- .load_acs_method_env()
  env$.simCompareTailgateThreshold <- function(...) {
    stop("tailgate worker failure")
  }

  expect_error(
    env$.acsCytofThresholdOne(
      method = "tailgate",
      xUns = c(0, 1),
      xStim = c(0, 2),
      settings = env$.acsCytofComparatorSettings("tailgate")
    ),
    "tailgate worker failure"
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

test_that("a failed comparator rerun keeps the previous results", {
  env <- .load_acs_method_env()
  path_base <- tempfile("acs-comparator-results-")
  withr::defer(unlink(path_base, recursive = TRUE))
  paths <- env$.acsCytofPopulationPaths(
    pop = "cd4",
    pathFcsBase = file.path(path_base, "fcs"),
    pathGsBase = file.path(path_base, "gs"),
    pathScratchBase = file.path(path_base, "scratch")
  )
  dir.create(paths$gs, recursive = TRUE)
  dir.create(dirname(paths$fbeta), recursive = TRUE)
  dir.create(dirname(paths$tailgate), recursive = TRUE)
  saveRDS("old fbeta", paths$fbeta)
  saveRDS("old tailgate", paths$tailgate)

  testthat::local_mocked_bindings(
    load_gs = function(...) as.list(seq_len(10L)),
    .package = "flowWorkspace"
  )
  run_with <- function(fail_method) {
    env$.acsCytofRunComparator <- function(gs, pop, method, ...) {
      if (identical(method, fail_method)) stop("boom")
      list(method = method, new = TRUE)
    }
    env$.acsCytofRunComparisonMethods(
      pop = "cd4",
      pathFcsBase = file.path(path_base, "fcs"),
      pathGsBase = file.path(path_base, "gs"),
      pathScratchBase = file.path(path_base, "scratch"),
      runMethods = TRUE
    )
  }

  expect_error(run_with("tailgate"), "boom")
  expect_identical(readRDS(paths$fbeta), "old fbeta")
  expect_identical(readRDS(paths$tailgate), "old tailgate")

  expect_no_error(run_with("none"))
  expect_true(readRDS(paths$fbeta)$new)
  expect_true(readRDS(paths$tailgate)$new)
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
  expect_false(grepl("Remove", runner_body, fixed = TRUE))
  expect_false(grepl("unlink", runner_body, fixed = TRUE))
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

  expect_true(grepl(".acsCytofMapPopulations(", chunk_body, fixed = TRUE))
  expect_true(grepl(
    ".acsCytofRunComparisonMethodsSafe",
    chunk_body,
    fixed = TRUE
  ))
  expect_true(grepl("seed = analysis_seed", chunk_body, fixed = TRUE))
  expect_true(grepl("pathFbeta = path_fbeta", chunk_body, fixed = TRUE))

  helper_body <- paste(
    deparse(body(.load_acs_method_env()$.acsCytofMapPopulations)),
    collapse = "\n"
  )
  expect_true(grepl("furrr::future_map", helper_body, fixed = TRUE))
  expect_true(grepl("future::multisession", helper_body, fixed = TRUE))
  expect_true(grepl("future::plan(oldPlan)", helper_body, fixed = TRUE))
})

test_that("FCS files missing from the manual sample map are dropped with a warning", {
  skip_if_not_installed("DataTidyACSCyTOFFAUST")
  env <- .load_acs_method_env()
  path_fcs <- tempfile("acs-fcs-")
  withr::defer(unlink(path_fcs, recursive = TRUE))
  dir.create(file.path(path_fcs, "cd4"), recursive = TRUE)
  fcs_names <- sprintf("file%02d.fcs", 1:10)
  file.create(file.path(path_fcs, "cd4", fcs_names))
  clean <- DataTidyACSCyTOFFAUST::clean_fcs_for_matching(fcs_names)
  lookup <- tibble::tibble(
    MatchFCSName = clean[-3L],
    SampleID = paste0("s", 1:9),
    Stim = "mtb"
  )

  expect_warning(
    out <- env$.acsCytofManualSampleMapFromFcs(path_fcs, "cd4", lookup),
    "file03.fcs"
  )
  expect_equal(nrow(out), 9L)
  expect_false("3" %in% out$ind)
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


test_that("manual comparison save preserves the last good RDS on pre-save failure", {
  env <- .load_acs_method_env()
  comparison_tbl <- tibble::tibble(
    method = c("stimgate", "stimgate"),
    pop = c("CD4 T cells", "CD4 T cells"),
    cyt = c("IFNg", "IFNg"),
    stim = c("mtb", "ebv"),
    freq_bs_auto = c(1.0, 2.0),
    freq_bs_man = c(1.1, 1.8)
  ) |>
    dplyr::mutate(
      diff = .data$freq_bs_auto - .data$freq_bs_man,
      abs_diff = abs(.data$diff),
      rel_error = .data$diff / .data$freq_bs_man,
      abs_rel_error = abs(.data$rel_error)
    )

  path_dir <- tempfile("acs-manual-output-")
  withr::defer(unlink(path_dir, recursive = TRUE))

  expect_no_error(env$.acsCytofManualSave(
    comparisonTbl = comparison_tbl,
    pathDirSave = path_dir,
    savePlots = FALSE
  ))
  path_rds <- file.path(path_dir, "manual-comparison.rds")
  expect_equal(readRDS(path_rds), comparison_tbl)

  expect_error(env$.acsCytofManualSave(
    comparisonTbl = tibble::tibble(),
    pathDirSave = path_dir,
    savePlots = FALSE
  ))
  expect_equal(readRDS(path_rds), comparison_tbl)
})

test_that("analysis 9 builds before saving to the canonical manual output", {
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("path_dir_save = NULL", content, fixed = TRUE))
  expect_true(grepl(".acsCytofManualSave(", content, fixed = TRUE))
  expect_true(grepl(
    'path_manual_output <- projr::projr_path_get_dir(',
    content,
    fixed = TRUE
  ))
  expect_false(grepl(
    'path_manual_output <- file.path(\n  path_scratch_base',
    content,
    fixed = TRUE
  ))
  expect_false(grepl(
    "unlink(path_manual_output, recursive = TRUE)",
    content,
    fixed = TRUE
  ))

  save_body <- paste(
    deparse(body(.load_acs_method_env()$.acsCytofManualSave)),
    collapse = "\n"
  )
  expect_true(grepl(".write_rds_atomic(", save_body, fixed = TRUE))
})
