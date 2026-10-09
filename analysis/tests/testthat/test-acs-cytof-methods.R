root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
script_runtime <- file.path(root_dir, "scripts", "r", "analysis-runtime.R")
script_helper <- file.path(root_dir, "scripts", "r", "acs_cytof-helper.R")
script_gate <- file.path(root_dir, "scripts", "r", "acs_cytof-gate.R")
script_methods <- file.path(root_dir, "scripts", "r", "acs_cytof-methods.R")
script_manual <- file.path(root_dir, "scripts", "r", "acs_cytof-manual.R")
script_style <- file.path(root_dir, "scripts", "r", "analysis-plot-style.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")
script_plot_cyt <- file.path(root_dir, "scripts", "r", "acs_cytof-plot_cyt.R")
script_compare <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")
script_fbeta <- file.path(root_dir, "scripts", "python", "fbeta.py")
qmd_path <- file.path(root_dir, "analysis", "9-real-compare-acs-cytof.qmd")

.load_acs_method_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_runtime, local = env)
  source(script_style, local = env)
  source(script_bw_plot, local = env)
  source(script_helper, local = env)
  source(script_gate, local = env)
  source(script_methods, local = env)
  source(script_manual, local = env)
  source(script_plot_cyt, local = env)
  env
}

test_that("ACS comparator settings pin the tuned and default comparators", {
  env <- .load_acs_method_env()
  expect_setequal(
    env$.acsCytofComparatorMethods(),
    c("fbeta", "tailgate", "fbeta_default", "tailgate_default")
  )

  for (method in c("fbeta", "fbeta_default")) {
    fbeta <- env$.acsCytofComparatorSettings(method)
    expect_identical(fbeta$family, "fbeta")
    expect_equal(fbeta$params$beta, 0.8)
    expect_equal(fbeta$params$theta, 2)
    expect_equal(fbeta$params$width, 10L)
    expect_null(fbeta$params$numBins)
  }
  expect_true(env$.acsCytofComparatorSettings("fbeta")$params$removeZero)
  expect_false(env$.acsCytofComparatorSettings("fbeta_default")$params$removeZero)
  expect_equal(env$.acsCytofComparatorSettings("fbeta")$cacheVersion, 3L)
  expect_equal(env$.acsCytofComparatorSettings("fbeta_default")$cacheVersion, 2L)

  for (method in c("tailgate", "tailgate_default")) {
    tailgate <- env$.acsCytofComparatorSettings(method)
    expect_identical(tailgate$family, "tailgate")
    expect_identical(tailgate$params$tailgateX, "stim")
    expect_null(tailgate$params$bandwidth)
    expect_identical(tailgate$params$derivativeMethod, "firstDeriv")
  }
  tuned <- env$.acsCytofComparatorSettings("tailgate")$params
  expect_false(tuned$autoTol)
  expect_equal(tuned$tol, 2e-5)
  expect_equal(tuned$bias, 0.2)
  expect_true(tuned$removeZero)
  expect_equal(env$.acsCytofComparatorSettings("tailgate")$cacheVersion, 2L)

  default <- env$.acsCytofComparatorSettings("tailgate_default")$params
  expect_true(default$autoTol)
  expect_equal(default$bias, 0)
  expect_false(default$removeZero)
  expect_equal(env$.acsCytofComparatorSettings("tailgate_default")$cacheVersion, 1L)
})

test_that("Tailgate moves its cutpoint up by the bias and can drop exact zeros", {
  skip_if_not_installed("cytoUtils")
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_compare, local = env)
  set.seed(1)
  x <- c(rep(0, 3000), stats::rnorm(2000, 1), stats::rnorm(40, 5))

  base <- env$.simCompareTailgateThreshold(x, tol = 1e-3)
  shifted <- env$.simCompareTailgateThreshold(x, tol = 1e-3, bias = 0.2)
  expect_equal(shifted$threshold, base$threshold + 0.2)

  no_zero <- env$.simCompareTailgateThreshold(x, tol = 1e-3, removeZero = TRUE)
  expect_equal(
    no_zero$threshold,
    env$.simCompareTailgateThreshold(x[x != 0], tol = 1e-3)$threshold
  )
  expect_error(env$.simCompareTailgateThreshold(x, bias = NA), "bias")
})

test_that("F-beta without zeros scales each pdf to its retained fraction", {
  skip_if_not_installed("reticulate")
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_compare, local = env)
  fbeta_env <- env$.simCompareFbetaEnvironment(pathFbeta = script_fbeta)
  set.seed(2)
  x_uns <- c(rep(0, 600), stats::rnorm(400, 1))
  x_stim <- c(rep(0, 300), stats::rnorm(650, 1), stats::rnorm(50, 4))

  out <- env$.simCompareFbetaThreshold(
    xUns = x_uns, xStim = x_stim, fbetaEnv = fbeta_env, removeZero = TRUE
  )
  direct <- fbeta_env$get_positivity_threshold(
    neg = matrix(x_uns[x_uns != 0], ncol = 1L),
    pos = matrix(x_stim[x_stim != 0], ncol = 1L),
    channelIndex = 0L, beta = 0.8, theta = 2, width = 10L, numBins = NULL,
    negScale = 0.4, posScale = 0.7
  )
  expect_equal(out$threshold, as.numeric(direct$threshold))
  # The retained pdfs integrate to the retained fractions, not to one.
  unscaled <- fbeta_env$get_positivity_threshold(
    neg = matrix(x_uns[x_uns != 0], ncol = 1L),
    pos = matrix(x_stim[x_stim != 0], ncol = 1L),
    channelIndex = 0L, beta = 0.8, theta = 2, width = 10L, numBins = NULL
  )
  expect_equal(as.numeric(out$fbeta$pdfneg), 0.4 * as.numeric(unscaled$pdfneg))
  expect_equal(as.numeric(out$fbeta$pdfpos), 0.7 * as.numeric(unscaled$pdfpos))

  # Without zeros, dropping them changes nothing.
  x_uns_pos <- x_uns[x_uns != 0]
  x_stim_pos <- x_stim[x_stim != 0]
  expect_equal(
    env$.simCompareFbetaThreshold(x_uns_pos, x_stim_pos, fbetaEnv = fbeta_env,
      removeZero = TRUE)$threshold,
    env$.simCompareFbetaThreshold(x_uns_pos, x_stim_pos, fbetaEnv = fbeta_env)$threshold
  )
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

test_that("comparator execution errors stay explicit missing outcomes, not fallback gates", {
  env <- .load_acs_method_env()
  env$.simCompareFbetaThreshold <- function(...) stop("reticulate worker failure")
  env$.simCompareTailgateThreshold <- function(...) stop("tailgate worker failure")
  # Tailgate checks its dependency before estimating.
  methods <- c("fbeta", if (requireNamespace("cytoUtils", quietly = TRUE)) "tailgate")
  for (method in methods) {
    expect_message(
      out <- env$.acsCytofThresholdOne(
        method = method,
        xUns = c(0, 1),
        xStim = c(0, 2),
        settings = env$.acsCytofComparatorSettings(method),
        fbetaEnv = new.env(parent = emptyenv())
      ),
      "worker failure"
    )
    expect_identical(out$thresholdOrigin, "runtime_error")
    expect_true(is.na(out$threshold))
    expect_false(out$thresholdFallbackUsed)
  }
})

test_that("a missing comparator implementation still stops the analysis", {
  env <- .load_acs_method_env()
  expect_error(
    env$.acsCytofThresholdOne(
      method = "tailgate",
      xUns = c(0, 1),
      xStim = c(0, 2),
      settings = env$.acsCytofComparatorSettings("tailgate")
    ),
    "Source scripts/r/sim-compare-freq_bs.R"
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
  env$.acsCytofReadPreprocessing <- function(...) list(sampleMap = data.frame(
    SampleID = rep(c("a", "b"), each = 5),
    stim = rep(c("uns", "p1", "mtb", "ebv", "p4"), 2),
    ind = 1:10
  ))
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

  expect_false(file.exists(paths$fbeta_default))
  expect_false(file.exists(paths$tailgate_default))

  expect_no_error(run_with("none"))
  for (method in env$.acsCytofComparatorMethods()) {
    expect_true(readRDS(paths[[method]])$new)
  }
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

test_that("FCS files missing from the manual sample map fail explicitly", {
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

  expect_error(
    env$.acsCytofManualSampleMapFromFcs(path_fcs, "cd4", lookup),
    "file03.fcs"
  )
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
    SampleID = c("a", "b"),
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
    'path_manual_output <- projr::projr_path_get(',
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
    deparse(body(.load_acs_method_env()$.acsCytofManualWrite)),
    collapse = "\n"
  )
  expect_true(grepl(".write_rds_atomic(", save_body, fixed = TRUE))
})

test_that("analysis 9 reads the canonical comparison without rebuilding raw inputs", {
  lines <- readLines(qmd_path, warn = FALSE)
  start <- which(lines == "#| label: format-manual-comparison")
  end <- which(lines == "```" & seq_along(lines) > start)[1L]
  code <- parse(text = lines[seq.int(start + 1L, end - 1L)])
  # Supply the configured cache path without requiring projr or external data.
  expect_identical(as.character(code[[1]][[2]]), "path_manual_output")
  env <- .load_acs_method_env()
  source(file.path(root_dir, "scripts", "r", "analysis-runtime.R"), local = env)
  env$analysis_qmd <- "analysis/9-real-compare-acs-cytof.qmd"
  env$path_manual_output <- tempfile("acs-manual-cache-")
  dir.create(env$path_manual_output)
  withr::defer(unlink(env$path_manual_output, recursive = TRUE))
  env$run_preprocessing <- FALSE
  env$run_stimgate <- FALSE
  env$run_comparators <- FALSE
  env$stimgate_loc_threshold_method <- "region"
  env$comp_against_manual_cyt <- function(...) stop("Raw-data rebuild was called")
  env$.acsCytofManualSave <- function(...) stop("Cache write was called")
  env$.acsCytofManualSummaryTable <- function(x) x
  cached <- tibble::tibble(
    method = "stimgate", freq_bs_auto = 0.1, thresholdFailed = FALSE,
    locThresholdMethod = "region"
  )
  attr(cached, "manifest") <- list(
    methods = list(cd4 = list(stimgate = list(
      context = list(gitSha = "abc"),
      settings = list(locThresholdMethod = "region"),
      channelSettings = list(IFNg = list(locThresholdMethod = "region"))
    ))),
    comparisonSettings = list(methods = "stimgate", locThresholdMethod = "region"),
    manualInputHash = "abc"
  )
  path <- file.path(env$path_manual_output, "manual-comparison.rds")
  saveRDS(cached, path)
  path_run_manifest <- file.path(env$path_manual_output, "manifest.rds")
  saveRDS(list(analysis_semantics_version = "acs-cytof-v8",
    stimgate_loc_threshold_method = "region"), path_run_manifest)
  withr::local_envvar(ANALYSIS_EXPECTED_RUN_ID = NA_character_)
  for (expr in as.list(code)[-1L]) eval(expr, env)
  expect_identical(env$manual_comparison_tbl, cached)
  expect_identical(env$manual_summary_tbl, cached)

  # Caches from before the threshold method was recorded are rejected.
  for (old in list(
    list(analysis_semantics_version = "acs-cytof-v3",
      stimgate_loc_threshold_method = "region"),
    list(analysis_semantics_version = "acs-cytof-v8")
  )) {
    saveRDS(old, path_run_manifest)
    expect_error(
      for (expr in as.list(code)[-1L]) eval(expr, env),
      "RUN_SIMULATIONS=true"
    )
  }
  saveRDS(list(analysis_semantics_version = "acs-cytof-v8",
    stimgate_loc_threshold_method = "region"), path_run_manifest)

  unlink(path)
  expect_error(
    for (expr in as.list(code)[-1L]) eval(expr, env),
    "RUN_SIMULATIONS=true RUN_PLOTS=false SIM_SIZE=final quarto render analysis/9"
  )
})

test_that("ACS manual-comparison plots use method colours and no titles", {
  env <- .load_acs_method_env()
  tbl <- tidyr::expand_grid(
    method = c("stimgate", "tailgate", "fbeta"),
    pop = c("CD4 T cells", "B cells"),
    cyt = c("IFNg", "IL2"),
    stim = c("mtb", "p1"),
    sample = 1:6
  ) |>
    dplyr::mutate(
      freq_bs_man = 0.1 * sample,
      freq_bs_auto = freq_bs_man * ifelse(sample %% 2 == 0, 1.5, 0.5),
      rel_error = (freq_bs_auto - freq_bs_man) / freq_bs_man,
      abs_rel_error = abs(rel_error)
    )
  plots <- list(
    scatter = env$.acsCytofManualPlotScatter(dplyr::filter(tbl, pop == "CD4 T cells")),
    relative = env$.acsCytofManualPlotRelativeError(tbl),
    signed = env$.acsCytofManualPlotSignedError(tbl)
  )
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_null(p$labels$title)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
  # Signed error summarises both directions for each method.
  expect_setequal(as.character(plots$signed$data$direction), c("over", "under"))
  expect_setequal(
    as.character(plots$signed$data$method),
    c("stimgate", "tailgate", "fbeta")
  )
  without_tg <- env$.acsCytofManualPlotSignedError(
    dplyr::filter(tbl, method != "tailgate")
  )
  expect_setequal(as.character(without_tg$data$method), c("stimgate", "fbeta"))
})

test_that("ACS method sets cover all methods and without Tailgate", {
  env <- .load_acs_method_env()
  sets <- env$.acsCytofMethodSets()
  expect_named(sets, c("all_methods", "no_tailgate"))
  expect_identical(sets$all_methods$label, "All methods")
  expect_identical(sets$no_tailgate$label, "Without Tailgate")
  expect_false("tailgate" %in% sets$no_tailgate$methods)
})

test_that("ACS finite fallback gates remain failures in comparison and coverage", {
  env <- .load_acs_method_env()
  single <- tibble::tibble(
    method = "stimgate", ind = c("2", "3", "4"), cyt = "IFNg",
    freq_stim_auto = c(1, 2, 0), freq_uns_auto = 0,
    freq_bs_auto = c(1, 2, 0), freq_bs_man = c(1, 2, 10),
    pop = "CD4 T cells", stim = "mtb"
  )
  thresholds <- tibble::tibble(
    ind = c("2", "3", "4"), cyt = "IFNg", threshold = c(1, 2, 100),
    thresholdOrigin = c("direct", "cluster", "high_value"),
    thresholdFallbackUsed = c(FALSE, FALSE, TRUE),
    locGenerated = c(TRUE, TRUE, FALSE),
    locGeneratedDirect = c(TRUE, FALSE, FALSE),
    locSource = thresholdOrigin, locReason = c(NA, NA, "no_threshold")
  )
  out <- env$.acsCytofJoinProvenance(single, thresholds) |>
    dplyr::mutate(diff = freq_bs_auto - freq_bs_man, abs_diff = abs(diff),
                  abs_rel_error = abs(diff / freq_bs_man))
  expect_true(is.na(out$freq_bs_auto[3]))
  expect_identical(out$locReason[3], "no_threshold")
  coverage <- env$.acsCytofManualSummaryTable(out)
  expect_equal(coverage$n_total, 3L)
  expect_equal(coverage$n_failed, 1L)
  expect_equal(coverage$n, 2L)
  expect_equal(coverage$pcc, 1)
  expect_error(env$.acsCytofJoinProvenance(single, thresholds[-1, ]), "Missing or duplicate")
})

test_that("ACS refuses mismatched data and preprocessing manifests", {
  env <- .load_acs_method_env()
  context <- list(gitSha = "abc", preprocessing = list(
    settings = list(transform = "asinh(x / 5)"), inputFileListHash = "files"
  ))
  manifest <- list(context = context, settings = list(clusterGates = TRUE))
  expect_no_error(env$.acsCytofValidateManifests(list(manifest, manifest)))
  expect_no_error(env$.acsCytofValidateManifests(list(manifest,
    list(context = list(gitSha = "other", preprocessing = context$preprocessing)))))
  for (changed in list(
    list(gitSha = "abc", preprocessing = list(inputFileListHash = "other")),
    list(gitSha = "abc", preprocessing = list(settings = list(transform = "none")))
  )) {
    expect_error(env$.acsCytofValidateManifests(list(manifest, list(context = changed))), "Mismatched ACS")
  }
  expect_error(env$.acsCytofValidateManifests(list(manifest, NULL)), "Mismatched ACS")
})

test_that("ACS compares identical cohorts and rejects missing or duplicate strata", {
  env <- .load_acs_method_env()
  rows <- tibble::tibble(method = rep(c("stimgate", "fbeta"), each = 2),
                         SampleID = rep(c("a", "b"), 2),
                         stim = "mtb", pop = "CD4 T cells", cyt = "IFNg")
  expect_no_error(env$.acsCytofValidateCohorts(rows, c("stimgate", "fbeta")))
  expect_error(env$.acsCytofValidateCohorts(rows[-4, ], c("stimgate", "fbeta")), "different or duplicate")
  expect_error(env$.acsCytofValidateCohorts(dplyr::bind_rows(rows, rows[4, ]), c("stimgate", "fbeta")), "different or duplicate")
})

test_that("ACS detects reordered saved GatingSet files", {
  env <- .load_acs_method_env()
  path <- tempfile("acs-manifest-")
  withr::defer(unlink(env$.acsCytofPreprocessingFile(path)))
  map <- data.frame(SampleID = "a", stim = c("uns", "p1", "mtb", "ebv", "p4"),
                    ind = 1:5, file = paste0(1:5, ".fcs"))
  saveRDS(list(sampleMap = map), env$.acsCytofPreprocessingFile(path))
  testthat::local_mocked_bindings(sampleNames = function(...) rev(map$file),
                                  .package = "flowWorkspace")
  expect_error(env$.acsCytofReadPreprocessing(path, list()), "reordered")
})

test_that("ACS preprocessing manifest does not stop load_gs() reading its GatingSet", {
  env <- .load_acs_method_env()
  path <- tempfile("acs-gs-")
  withr::defer(unlink(c(path, env$.acsCytofPreprocessingFile(path)), recursive = TRUE))
  fr <- flowCore::flowFrame(matrix(1:4, ncol = 2, dimnames = list(NULL, c("A", "B"))))
  fs <- flowCore::flowSet(list(s1 = fr))
  flowWorkspace::save_gs(flowWorkspace::GatingSet(flowWorkspace::flowSet_to_cytoset(fs)), path)
  saveRDS(list(sampleMap = NULL), env$.acsCytofPreprocessingFile(path))
  gs <- flowWorkspace::load_gs(path)
  expect_identical(flowWorkspace::sampleNames(gs), "s1")
  expect_true(file.exists(env$.acsCytofPreprocessingFile(path)))
})

test_that("ACS coverage excludes every failed estimate even when a fallback is finite", {
  env <- .load_acs_method_env()
  rows <- tibble::tibble(method = "fbeta", pop = "CD4 T cells", cyt = "IFNg", stim = "mtb",
                         freq_bs_auto = c(0, 0), freq_bs_man = c(1, 2),
                         thresholdFailed = TRUE, abs_diff = c(1, 2), abs_rel_error = 1)
  summary <- env$.acsCytofManualSummaryTable(rows)
  expect_equal(summary$n, 0L)
  expect_equal(summary$n_failed, 2L)
  expect_true(is.na(summary$pcc))
  expect_true(is.na(summary$mae))
})

test_that("ACS saves exclusions and manifests with the comparison transaction", {
  env <- .load_acs_method_env()
  rows <- tibble::tibble(SampleID = "a", method = "fbeta", pop = "CD4 T cells", cyt = "IFNg", stim = "mtb",
                         freq_bs_auto = 1, freq_bs_man = 1,
                         thresholdFailed = FALSE, abs_diff = 0, abs_rel_error = 0)
  attr(rows, "manifest") <- list(manualInputHash = "abc", methods = list())
  attr(rows, "exclusions") <- tibble::tibble(SampleID = "excluded", exclusionReason = "no_manual_key")
  path <- tempfile("acs-save-")
  withr::defer(unlink(path, recursive = TRUE))
  env$.acsCytofManualSave(rows, path, FALSE)
  expect_identical(readRDS(file.path(path, "acs-manifest.rds")), attr(rows, "manifest"))
  expect_identical(readRDS(file.path(path, "manual-comparison.rds")), rows)
  exclusions <- utils::read.csv(file.path(path, "manual-comparison-exclusions.csv"))
  expect_equal(exclusions$SampleID, "excluded")
  env$.acsCytofManualWrite <- function(...) stop("failed replacement")
  expect_error(env$.acsCytofManualSave(rows, path, FALSE), "failed replacement")
  expect_identical(readRDS(file.path(path, "manual-comparison.rds")), rows)
})

test_that("ACS cached comparisons reject legacy and mixed method manifests", {
  env <- .load_acs_method_env()
  table <- tibble::tibble(thresholdFailed = FALSE)
  expect_error(env$.acsCytofValidateComparisonManifest(table), "Legacy or incomplete")
  attr(table, "manifest") <- list(
    methods = list(cd4 = list(stimgate = list(context = list(gitSha = "a")),
                             fbeta = list(context = list(gitSha = "b", preprocessing = list(inputFileListHash = "other"))))),
    comparisonSettings = list(methods = c("stimgate", "fbeta"), locThresholdMethod = "region"),
    manualInputHash = "abc"
  )
  expect_error(env$.acsCytofValidateComparisonManifest(table), "Mismatched ACS")
})


test_that("ACS error denominators retain zero manual values only for absolute errors", {
  env <- .load_acs_method_env()
  rows <- tibble::tibble(SampleID = letters[1:4], method = "stimgate", pop = "CD4", cyt = "IFNg", stim = "mtb",
    freq_bs_auto = c(1, 2, 3, NA_real_), freq_bs_man = c(0, -1, 1, 1),
    abs_diff = c(1, 3, 2, NA_real_), abs_rel_error = c(NA_real_, NA_real_, 2, NA_real_), thresholdFailed = c(FALSE, FALSE, FALSE, TRUE))
  out <- env$.acsCytofManualSummaryTable(rows)
  expect_equal(out$n, 3)
  expect_equal(out$n_relative, 1)
  expect_equal(out$n_relative_excluded, 3)
  expect_equal(out$n_manual_nonpositive, 2)
  expect_equal(out$mae, 2)
  uncertainty <- env$.acsCytofManualUncertainty(rows, reps = 99)
  expect_equal(uncertainty$n_donors, 3)
  expect_equal(uncertainty$n_relative_donors, 1)
  expect_true(is.na(uncertainty$mean_abs_rel_error_lower))
})

test_that("ACS scientific manifests work when a diagnostic git revision is unavailable", {
  env <- .load_acs_method_env()
  env$.git_sha <- function(...) NA_character_
  preprocessing <- list(settings = list(transform = "asinh(x / 5)"))
  manifest <- env$.acsCytofManifest(preprocessing)
  expect_true(is.na(manifest$gitSha))
  expect_identical(manifest$preprocessing, preprocessing)
  expect_no_error(env$.acsCytofValidateManifests(list(
    list(context = manifest), list(context = manifest)
  )))
})

test_that("ACS bootstrap keeps common donor draws across tubes, methods and stimuli", {
  env <- .load_acs_method_env()
  rows <- tidyr::expand_grid(SampleID = letters[1:4], method = c("stimgate", "fbeta"), stim = c("mtb", "p4")) |>
    dplyr::mutate(pop = "CD4", cyt = "IFNg", freq_bs_man = 1,
      freq_bs_auto = match(SampleID, letters), thresholdFailed = FALSE)
  set.seed(43)
  seed <- .Random.seed
  out <- env$.acsCytofManualUncertainty(rows, reps = 99)
  expect_identical(.Random.seed, seed)
  expect_equal(dplyr::n_distinct(out$mae_lower), 1L)
  expect_equal(dplyr::n_distinct(out$mae_upper), 1L)
  duplicated <- env$.acsCytofManualUncertainty(dplyr::bind_rows(rows, rows), reps = 99)
  expect_equal(out, duplicated)
  expect_error(env$.acsCytofManualUncertainty(rows, reps = 0), "at least two")
})

test_that("ACS gate diagnostics remove panel titles before arranging plots", {
  env <- .load_acs_method_env()
  panel <- ggplot2::ggplot(data.frame(x = 1), ggplot2::aes(x, x)) +
    ggplot2::labs(title = "Sample 2", subtitle = "Diagnostic")
  testthat::local_mocked_bindings(
    plotStim = function(..., grid) {
      expect_false(grid)
      list(panel, panel)
    }, .package = "stimgate"
  )
  captured <- NULL
  testthat::local_mocked_bindings(
    plot_grid = function(plotlist, ...) {
      captured <<- plotlist
      ggplot2::ggplot()
    }, .package = "cowplot"
  )
  expect_no_error(env$.acsCytofPlotGateCheck(NULL, "unused"))
  expect_length(captured, 2L)
  for (p in captured) {
    expect_null(p$labels$title)
    expect_null(p$labels$subtitle)
    expect_identical(p$mapping, panel$mapping)
  }
  # The original package plot remains available with its labels.
  expect_identical(panel$labels$title, "Sample 2")
  qmd <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
  expect_match(qmd, "stimgate_check_sample_2.pdf", fixed = TRUE)
})

test_that("ACS log1p boxplot display preserves raw error statistics and zeros", {
  env <- .load_acs_method_env()
  tbl <- tidyr::expand_grid(
    method = c("stimgate", "fbeta"),
    pop = c("CD4 T cells", "CD8 T cells"),
    cyt = c("IFNg", "IL2"),
    sample = 1:8
  ) |>
    dplyr::mutate(
      abs_rel_error = c(0, 0.01, 0.02, 0.1, 0.2, 0.3, 1, 1000)[sample] *
        ifelse(cyt == "IL2", 0.01, 1)
    )
  p <- env$.acsCytofManualPlotRelativeError(tbl)
  built <- ggplot2::ggplot_build(p)
  # Replacing the display coordinates leaves every box statistic identical.
  linear <- suppressMessages(p + ggplot2::coord_cartesian())
  expect_equal(built$data[[1]], ggplot2::ggplot_build(linear)$data[[1]])
  expect_equal(p$coordinates$trans$y$transform(c(0, 1, 1000)), log1p(c(0, 1, 1000)))
  expect_equal(length(unique(built$layout$layout$SCALE_Y)), 4L)
  expect_equal(p$data$abs_rel_error, tbl$abs_rel_error)
  expect_match(p$labels$y, "log1p scale; ticks in original units", fixed = TRUE)
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("ACS manual scatter has smaller stimulus points and independent cytokine axes", {
  env <- .load_acs_method_env()
  tbl <- tidyr::expand_grid(
    method = c("stimgate", "fbeta"), pop = "CD4 T cells",
    cyt = c("IFNg", "IL2"), stim = c("mtb", "p1"), sample = 1:4
  ) |>
    dplyr::mutate(freq_bs_man = sample / 10, freq_bs_auto = sample / 8)
  p <- env$.acsCytofManualPlotScatter(tbl)
  built <- ggplot2::ggplot_build(p)
  expect_setequal(built$data[[2]]$shape, c(16, 17))
  expect_true(all(built$data[[2]]$alpha == 0.5))
  expect_true(all(built$data[[2]]$size == 0.9))
  expect_equal(length(unique(built$layout$layout$SCALE_X)), 2L)
  expect_equal(length(unique(built$layout$layout$SCALE_Y)), 2L)
  expect_equal(p$facet$params$ncol, 3)
  expect_equal(p$data$freq_bs_auto, tbl$freq_bs_auto)
})

test_that("analysis 9 scatter chunk prints and saves each method set and population", {
  env <- .load_acs_method_env()
  lines <- readLines(qmd_path, warn = FALSE)
  start <- which(lines == "#| label: manual-comparison")
  end <- which(lines == "```" & seq_along(lines) > start)[1L]
  code <- parse(text = lines[seq.int(start + 1L, end - 1L)])
  env$manual_comparison_tbl <- tidyr::expand_grid(
    method = c("stimgate", "tailgate", "fbeta"),
    pop = c("CD4 T cells", "CD8 T cells"), cyt = "IFNg", stim = "mtb"
  ) |>
    dplyr::mutate(freq_bs_man = 0.1, freq_bs_auto = 0.12)
  plots <- paths <- list()
  env$.analysis_print_save_fig <- function(p, path, ...) {
    plots[[length(plots) + 1L]] <<- p
    paths[[length(paths) + 1L]] <<- path
  }
  env$.analysis_fig_dir <- function(parts, ...) do.call(file.path, as.list(parts))
  env$fig_key <- "analysis9"
  env$root_dir <- root_dir
  env$run_plots <- TRUE
  out <- capture.output(eval(code, env))
  expect_length(plots, 4L)
  expect_equal(length(unique(unlist(paths))), 4L)
  expect_true(any(grepl("##### Population: CD4 T cells", out, fixed = TRUE)))
  expect_true(all(vapply(plots, function(p) dplyr::n_distinct(p$data$pop) == 1L, logical(1))))
  expect_setequal(as.character(plots[[1]]$data$method), c("stimgate", "tailgate", "fbeta"))
  expect_setequal(as.character(plots[[3]]$data$method), c("stimgate", "fbeta"))
  expect_true(all(grepl("CD[48]_T_cells", unlist(paths))))
  plots <- list()
  env$run_plots <- FALSE
  capture.output(eval(code, env))
  expect_length(plots, 0L)
})

test_that("ACS expression diagnostics back-transform only asinh-transformed channels", {
  skip_if_not_installed("UtilsCytoRSV")
  env <- .load_acs_method_env()
  source(file.path(dirname(script_helper), "acs_cytof-preprocess.R"), local = env)
  path <- tempfile("acs-gs-check-")
  withr::defer(unlink(path, recursive = TRUE))
  # Time is stored raw: 5 * sinh(700) is finite but far too large to bin,
  # and 5 * sinh(1000) is infinite.
  m <- cbind(Dy161Di = c(0, 1, 2, 3), Time = c(700, 800, 900, 1000))
  fs <- flowCore::flowSet(list(s1 = flowCore::flowFrame(m)))
  flowWorkspace::save_gs(flowWorkspace::GatingSet(flowWorkspace::flowSet_to_cytoset(fs)), path)
  plots <- env$.acsCytofPlotGatingSetCheck(path)
  expect_setequal(names(plots), c("asinh", "none"))
  none <- plots$none$data
  expect_equal(sort(none$expr[none$chnl == "Time"]), c(700, 800, 900, 1000))
  expect_equal(sort(none$expr[none$chnl == "Dy161Di"]), 5 * sinh(0:3))
  expect_equal(nrow(plots$asinh$data), 8L)
  for (p in plots) expect_no_error(ggplot2::ggplot_build(p))
})
