root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

.compare_method_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "sim-misc.R", "sim-bandwidth.R",
    "sim-compare-freq_bs.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  env
}

.compare_method_qmds <- c(
  "7-sim-compare-freq_bs.qmd" = "corrected-comparison-v22",
  "8-sim-compare-freq_bs-batch.qmd" = "batch-mismatch-comparison-v22"
)

test_that("comparison QMDs set and record the threshold method", {
  for (qmd in names(.compare_method_qmds)) {
    lines <- readLines(file.path(root_dir, "analysis", qmd), warn = FALSE)
    txt <- paste(lines, collapse = "\n")
    expect_match(
      txt, 'stimgate_loc_threshold_method <- "cap"',
      fixed = TRUE, info = qmd
    )
    expect_match(
      txt,
      paste0(
        'comparison_semantics_version <- "', .compare_method_qmds[[qmd]], '"'
      ),
      fixed = TRUE, info = qmd
    )
    # Every StimGate wrapper call passes the method explicitly.
    n_calls <- sum(grepl(
      "locEnforceShapeThreshold = loc_enforce_shape_threshold,", lines,
      fixed = TRUE
    ))
    n_method <- sum(grepl(
      "locThresholdMethod = stimgate_loc_threshold_method,", lines,
      fixed = TRUE
    ))
    expect_gt(n_calls, 0L)
    expect_identical(n_method, n_calls, info = qmd)

    # The method is part of the settings recorded in the manifest and
    # required by canonical reads.
    start <- which(grepl("^analysis_result_params <- list", lines))
    expect_length(start, 1L)
    end <- start + which(lines[(start + 1L):length(lines)] == ")")[[1L]]
    expr <- parse(text = lines[start:end])[[1L]][[3L]]
    expect_identical(
      as.list(expr)$stimgate_loc_threshold_method,
      quote(stimgate_loc_threshold_method),
      info = qmd
    )
    expect_match(
      txt, "required_params = analysis_result_params",
      fixed = TRUE, info = qmd
    )
  }
})

test_that("canonical comparison results without the method are rejected", {
  env <- .compare_method_env()
  current <- withr::local_tempdir()
  file.create(file.path(current, "COMPLETE"))
  saveRDS(1L, file.path(current, "result.rds"))
  ctx <- list(
    analysis_key = c("sim", "compare", "freq_bs"), current_dir = current,
    qmd_path = "analysis/7-sim-compare-freq_bs.qmd"
  )
  withr::local_envvar(ANALYSIS_EXPECTED_RUN_ID = NA_character_)
  old <- list(
    comparison_semantics_version = "corrected-comparison-v17",
    cluster_gates = TRUE
  )
  saveRDS(
    list(analysis_key = ctx$analysis_key, params = old),
    file.path(current, "manifest.rds")
  )
  required <- list(
    comparison_semantics_version = "corrected-comparison-v22",
    cluster_gates = TRUE,
    stimgate_loc_threshold_method = "region"
  )
  expect_error(
    env$.analysis_current_file(ctx, "result.rds", required),
    "comparison_semantics_version, stimgate_loc_threshold_method",
    fixed = TRUE
  )
  # Only the version bump: the missing method alone still rejects the cache.
  old$comparison_semantics_version <- "corrected-comparison-v22"
  saveRDS(
    list(analysis_key = ctx$analysis_key, params = old),
    file.path(current, "manifest.rds")
  )
  expect_error(
    env$.analysis_current_file(ctx, "result.rds", required),
    "parameters: stimgate_loc_threshold_method", fixed = TRUE
  )
  saveRDS(
    list(analysis_key = ctx$analysis_key, params = required),
    file.path(current, "manifest.rds")
  )
  expect_true(file.exists(
    env$.analysis_current_file(ctx, "result.rds", required)
  ))
})

test_that("comparison wrappers default to and forward the region method", {
  env <- .compare_method_env()
  for (wrapper in c(
    ".simCompareStimgateRows", ".simCompareFreqBs",
    ".simCompareRunScenario", ".simCompareFreqBsGrid"
  )) {
    expect_identical(
      formals(env[[wrapper]])$locThresholdMethod, "region",
      info = wrapper
    )
  }

  # .simCompareStimgateRows -> stimControl()
  seen <- NULL
  capture_gate <- function(..., control) {
    seen <<- control$locThresholdMethod
    stop("captured threshold method")
  }
  testthat::local_mocked_bindings(
    gateStim = capture_gate, .package = "stimgate"
  )
  env$.simCompareTruthTable <- function(...) {
    tibble::tibble(sample = "1", ind = "2", chnl = "F1")
  }
  for (method in c("region", "match")) {
    seen <- NULL
    rows <- env$.simCompareStimgateRows(
      gs = NULL, labelsList = NULL, pathProject = withr::local_tempdir(),
      nSample = 1L, nCondition = 2L, nMarker = 1L,
      biasUns = 0, bw = 0.1, locThresholdMethod = method
    )
    expect_identical(seen, method)
    # Failed StimGate rows still say which method was requested.
    expect_identical(rows$method, "stimgate_error")
    expect_identical(rows$locThresholdMethod, method)
  }

  # .simCompareFreqBsGrid -> .simCompareRunScenario -> .simCompareFreqBs
  env$.simCompareEnsureCurrentCheckout <- function(...) invisible(TRUE)
  env$.simCompareFreqBs <- function(...) {
    tibble::tibble(
      approach = "stimgate", method = "stimgate",
      locThresholdMethod = list(...)$locThresholdMethod,
      error = NA_character_
    )
  }
  row <- tibble::tibble(sim_id = 1L, sim_seed = 11L)
  out <- env$.simCompareRunScenario(row, nSample = 1, nIter = 1)
  expect_identical(out$locThresholdMethod, "region")
  out <- env$.simCompareRunScenario(
    row, nSample = 1, nIter = 1, locThresholdMethod = "match"
  )
  expect_identical(out$locThresholdMethod, "match")
  out <- env$.simCompareFreqBsGrid(
    row, nSample = 1, nIter = 1, parallel = FALSE, progress = FALSE,
    locThresholdMethod = "match"
  )
  expect_identical(out$locThresholdMethod, "match")
})

test_that("StimGate rows take the method from the saved channel settings", {
  env <- .compare_method_env()
  path <- withr::local_tempdir()
  dir.create(file.path(path, "metaData"))
  save_settings <- function(settings) {
    saveRDS(settings, file.path(path, "metaData", "chnlSettings.rds"))
  }
  save_settings(list(
    MarkerF1 = list(chnlCut = "F1", locThresholdMethod = "match")
  ))
  expect_identical(
    env$.simCompareStimgateLocThresholdMethod(path, "MarkerF1"), "match"
  )
  save_settings(list(
    Other = list(chnlCut = "F1", locThresholdMethod = "region")
  ))
  expect_identical(
    env$.simCompareStimgateLocThresholdMethod(path, "MarkerF1"), "region"
  )
  save_settings(list(MarkerF1 = list(chnlCut = "F1")))
  expect_error(
    env$.simCompareStimgateLocThresholdMethod(path, "MarkerF1"),
    "no valid locThresholdMethod"
  )
})

test_that("scenario caches without the requested method are not reused", {
  env <- .compare_method_env()
  cached <- tibble::tibble(
    sim_id = 1L, approach = c("stimgate", "stimgate", "fbeta"),
    method = c("stimgate", "stimgate_loc_condition", "fbeta"),
    locThresholdMethod = c("region", "region", NA_character_),
    error = NA_character_
  )
  row <- data.frame(sim_id = 1L)
  validate <- function(x, method) {
    env$.simCompareValidateScenarioCache(x, row, locThresholdMethod = method)
  }
  expect_true(validate(cached, "region"))
  expect_false(validate(cached, "match"))
  old <- cached[, setdiff(names(cached), "locThresholdMethod")]
  expect_false(validate(old, "region"))
  partial <- cached
  partial$locThresholdMethod[[2L]] <- NA_character_
  expect_false(validate(partial, "region"))
  # Callers that do not request a method keep the previous behaviour.
  expect_true(validate(old, NULL))

  # Resume reruns a scenario whose cache predates the recorded method.
  dir_cache <- withr::local_tempdir()
  env$.simCompareEnsureCurrentCheckout <- function(...) invisible(TRUE)
  calls <- 0L
  env$.simCompareFreqBs <- function(...) {
    calls <<- calls + 1L
    tibble::tibble(
      approach = "stimgate", method = "stimgate",
      locThresholdMethod = list(...)$locThresholdMethod,
      error = NA_character_
    )
  }
  # Completeness checks need nSample/nIter; leave them NULL here so only the
  # method decides whether the cache is reused.
  env$.simCompareValidateScenarioCache <- local({
    original <- env$.simCompareValidateScenarioCache
    function(cached, row, nSample = NULL, nIter = NULL, ...) {
      original(cached, row, ...)
    }
  })
  path_out <- env$.simCompareScenarioOutputPath(
    sim_id = 1L, dirCache = dir_cache
  )
  dir.create(dirname(path_out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(old, path_out)
  scenario_row <- tibble::tibble(sim_id = 1L)
  out <- env$.simCompareRunScenario(
    scenario_row, nSample = 1, nIter = 1, dirCache = dir_cache
  )
  expect_identical(calls, 1L)
  expect_identical(out$locThresholdMethod, "region")
  out <- env$.simCompareRunScenario(
    scenario_row, nSample = 1, nIter = 1, dirCache = dir_cache
  )
  expect_identical(calls, 1L)
  out <- env$.simCompareRunScenario(
    scenario_row, nSample = 1, nIter = 1, dirCache = dir_cache,
    locThresholdMethod = "match"
  )
  expect_identical(calls, 2L)
  expect_identical(out$locThresholdMethod, "match")
})

test_that("comparison StimGate rows record the method StimGate used", {
  testthat::skip_if_not_installed("simcyto")
  env <- .compare_method_env()
  withr::local_preserve_seed()
  set.seed(1)
  for (method in c("region", "match")) {
    res <- suppressWarnings(env$.simCompareFreqBs(
      nSample = 1L, nMarker = 1L, nCondition = 2L, nCluster = 2L,
      nIter = 1L, biasUns = 0, bw = 0.1, bwMtd = "hpi1",
      nCellStim = 50L, probResponse = 0.1, meanPos = 5,
      transformation = "gaussian", samplePerturbationSd = 0,
      conditionPerturbationSd = 0, clusterPerturbationSd = 0,
      backgroundRelativeToResponse = 0.1, ncellUnsRelativeToStim = 1,
      locThresholdMethod = method
    ))
    stimgate_rows <- res[res$approach %in% "stimgate", , drop = FALSE]
    expect_gt(nrow(stimgate_rows), 0L)
    expect_true(all(is.na(stimgate_rows$error)))
    expect_true(all(stimgate_rows$locThresholdMethod == method))
    expect_true(all(is.na(
      res$locThresholdMethod[!res$approach %in% "stimgate"]
    )))
  }
})
