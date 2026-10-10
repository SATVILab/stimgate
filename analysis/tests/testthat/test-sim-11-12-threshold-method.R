.threshold_method_env <- function(files) {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."),
    winslash = "/")
  env <- new.env(parent = getNamespace("stimgate"))
  for (f in c("analysis-runtime.R", "analysis-plot-style.R", files)) {
    source(file.path(root, "scripts", "r", f), local = env)
  }
  env
}

.threshold_method_qmd_chunk <- function(qmd, label) {
  start <- match(paste0("#| label: ", label), qmd)
  if (is.na(start)) stop("Missing chunk: ", label)
  end <- start + which(qmd[(start + 1L):length(qmd)] == "```")[[1L]]
  qmd[(start + 1L):(end - 1L)]
}

.threshold_method_read_qmd <- function(name) {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."),
    winslash = "/")
  readLines(file.path(root, "analysis", name), warn = FALSE)
}

# A project directory holding only the saved channel settings.
.threshold_method_project <- function(methods) {
  path <- withr::local_tempdir(.local_envir = parent.frame())
  dir.create(file.path(path, "metaData"))
  settings <- lapply(methods, function(m) list(locThresholdMethod = m))
  names(settings) <- paste0("M", seq_along(methods))
  saveRDS(settings, file.path(path, "metaData", "chnlSettings.rds"))
  path
}

test_that("saved threshold methods must equal the requested method", {
  env <- .threshold_method_env(c("sim-low-separation.R", "sim-cluster-lab.R",
    "sim-compare-freq_bs.R", "sim-cluster-weak.R"))
  helpers <- c(".simLowSepThresholdMethod", ".simClusterThresholdMethod",
    ".simClusterWeakThresholdMethod")
  ok <- .threshold_method_project(c("region", "region"))
  mixed <- .threshold_method_project(c("region", "match"))
  missing <- .threshold_method_project(NA_character_)
  for (helper in helpers) {
    expect_identical(env[[helper]](ok, "region"), "region", info = helper)
    expect_error(env[[helper]](ok, "match"), "locThresholdMethod",
      info = helper)
    expect_error(env[[helper]](mixed, "region"), "locThresholdMethod",
      info = helper)
    expect_error(env[[helper]](missing, "region"), "locThresholdMethod",
      info = helper)
  }
})

test_that("Analysis 11 forwards the threshold method to stimControl()", {
  testthat::skip_if_not_installed("simcyto")
  env <- .threshold_method_env("sim-low-separation.R")
  withr::local_preserve_seed()
  seen <- NULL
  env$gateStim <- function(..., control) {
    seen <<- control$locThresholdMethod
    stop("captured threshold control")
  }
  row <- env$.simLowSepGrid(55700L)[1, ]
  row$n_cell <- 200
  for (method in c("region", "match")) {
    seen <- NULL
    settings <- env$.simLowSepMainSettings(loc_threshold_method = method)
    expect_identical(settings$loc_threshold_method, method)
    expect_error(env$.simLowSepRunScenario(row, 1L, settings),
      "captured threshold control")
    expect_identical(seen, method)
  }
  settings$loc_threshold_method <- NULL
  expect_error(env$.simLowSepRunScenario(row, 1L, settings),
    "set explicitly")
})

test_that("Analysis 11 rejects caches without the threshold method", {
  env <- .threshold_method_env("sim-low-separation.R")
  withr::local_envvar(ANALYSIS_EXPECTED_RUN_ID = NA)
  path <- withr::local_tempfile(fileext = ".rds")
  grid <- env$.simLowSepGrid(1L)
  main <- env$.simLowSepMainSettings(loc_threshold_method = "region")
  settings <- env$.simLowSepCacheSettings(grid, 2L, main, 1L,
    list(sim_size = "final"))
  expect_identical(settings$analysis_semantics_version,
    "sim-low-separation-v7")
  expect_identical(settings$settings$loc_threshold_method, "region")

  # A cache made before the change: v1, no method in settings or gates.
  old_main <- main[setdiff(names(main), "loc_threshold_method")]
  old <- env$.simLowSepCacheSettings(grid, 2L, old_main, 1L,
    list(sim_size = "final"))
  old$analysis_semantics_version <- "sim-low-separation-v1"
  env$.simLowSepWriteCache(list(gates = tibble::tibble(gate = 1)), old, path)
  expect_error(env$.simLowSepReadCache(path, settings, c("sim", "test")),
    "different settings")

  # Matching settings but gate rows without the method are also rejected.
  env$.simLowSepWriteCache(list(gates = tibble::tibble(gate = 1)), settings,
    path)
  expect_error(env$.simLowSepReadCache(path, settings, c("sim", "test")),
    "locThresholdMethod")
  env$.simLowSepWriteCache(
    list(gates = tibble::tibble(gate = 1, locThresholdMethod = "match")),
    settings, path
  )
  expect_error(env$.simLowSepReadCache(path, settings, c("sim", "test")),
    "locThresholdMethod")
})

test_that("QMD 11 sets the threshold method explicitly to cap", {
  qmd <- .threshold_method_read_qmd("11-sim-low-separation-cyt-pos.qmd")
  code <- paste(.threshold_method_qmd_chunk(qmd, "main-settings"),
    collapse = "\n")
  expect_match(code,
    '.simLowSepMainSettings(loc_threshold_method = "cap")', fixed = TRUE)
  cache <- paste(.threshold_method_qmd_chunk(qmd, "simulate"), collapse = "\n")
  expect_match(cache, "settings = main_settings", fixed = TRUE)
})

test_that("Analysis 12 wrappers forward the threshold method", {
  testthat::skip_if_not_installed("simcyto")
  env <- .threshold_method_env(c("sim-cluster-lab.R",
    "sim-compare-freq_bs.R", "sim-cluster-weak.R"))
  withr::local_preserve_seed()
  seen <- NULL
  capture_gate <- function(..., control) {
    seen <<- control$locThresholdMethod
    stop("captured threshold control")
  }
  testthat::local_mocked_bindings(gateStim = capture_gate,
    .package = "stimgate")
  weak <- list(seed = 559L, n_strong = 1L, n_weak = 1L,
    strong_prob = 0.01, weak_prob = 0.002, n_cell = 500L,
    mean_pos = 4.5, variance = 1.5, bw = 0.25, bias_uns = 0.25)
  for (method in c("region", "match")) {
    seen <- NULL
    expect_error(env$.simClusterLabRun(558L, 2L, 50L, 4,
      loc_threshold_method = method), "captured threshold control")
    expect_identical(seen, method)
    seen <- NULL
    expect_error(env$.simClusterWeakRun(utils::modifyList(weak,
      list(loc_threshold_method = method))), "captured threshold control")
    expect_identical(seen, method)
  }
  expect_identical(formals(env$.simClusterLabRun)$loc_threshold_method,
    "region")
  expect_error(env$.simClusterWeakRun(weak), "set explicitly")
})

test_that("Analysis 12 weak-response validation rejects old caches", {
  env <- .threshold_method_env(c("sim-compare-freq_bs.R",
    "sim-cluster-weak.R"))
  settings <- list(seed = 559L, n_strong = 1L, n_weak = 1L,
    strong_prob = 0.01, weak_prob = 0.002, n_cell = 5L,
    mean_pos = 4.5, variance = 1.5, bw = 0.25, bias_uns = 0.25,
    loc_threshold_method = "region")
  old_settings <- settings[setdiff(names(settings), "loc_threshold_method")]
  scores <- tibble::tibble(ind = "2", stage = c("Before", "After"))
  expect_identical(env$.simClusterWeakSemantics, "cluster-weak-v6")
  expect_error(env$.simClusterWeakValidate(list(semantics = "cluster-weak-v1",
    settings = old_settings, scores = scores), settings), "settings changed")
  expect_error(env$.simClusterWeakValidate(list(semantics = "cluster-weak-v6",
    settings = settings, scores = scores), settings), "locThresholdMethod")
  scores$locThresholdMethod <- "match"
  expect_error(env$.simClusterWeakValidate(list(semantics = "cluster-weak-v6",
    settings = settings, scores = scores), settings), "locThresholdMethod")
})

test_that("QMD 12 sets the method to cap and rejects old lab caches", {
  env <- .threshold_method_env("sim-cluster-lab.R")
  qmd <- .threshold_method_read_qmd("12-sim-cluster-gates.qmd")
  weak_code <- paste(.threshold_method_qmd_chunk(qmd, "weak-settings"),
    collapse = "\n")
  expect_match(weak_code, 'loc_threshold_method = "cap"', fixed = TRUE)

  dir_cache <- withr::local_tempdir()
  env$analysis_quick <- FALSE
  env$analysis_dev <- FALSE
  env$analysis_key <- c("sim", "cluster-gates")
  env$analysis_qmd <- "analysis/12-sim-cluster-gates.qmd"
  env$.analysis_cache_dir <- function(...) dir_cache
  eval(parse(text = .threshold_method_qmd_chunk(qmd, "lab-settings")),
    envir = env)
  expect_identical(env$lab_settings$loc_threshold_method, "cap")
  expect_identical(env$lab_semantics, "cluster-lab-v6")
  # The lab settings are passed to .simClusterLabRun() through do.call().
  expect_true(all(names(env$lab_settings) %in%
    names(formals(env$.simClusterLabRun))))

  read_code <- parse(text = .threshold_method_qmd_chunk(qmd, "lab-read"))
  old_settings <- env$lab_settings[
    setdiff(names(env$lab_settings), "loc_threshold_method")
  ]
  result <- list(allocations = tibble::tibble(ind = "2"))
  saveRDS(list(semantics = "cluster-lab-v1", settings = old_settings,
    result = result), env$lab_cache_file)
  expect_error(eval(read_code, envir = env), "different scientific settings")

  saveRDS(list(semantics = "cluster-lab-v6", settings = env$lab_settings,
    result = result), env$lab_cache_file)
  expect_error(eval(read_code, envir = env), "different scientific settings")

  result$allocations$locThresholdMethod <- "cap"
  saveRDS(list(semantics = "cluster-lab-v6", settings = env$lab_settings,
    result = result), env$lab_cache_file)
  expect_no_error(eval(read_code, envir = env))
  expect_identical(env$lab_result, result)
})
