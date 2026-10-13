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

test_that("Analysis 11 forwards the threshold and cytokine-positive settings", {
  testthat::skip_if_not_installed("simcyto")
  env <- .threshold_method_env("sim-low-separation.R")
  withr::local_preserve_seed()
  seen <- NULL
  env$gateStim <- function(..., calcCytPosGates, control) {
    seen <<- list(method = control$locThresholdMethod, plateau = control$locPlateauStep, cyt_pos = calcCytPosGates,
      cyt_pos_method = control$cytPosMethod)
    stop("captured threshold control")
  }
  env$stimControl <- function(...) list(...)
  row <- env$.simLowSepGrid(55700L)[1, ]
  row$n_cell <- 200
  for (method in c("region", "match")) {
    seen <- NULL
    settings <- env$.simLowSepMainSettings(loc_threshold_method = method, n_sample = 1L, locPlateauStep = TRUE)
    expect_identical(settings$loc_threshold_method, method)
    expect_error(env$.simLowSepScenario(row, settings), "captured threshold control")
    expect_identical(seen, list(method = method, plateau = TRUE, cyt_pos = TRUE, cyt_pos_method = "antimode"))
  }
  settings$loc_threshold_method <- NULL
  expect_error(env$.simLowSepScenario(row, settings), "set explicitly")
  settings <- env$.simLowSepMainSettings(loc_threshold_method = "cap", n_sample = 1L)
  settings$cyt_pos_method <- "coexpression"
  expect_error(env$.simLowSepScenario(row, settings), "antimode")
})

test_that("Analysis 11 validation rejects gates without the requested method", {
  env <- .threshold_method_env("sim-low-separation.R")
  gates <- function(method) {
    tibble::tibble(sample = 1L, marker = c("TNF", "IFNg"), gate = 1, gate_cyt = 1,
      locThresholdMethod = method)
  }
  tbl <- tibble::tibble(sim_id = 1:2, n_sample_run = 1L,
    gates = list(gates("cap"), gates("region")),
    coex_gates = list(tibble::tibble(x = 1:4), tibble::tibble(x = 1:4)),
    counts = list(tibble::tibble(x = 1:128), tibble::tibble(x = 1:128)),
    cells = list(NULL, NULL))
  problems <- env$.simLowSepValidator(1L, "cap")(tbl)
  expect_length(problems, 1L)
  expect_match(problems, "locThresholdMethod for sim_id: 2")
  tbl$gates[[2]]$locThresholdMethod <- NULL
  expect_match(env$.simLowSepValidator(1L, "cap")(tbl), "sim_id: 2")
})

test_that("QMDs 11a and 11b set the threshold method explicitly to cap", {
  for (name in c("11a-sim-low-separation-cyt-pos.qmd", "11b-sim-low-separation-batch.qmd")) {
    qmd <- .threshold_method_read_qmd(name)
    code <- paste(.threshold_method_qmd_chunk(qmd, "main-settings"),
      collapse = "\n")
    expect_match(code,
      'loc_threshold_method = "cap", n_sample = sizes$n_sample,',
      fixed = TRUE, info = name)
    expect_match(code, "loc_plateau_step <- TRUE", fixed = TRUE)
    expect_match(code, "locPlateauStep = loc_plateau_step", fixed = TRUE)
    grid <- paste(.threshold_method_qmd_chunk(qmd, "grid"), collapse = "\n")
    expect_match(grid, "scenario_settings = main_settings", fixed = TRUE, info = name)
    sim <- paste(.threshold_method_qmd_chunk(qmd, "simulate"), collapse = "\n")
    expect_match(sim, "settings = main_settings", fixed = TRUE, info = name)
    expect_match(sim, "main_settings$loc_threshold_method", fixed = TRUE, info = name)
  }
})

test_that("Analysis 12 wrappers forward the threshold method", {
  testthat::skip_if_not_installed("simcyto")
  env <- .threshold_method_env(c("sim-cluster-lab.R",
    "sim-compare-freq_bs.R", "sim-cluster-weak.R"))
  withr::local_preserve_seed()
  seen <- NULL
  capture_gate <- function(..., control) {
    seen <<- list(method = control$locThresholdMethod, plateau = control$locPlateauStep)
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
      loc_threshold_method = method, locPlateauStep = TRUE), "captured threshold control")
    expect_identical(seen, list(method = method, plateau = TRUE))
    seen <- NULL
    expect_error(env$.simClusterWeakRun(utils::modifyList(weak,
      list(loc_threshold_method = method, locPlateauStep = TRUE))), "captured threshold control")
    expect_identical(seen, list(method = method, plateau = TRUE))
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
  expect_identical(env$.simClusterWeakSemantics, "cluster-weak-v7")
  expect_error(env$.simClusterWeakValidate(list(semantics = "cluster-weak-v1",
    settings = old_settings, scores = scores), settings), "settings changed")
  expect_error(env$.simClusterWeakValidate(list(semantics = "cluster-weak-v7",
    settings = settings, scores = scores), settings), "locThresholdMethod")
  scores$locThresholdMethod <- "match"
  expect_error(env$.simClusterWeakValidate(list(semantics = "cluster-weak-v7",
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
  expect_true(env$lab_settings$locPlateauStep)
  expect_match(weak_code, "locPlateauStep = loc_plateau_step", fixed = TRUE)
  expect_identical(env$lab_semantics, "cluster-lab-v7")
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

  saveRDS(list(semantics = "cluster-lab-v7", settings = env$lab_settings,
    result = result), env$lab_cache_file)
  expect_error(eval(read_code, envir = env), "different scientific settings")

  result$allocations$locThresholdMethod <- "cap"
  result$allocations$locPlateauStep <- TRUE
  saveRDS(list(semantics = "cluster-lab-v7", settings = env$lab_settings,
    result = result), env$lab_cache_file)
  expect_no_error(eval(read_code, envir = env))
  expect_identical(env$lab_result, result)
})
