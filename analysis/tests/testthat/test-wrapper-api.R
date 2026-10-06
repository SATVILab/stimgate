root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_cyt <- file.path(root_dir, "scripts", "r", "functionsForBenchmarking-Cyt.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_bw_io <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")
script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

.load_analysis_env <- function() {
  for (f in c(script_misc, script_cyt, script_bw, script_bw_io, script_bw_plot, script_comp)) {
    if (!file.exists(f)) stop("Expected analysis helper not found: ", f)
  }
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_cyt, local = env)
  source(script_bw, local = env)
  source(script_bw_io, local = env)
  source(script_bw_plot, local = env)
  source(script_comp, local = env)
  env
}

test_that("analysis wrapper functions do not expose removed calcSinglePosGates argument", {
  env <- .load_analysis_env()

  expect_false("calcSinglePosGates" %in% names(formals(env$.simBandwidthBsFreq)))
  expect_false("calcSinglePosGates" %in% names(formals(env$.simCompareStimgateRows)))
  expect_false("calcSinglePosGates" %in% names(formals(env$.simCompareFreqBs)))
})

test_that(".simBandwidthBsFreq exposes adaptive bandwidth settings directly", {
  env <- .load_analysis_env()

  adaptive_nm <- c(
    "bwAdaptive",
    "bwAdaptiveCore",
    "bwAdaptiveExtra",
    "bwAdaptiveCrossover",
    "bwAdaptiveTransitionWidth"
  )

  expect_true(all(adaptive_nm %in% names(formals(env$.simBandwidthBsFreq))))
})

test_that("analysis calls use the gateStim and stimControl argument contracts", {
  gate_args <- names(formals(stimgate::gateStim))
  control_args <- names(formals(stimgate::stimControl))
  api_args <- union(gate_args, control_args)
  removed_args <- c(paste0("tol", "Clust"), "gateQuant", "maxPosProbX")
  expect_false(any(removed_args %in% api_args))

  check_calls <- function(node) {
    if (!is.call(node) && !is.expression(node)) return(invisible(NULL))
    if (is.call(node)) {
      target <- paste(deparse(node[[1L]]), collapse = "")
      if (target %in% c("gateStim", "stimgate::gateStim", "stimgate::stimControl")) {
        forwarded <- names(as.list(node)[-1L])
        expect_true(all(forwarded %in% api_args), info = target)
        allowed <- if (target == "stimgate::stimControl") control_args else gate_args
        expect_true(all(forwarded %in% allowed), info = target)
        if (target == "stimgate::stimControl" && "clusterGates" %in% forwarded) {
          expect_true(
            identical(as.list(node)$clusterGates, quote(clusterGates)) ||
              identical(as.list(node)$clusterGates, FALSE) ||
              identical(as.list(node)$clusterGates, TRUE)
          )
        }
      }
    }
    for (i in seq_along(node)) {
      if (identical(node[[i]], quote(expr = ))) next
      if (is.call(node[[i]]) || is.expression(node[[i]])) check_calls(node[[i]])
    }
    invisible(NULL)
  }

  scripts <- list.files(file.path(root_dir, "scripts", "r"), pattern = "\\.R$", full.names = TRUE)
  for (script in scripts) check_calls(parse(script))
})

test_that("scripts pass bw to gateStim and calcCytPosGates/minCell to stimControl", {
  gate_args <- names(formals(stimgate::gateStim))
  control_args <- names(formals(stimgate::stimControl))
  expect_true("bw" %in% gate_args)
  expect_false(any(c("calcCytPosGates", "minCell") %in% gate_args))
  expect_true(all(c("calcCytPosGates", "minCell") %in% control_args))
  expect_false("bw" %in% control_args)

  check_split <- function(node) {
    if (!is.call(node) && !is.expression(node)) return(invisible(NULL))
    if (is.call(node)) {
      target <- paste(deparse(node[[1L]]), collapse = "")
      forwarded <- names(as.list(node)[-1L])
      if (target %in% c("gateStim", "stimgate::gateStim")) {
        expect_false(
          any(c("calcCytPosGates", "minCell") %in% forwarded),
          info = target
        )
      }
      if (target == "stimgate::stimControl") {
        expect_false("bw" %in% forwarded, info = target)
      }
    }
    for (i in seq_along(node)) {
      if (identical(node[[i]], quote(expr = ))) next
      if (is.call(node[[i]]) || is.expression(node[[i]])) check_split(node[[i]])
    }
    invisible(NULL)
  }

  scripts <- list.files(
    file.path(root_dir, "scripts", "r"),
    pattern = "\\.R$",
    full.names = TRUE
  )
  for (script in scripts) check_split(parse(script))
})

test_that("analysis wrappers retain legacy arguments for callers and manifests", {
  env <- .load_analysis_env()
  legacy_args <- c("gateQuant", "maxPosProbX")
  for (wrapper in c(".simBandwidthBsFreq", ".simCompareStimgateRows", ".simCompareFreqBs")) {
    expect_true(all(legacy_args %in% names(formals(env[[wrapper]]))), info = wrapper)
  }
})

test_that(".simBandwidthBsFreq forwards to gateStim without unknown-argument error", {
  env <- .load_analysis_env()

  res <- env$.simBandwidthBsFreq(
    nSample = 1L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0,
    bwAdaptive = TRUE,
    bwAdaptiveCore = 0.1,
    bwAdaptiveExtra = 0.25,
    bwAdaptiveCrossover = NULL,
    bwAdaptiveTransitionWidth = 0,
    bwFallback = 0.2,
    bwMin = "none",
    bwMax = "none",
    bwMtd = "hpi1",
    nCellStim = 50L,
    probResponse = 0.1,
    meanPos = 5,
    transformation = "gaussian",
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.1,
    ncellUnsRelativeToStim = 1
  )

  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
  expect_type(res$clusterGates, "logical")
  expect_true(all(!res$clusterGates))
})

test_that(".simCompareFreqBs forwards to gateStim via .simCompareStimgateRows without unknown-argument error", {
  env <- .load_analysis_env()

  res <- env$.simCompareFreqBs(
    nSample = 1L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0,
    bw = 0.1,
    bwMtd = "hpi1",
    nCellStim = 50L,
    probResponse = 0.1,
    meanPos = 5,
    transformation = "gaussian",
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.1,
    ncellUnsRelativeToStim = 1
  )

  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
  expect_type(res$clusterGates, "logical")
  expect_true(all(!res$clusterGates))

  stimgate_res <- res[res[["approach"]] == "stimgate", , drop = FALSE]
  expect_true(nrow(stimgate_res) > 0)
  expect_false(any(stimgate_res[["method"]] == "stimgate_error"))
  expect_true(all(is.na(stimgate_res[["error"]])))
})

test_that(
  ".simCompareFreqBs and .simCompareStimgateRows use the stimControl gateCombn default",
  {
    env <- .load_analysis_env()

    expect_false("gateCombn" %in% names(formals(env$.simCompareStimgateRows)))
    expect_false("gateCombn" %in% names(formals(env$.simCompareFreqBs)))
  }
)


test_that("threshold clustering accepts only one non-missing logical value", {
  env <- .load_analysis_env()
  wrappers <- c(
    ".simBandwidthBsFreq", ".simCompareStimgateRows", ".simCompareFreqBs",
    ".simCompareRunScenario", ".simCompareFreqBsGrid"
  )
  for (wrapper in wrappers) {
    expect_identical(formals(env[[wrapper]])$clusterGates, FALSE)
    for (value in list(NULL, logical(), NA, NA_real_, 0, 1e-7, "TRUE", c(TRUE, FALSE))) {
      expect_error(
        do.call(env[[wrapper]], list(clusterGates = value)),
        "clusterGates must be a single TRUE or FALSE", info = wrapper
      )
    }
  }
  for (wrapper in c(".simBandwidthEstBwDirect", ".simBandwidthEstBwDirectAdaptive")) {
    expect_false("clusterGates" %in% names(formals(env[[wrapper]])))
  }
})

test_that("gating wrappers pass the logical clustering toggle to stimControl", {
  skip_if_not_installed("simcyto")
  env <- .load_analysis_env()
  seen <- NULL
  capture_gate <- function(..., control) {
    seen <<- control$clusterGates
    stop("captured clustering control")
  }
  env$gateStim <- capture_gate
  testthat::local_mocked_bindings(gateStim = capture_gate, .package = "stimgate")
  env$.simCompareTruthTable <- function(...) tibble::tibble()
  for (toggle in c(FALSE, TRUE)) {
    seen <- NULL
    expect_error(env$.simBandwidthBsFreq(
      nSample = 1L, nMarker = 1L, nCondition = 2L, nCluster = 2L,
      nIter = 1L, biasUns = 0, bw = 0.1,
      nCellStim = 50L, probResponse = 0.1, meanPos = 5,
      transformation = "gaussian", samplePerturbationSd = 0,
      conditionPerturbationSd = 0, clusterPerturbationSd = 0,
      backgroundRelativeToResponse = 0.1, ncellUnsRelativeToStim = 1,
      clusterGates = toggle
    ), "captured clustering control")
    expect_identical(seen, toggle)

    seen <- NULL
    env$.simCompareStimgateRows(
      gs = NULL, labelsList = NULL, pathProject = tempdir(),
      nSample = 1L, nCondition = 2L, nMarker = 1L,
      biasUns = 0, bw = 0.1, clusterGates = toggle
    )
    expect_identical(seen, toggle)
  }
})

test_that("analysis code no longer uses the retired clustering tolerance names", {
  retired <- paste0("tol", c("Clust", "_clust"))
  skip_if(Sys.which("git") == "", "Source audit requires git")
  # Render intermediates can retain retired API names; audit source files only.
  source_paths <- system2("git", c(
    "-C", shQuote(root_dir), "ls-files", "--cached", "--others",
    "--exclude-standard", "--", "analysis", "scripts/r"
  ), stdout = TRUE)
  expect_null(attr(source_paths, "status"))
  files <- file.path(root_dir, unique(source_paths[
    grepl("\\.(R|qmd|md|sh|ya?ml)$", source_paths)
  ]))
  expect_gt(length(files), 0L)
  for (file in files) {
    text <- readLines(file, warn = FALSE)
    expect_false(any(vapply(retired, function(name) {
      any(grepl(name, text, fixed = TRUE))
    }, logical(1))), info = file)
  }
})
