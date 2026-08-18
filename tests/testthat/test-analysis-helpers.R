test_that("Quarto analysis files do not call scripts/r helpers through stimgate:::", {
  root_dir <- testthat::test_path("../..")
  analysis_dir <- file.path(root_dir, "analysis")
  skip_if_not(dir.exists(analysis_dir), "analysis directory not found")

  qmd_files <- list.files(
    analysis_dir,
    pattern = "\\.qmd$",
    recursive = TRUE,
    full.names = TRUE
  )
  expect_true(length(qmd_files) > 0)

  scripts_r_dir <- file.path(root_dir, "scripts", "r")
  script_files <- list.files(
    scripts_r_dir,
    pattern = "\\.R$",
    recursive = TRUE,
    full.names = TRUE
  )

  script_fns <- character(0)
  for (f in script_files) {
    lines <- readLines(f, warn = FALSE)
    fn_matches <- grep("^\\s*(\\.?[A-Za-z0-9_]+)\\s*<-\\s*function", lines, value = TRUE)
    extracted <- sub("^\\s*(\\.?[A-Za-z0-9_]+)\\s*<-\\s*function.*$", "\\1", fn_matches)
    script_fns <- c(script_fns, extracted)
  }
  script_fns <- unique(script_fns)

  violations <- character(0)
  for (qmd_path in qmd_files) {
    qmd_rel <- file.path("analysis", sub(".*[/\\]analysis[/\\]", "", qmd_path))
    lines <- readLines(qmd_path, warn = FALSE)
    for (i in seq_along(lines)) {
      line <- lines[i]
      if (grepl("stimgate:::", line)) {
        for (fn_name in script_fns) {
          pattern <- paste0("stimgate:::", gsub(".", "\\.", fn_name, fixed = TRUE))
          if (grepl(pattern, line)) {
            violations <- c(
              violations,
              sprintf("%s (line %d): %s", qmd_rel, i, trimws(line))
            )
          }
        }
        if (grepl("stimgate:::\\.sim", line)) {
          violations <- c(
            violations,
            sprintf("%s (line %d): %s", qmd_rel, i, trimws(line))
          )
        }
      }
    }
  }

  expect_equal(
    unique(violations),
    character(0),
    info = "Found QMD files calling scripts/r helpers via stimgate:::"
  )
})

test_that("scripts/r helper sets source cleanly and resolve required dependencies", {
  root_dir <- testthat::test_path("../..")
  scripts_r_dir <- file.path(root_dir, "scripts", "r")
  skip_if_not(dir.exists(scripts_r_dir), "scripts/r directory not found")

  helper_files <- c(
    "sim-misc.R",
    "sim-trans.R",
    "sim-bandwidth.R",
    "sim-compare-freq_bs.R",
    "functionsForBenchmarking-Pheno.R"
  )

  test_env <- new.env(parent = globalenv())

  for (f in helper_files) {
    f_path <- file.path(scripts_r_dir, f)
    skip_if_not(file.exists(f_path), paste(f, "not found"))
    expect_silent(sys.source(f_path, envir = test_env))
  }

  # Verify core functions exist in test_env
  expect_true(exists(".simMiscGetTrans", envir = test_env, mode = "function"))
  expect_true(exists(".simBandwidthReadLocDetails", envir = test_env, mode = "function"))
  expect_true(exists(".simBandwidthGetTrans", envir = test_env, mode = "function"))
  expect_true(exists(".simCompareReadLocDetails", envir = test_env, mode = "function"))
  expect_true(exists(".simCompareGetTrans", envir = test_env, mode = "function"))

  # Test dependency resolution across helpers
  # .simCompareGetTrans should find .simBandwidthGetTrans or .simMiscGetTrans
  get_trans_fn <- get(".simCompareGetTrans", envir = test_env)
  expect_no_error(get_trans_fn("gaussian"))

  # .simCompareReadLocDetails should delegate to .simBandwidthReadLocDetails
  read_details_fn <- get(".simCompareReadLocDetails", envir = test_env)
  dummy_dir <- file.path(tempdir(), "nonexistent_proj")
  expect_equal(read_details_fn(dummy_dir, nSample = 1, nCondition = 2), tibble::tibble())
})

test_that("wrapper functions forward parameters that exist in gateStim formals", {
  skip_if_not_installed("stimgate")

  root_dir <- testthat::test_path("../..")
  script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
  script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
  script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

  skip_if_not(file.exists(script_bw), "sim-bandwidth.R not found")

  source(script_misc, local = TRUE)
  source(script_bw, local = TRUE)
  source(script_comp, local = TRUE)

  gate_stim_formals <- names(formals(stimgate::gateStim))

  wrappers <- list(
    .simBandwidthBsFreq = .simBandwidthBsFreq,
    .simBandwidthEstBw = .simBandwidthEstBw,
    .simPlotGate = .simPlotGate,
    .simCompareStimgateRows = .simCompareStimgateRows
  )

  extract_gate_stim_args <- function(fn) {
    fn_body <- body(fn)
    found_args <- character(0)

    find_calls <- function(expr) {
      if (is.call(expr)) {
        call_name <- as.character(expr[[1]])
        if (length(call_name) > 0 && call_name[length(call_name)] == "gateStim") {
          call_args <- names(expr)[-1]
          call_args <- call_args[nzchar(call_args)]
          found_args <<- c(found_args, call_args)
        }
        for (i in seq_along(expr)) {
          find_calls(expr[[i]])
        }
      }
    }

    find_calls(fn_body)
    unique(found_args)
  }

  for (wrap_name in names(wrappers)) {
    forwarded <- extract_gate_stim_args(wrappers[[wrap_name]])
    expect_true(
      length(forwarded) > 0,
      info = paste("No gateStim call found in", wrap_name)
    )

    missing_params <- setdiff(forwarded, gate_stim_formals)
    expect_equal(
      missing_params,
      character(0),
      info = sprintf(
        "Wrapper %s forwards parameters not present in gateStim formals: %s",
        wrap_name,
        paste(missing_params, collapse = ", ")
      )
    )
  }
})

test_that("analysis helpers do not expose or forward calcSinglePosGates", {
  skip_if_not_installed("stimgate")

  root_dir <- testthat::test_path("../..")
  script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
  script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
  script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

  skip_if_not(file.exists(script_bw), "sim-bandwidth.R not found")

  source(script_misc, local = TRUE)
  source(script_bw, local = TRUE)
  source(script_comp, local = TRUE)

  expect_false("calcSinglePosGates" %in% names(formals(.simBandwidthBsFreq)))
  expect_false("calcSinglePosGates" %in% names(formals(.simCompareStimgateRows)))
  expect_false("calcSinglePosGates" %in% names(formals(.simCompareFreqBs)))
})

test_that("simulation wrappers execute tiny representative calls without errors", {
  skip_if_not_installed("stimgate")

  root_dir <- testthat::test_path("../..")
  script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
  script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")

  skip_if_not(file.exists(script_bw), "sim-bandwidth.R not found")

  source(script_misc, local = TRUE)
  source(script_bw, local = TRUE)

  # Tiny call for .simBandwidthBsFreq
  res_bs <- .simBandwidthBsFreq(
    nSample = 1L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0,
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

  expect_s3_class(res_bs, "data.frame")
  expect_true(nrow(res_bs) > 0)

  # Tiny call for .simBandwidthEstBw
  res_est <- .simBandwidthEstBw(
    nSample = 1L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0,
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

  expect_s3_class(res_est, "data.frame")
  expect_true(nrow(res_est) > 0)
})

test_that("direct bandwidth simulation helpers agree numerically with package implementation", {
  skip_if_not_installed("stimgate")

  root_dir <- testthat::test_path("../..")
  script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
  script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")

  skip_if_not(file.exists(script_bw), "sim-bandwidth.R not found")

  source(script_misc, local = TRUE)
  source(script_bw, local = TRUE)

  set.seed(42)
  x_sample <- rnorm(200, mean = 5, sd = 1.5)

  # Check standard normalized bandwidth via .simBandwidthBwOne vs .bwCalcOne
  bw_direct_std <- .simBandwidthBwOne(
    x = x_sample,
    bwMtd = "hpi1Norm",
    bwMin = 1e-10,
    bwMax = 1e10,
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L,
    bwFallback = NULL,
    normPeakFrac = 0.15,
    normExtraFrac = 0.25
  )

  bw_pkg_std <- stimgate:::.bwCalcOne(
    x = x_sample,
    bwMtd = "hpi1Norm",
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L,
    normPeakFrac = 0.15,
    normExtraFrac = 0.25,
    adaptive = FALSE
  )

  expect_equal(as.numeric(bw_direct_std), as.numeric(bw_pkg_std))
})
