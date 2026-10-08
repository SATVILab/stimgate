root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)
scripts_r_dir <- file.path(root_dir, "scripts", "r")
qmd9_path <- file.path(root_dir, "analysis", "9-real-compare-acs-cytof.qmd")
qmd10_path <- file.path(
  root_dir, "analysis", "10-real-compare-acs-cytof-validation.qmd"
)

.load_acs_threshold_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c(
    "analysis-runtime.R", "analysis-plot-style.R",
    "sim-bandwidth-analysis-plot.R", "acs_cytof-helper.R",
    "acs_cytof-gate.R", "acs_cytof-methods.R", "acs_cytof-manual.R",
    "acs_cytof-plot_cyt.R"
  )) {
    source(file.path(scripts_r_dir, file), local = env)
  }
  env
}

.acs_stimgate_manifest <- function(method = "region", resolved = method) {
  list(
    context = list(gitSha = "abc", preprocessing = list(hash = "files")),
    settings = list(clusterGates = TRUE, locThresholdMethod = method),
    channelSettings = list(IFNg = list(locThresholdMethod = resolved))
  )
}

# Run the ACS population runner with gating mocked, returning the control
# passed to gateStim() and the population paths.
.run_acs_runner_mocked <- function(env, pathScratch, requested, resolved) {
  pathGs <- file.path(pathScratch, "gs")
  dir.create(file.path(pathGs, "cd4"), recursive = TRUE)
  env$.acsCytofEnsureCurrentCheckout <- function(...) invisible(TRUE)
  env$.acsCytofSetDebug <- function() function() invisible(NULL)
  env$.acsCytofManifest <- function(preprocessing) list(version = 1L)
  env$.acsCytofReadPreprocessing <- function(...) {
    list(sampleMap = data.frame(
      SampleID = "a", stim = c("uns", "p1", "mtb", "ebv", "p4"), ind = 1:5
    ))
  }
  seen <- new.env()
  testthat::local_mocked_bindings(
    load_gs = function(...) as.list(1:5), .package = "flowWorkspace"
  )
  testthat::local_mocked_bindings(
    gateStim = function(pathProject, ..., control) {
      seen$control <- control
      invisible(pathProject)
    },
    stimgateMetaReadSettingsChnls = function(...) {
      list(IFNg = list(locThresholdMethod = resolved))
    },
    .package = "stimgate"
  )
  args <- list(
    pop = "cd4", pathFcsBase = "raw", pathGsBase = pathGs,
    pathScratchBase = pathScratch, runPreprocessing = FALSE,
    runMethods = TRUE, runPlots = FALSE
  )
  if (!is.null(requested)) args$locThresholdMethod <- requested
  out <- do.call(env$.acsCytofRunPopulation, args)
  list(control = seen$control, paths = out$paths)
}

test_that("the ACS runner forwards and records locThresholdMethod", {
  env <- .load_acs_threshold_env()
  expect_identical(
    formals(env$.acsCytofRunPopulation)$locThresholdMethod, "region"
  )
  for (method in list(NULL, "region", "match")) {
    expected <- if (is.null(method)) "region" else method
    pathScratch <- tempfile("acs-threshold-")
    withr::defer(unlink(pathScratch, recursive = TRUE))
    run <- .run_acs_runner_mocked(env, pathScratch, method, expected)
    expect_s3_class(run$control, "stimControl")
    expect_identical(run$control$locThresholdMethod, expected)
    expect_true(run$control$clusterGates)
    expect_true(run$control$calcCytPosGates)
    manifest <- readRDS(file.path(run$paths$stimgate, "acs-manifest.rds"))
    expect_identical(manifest$settings$locThresholdMethod, expected)
    expect_no_error(env$.acsCytofValidateStimGateMethod(manifest, expected))
  }
})

test_that("the ACS runner keeps old output when a channel resolves otherwise", {
  env <- .load_acs_threshold_env()
  pathScratch <- tempfile("acs-threshold-")
  withr::defer(unlink(pathScratch, recursive = TRUE))
  pathStimgate <- env$.acsCytofPopulationPaths(
    "cd4", "raw", file.path(pathScratch, "gs"), pathScratch
  )$stimgate
  dir.create(pathStimgate, recursive = TRUE)
  writeLines("previous", file.path(pathStimgate, "previous.txt"))
  expect_error(
    .run_acs_runner_mocked(env, pathScratch, "region", "match"),
    "locThresholdMethod = 'region'"
  )
  expect_true(file.exists(file.path(pathStimgate, "previous.txt")))
  expect_false(file.exists(file.path(pathStimgate, "acs-manifest.rds")))
})

test_that("saved StimGate manifests without the requested method are rejected", {
  env <- .load_acs_threshold_env()
  validate <- env$.acsCytofValidateStimGateMethod
  expect_no_error(validate(.acs_stimgate_manifest(), "region"))
  legacy <- .acs_stimgate_manifest()
  legacy$settings$locThresholdMethod <- NULL
  legacy$channelSettings$IFNg$locThresholdMethod <- NULL
  expect_error(validate(legacy, "region"), "re-run all methods")
  expect_error(validate(.acs_stimgate_manifest("match"), "region"), "region")
  expect_error(
    validate(.acs_stimgate_manifest("region", "match"), "region"), "region"
  )
  noChannels <- .acs_stimgate_manifest()
  noChannels$channelSettings <- list()
  expect_error(validate(noChannels, "region"), "region")

  path <- tempfile("acs-stimgate-")
  withr::defer(unlink(path, recursive = TRUE))
  dir.create(file.path(path, "metaData"), recursive = TRUE)
  pathChnl <- file.path(path, "metaData", "chnlSettings.rds")
  saveRDS(legacy$channelSettings, pathChnl)
  saveRDS(legacy, file.path(path, "acs-manifest.rds"))
  expect_error(env$.acsCytofReadStimGateManifest(path, "region"), "region")
  current <- .acs_stimgate_manifest()
  saveRDS(current$channelSettings, pathChnl)
  saveRDS(current, file.path(path, "acs-manifest.rds"))
  expect_identical(env$.acsCytofReadStimGateManifest(path, "region"), current)
})

test_that("cached ACS comparisons must record the StimGate threshold method", {
  env <- .load_acs_threshold_env()
  table <- tibble::tibble(
    method = c("stimgate", "fbeta"), thresholdFailed = FALSE,
    locThresholdMethod = c("region", NA_character_)
  )
  fbeta <- list(context = .acs_stimgate_manifest()$context)
  manifest <- list(
    methods = list(
      cd4 = list(stimgate = .acs_stimgate_manifest(), fbeta = fbeta)
    ),
    comparisonSettings = list(
      methods = c("stimgate", "fbeta"), locThresholdMethod = "region"
    ),
    manualInputHash = "abc"
  )
  attr(table, "manifest") <- manifest
  expect_no_error(env$.acsCytofValidateComparisonManifest(table))
  expect_error(
    env$.acsCytofValidateComparisonManifest(
      table, locThresholdMethod = "match"
    ),
    "match"
  )

  legacy <- table
  attr(legacy, "manifest")$comparisonSettings$locThresholdMethod <- NULL
  expect_error(env$.acsCytofValidateComparisonManifest(legacy), "region")

  oldStimGate <- table
  attr(oldStimGate, "manifest")$methods$cd4$stimgate$settings <- list()
  expect_error(
    env$.acsCytofValidateComparisonManifest(oldStimGate), "region"
  )

  noColumn <- table
  noColumn$locThresholdMethod <- NULL
  attr(noColumn, "manifest") <- manifest
  expect_error(env$.acsCytofValidateComparisonManifest(noColumn), "rows")
})

test_that("the analysis 9 run manifest must record the threshold method", {
  env <- .load_acs_threshold_env()
  check <- env$.acsCytofCheckRunManifestMethod
  qmd <- "analysis/9-real-compare-acs-cytof.qmd"
  expect_no_error(check(
    list(stimgate_loc_threshold_method = "region"), "region", qmd
  ))
  expect_error(
    check(list(analysis_semantics_version = "acs-cytof-v4"), "region", qmd),
    "stimgate_loc_threshold_method"
  )
  expect_error(
    check(list(stimgate_loc_threshold_method = "match"), "region", qmd),
    "stimgate_loc_threshold_method"
  )
})

test_that("ACS StimGate output rows carry the threshold method", {
  env <- .load_acs_threshold_env()
  expect_identical(
    formals(env$comp_against_manual_cyt)$loc_threshold_method, "region"
  )
  expect_identical(
    formals(env$.acsCytofManualComparisonTable)$locThresholdMethod, "region"
  )
  expect_identical(
    formals(env$.acsCytofManualAutoTable)$locThresholdMethod, "region"
  )
  single <- tibble::tibble(
    method = "stimgate", ind = "2", cyt = "IFNg",
    freq_stim_auto = 1, freq_uns_auto = 0, freq_bs_auto = 1
  )
  thresholds <- tibble::tibble(
    ind = "2", cyt = "IFNg", threshold = 1, thresholdOrigin = "direct",
    thresholdFallbackUsed = FALSE, locThresholdMethod = "region"
  )
  out <- env$.acsCytofJoinProvenance(single, thresholds)
  expect_identical(out$locThresholdMethod, "region")
})

test_that("analyses 9 and 10 set and record the region threshold method", {
  qmd9 <- paste(readLines(qmd9_path, warn = FALSE), collapse = "\n")
  qmd10 <- paste(readLines(qmd10_path, warn = FALSE), collapse = "\n")
  setting <- 'stimgate_loc_threshold_method <- "region"'
  expect_true(grepl(setting, qmd9, fixed = TRUE))
  expect_true(grepl(setting, qmd10, fixed = TRUE))
  forwarded <- "locThresholdMethod = stimgate_loc_threshold_method"
  # Tester run, population runs and the cached-comparison validation.
  expect_gte(lengths(regmatches(
    qmd9, gregexpr(forwarded, qmd9, fixed = TRUE)
  )), 3L)
  expect_true(grepl(
    "loc_threshold_method = stimgate_loc_threshold_method", qmd9, fixed = TRUE
  ))
  expect_true(grepl(
    "stimgate_loc_threshold_method = stimgate_loc_threshold_method",
    qmd9, fixed = TRUE
  ))
  for (qmd in list(qmd9, qmd10)) {
    expect_true(grepl(".acsCytofCheckRunManifestMethod(", qmd, fixed = TRUE))
    expect_true(grepl(forwarded, qmd, fixed = TRUE))
    expect_true(grepl('"acs-cytof-v4"', qmd, fixed = TRUE))
    expect_false(grepl('"acs-cytof-v3"', qmd, fixed = TRUE))
  }
})
