root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

.bw_method_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "analysis-plot-style.R", "sim-misc.R",
    "sim-bandwidth.R", "sim-bandwidth-analysis-io.R",
    "sim-bandwidth-analysis-run.R", "sim-debug-loc.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  env
}

test_that("bandwidth QMDs set the threshold method explicitly", {
  for (qmd in c(
    "2a-sim-bw-freq_bs-global.qmd", "2b-sim-bias_uns-freq_bs.qmd",
    "6-sim-bw-freq_bs-adaptive.qmd"
  )) {
    txt <- paste(
      readLines(file.path(root_dir, "analysis", qmd)),
      collapse = "\n"
    )
    expect_match(txt, 'loc_threshold_method <- "cap"', fixed = TRUE)
    # scenario_settings is recorded in the manifest's required parameters.
    expect_match(
      txt, "locThresholdMethod = loc_threshold_method\n)",
      fixed = TRUE, info = qmd
    )
    expect_match(txt, "scenario_settings = scenario_settings", fixed = TRUE)
  }
  txt <- paste(
    readLines(file.path(root_dir, "analysis", "2c-sim-test.qmd")),
    collapse = "\n"
  )
  expect_match(txt, 'locThresholdMethod = "region"\n)', fixed = TRUE)
})

test_that("results recorded without the threshold method are rejected", {
  env <- .bw_method_env()
  current <- withr::local_tempdir()
  file.create(file.path(current, "COMPLETE"))
  saveRDS(1L, file.path(current, "result.rds"))
  ctx <- list(
    analysis_key = c("sim", "test"), current_dir = current,
    qmd_path = "analysis/test.qmd"
  )
  settings <- list(nSample = 2L, clusterGates = FALSE)
  saveRDS(
    list(
      analysis_key = ctx$analysis_key,
      params = list(analysis_semantics_version = "v1",
                    scenario_settings = settings)
    ),
    file.path(current, "manifest.rds")
  )
  required <- list(
    analysis_semantics_version = "v1",
    scenario_settings = c(settings, list(locThresholdMethod = "region"))
  )
  expect_error(
    env$.analysis_current_file(ctx, "result.rds", required),
    "required manifest parameters: scenario_settings", fixed = TRUE
  )
  required$scenario_settings$locThresholdMethod <- NULL
  expect_true(file.exists(
    env$.analysis_current_file(ctx, "result.rds", required)
  ))
})

test_that("bandwidth outputs and debug views record the threshold method", {
  testthat::skip_if_not_installed("simcyto")
  env <- .bw_method_env()
  expect_identical(
    formals(env$.simBandwidthBsFreq)$locThresholdMethod, "region"
  )
  row <- tibble::tibble(
    transformation = "skew", prob_response = 0.02, n_cell = 2000,
    mean_pos_setting = "high",
    mean_pos = env$.simMiscGetMeanPosTbl()$mean_pos[1],
    bw = 0.25, bias_uns_setting = "low", bias_uns = 0.05,
    sample_perturbation_sd = 0, condition_perturbation_sd = 0,
    cluster_perturbation_sd = 0, background_relative_to_response = 0.2,
    n_cell_uns_relative_to_stim = 1, sim_seed = 123L, sim_id = 1L
  )
  settings <- list(
    nSample = 2L, nMarker = 1, nCondition = 2, nCluster = 2, nIter = 1,
    bwMin = "none", bwMax = "none", probExact = TRUE, covEvMin = 1.5,
    covEvMax = 1.5, clusterGates = FALSE, locEnforceShapeThreshold = FALSE,
    calcCytPosGates = FALSE
  )
  quiet <- function(code) {
    out <- NULL
    utils::capture.output(out <- suppressMessages(code))
    out
  }
  run <- function(method) {
    quiet(env$.simDebugLoc(
      env$.simBandwidthRunRow(
        row, env$.simBandwidthFreqBsGlobalScenario,
        c(settings, list(locThresholdMethod = method))
      ),
      stopAfter = FALSE
    ))
  }
  region <- run("region")
  match <- run("match")

  res <- region$result
  expect_true(all(res$locThresholdMethod == "region"))
  cond <- res[res$method == "loc_condition", , drop = FALSE]
  expect_true(all(is.finite(cond$locRegionX)))
  direct <- cond$locGeneratedDirect %in% TRUE
  expect_equal(cond$threshold[direct], cond$locRegionX[direct])
  expect_true(all(match$result$locThresholdMethod == "match"))

  # Same simulated cells and filtering; only the gate choice differs.
  sumRegion <- env$.simDebugLocSummary(region)
  sumMatch <- env$.simDebugLocSummary(match)
  expect_identical(sumRegion$locThresholdMethod, "region")
  expect_identical(sumMatch$locThresholdMethod, "match")
  expect_equal(sumRegion$locRegionX, sumMatch$locRegionX)
  expect_identical(sumRegion$nGpStim, sumMatch$nGpStim)
  if (isTRUE(region$cp$locGeneratedDirect)) {
    expect_equal(sumRegion$threshold, sumRegion$locRegionX)
  }

  lines <- env$.simDebugLocLines(region)
  expect_true("gate (region method)" %in% lines$label)
  expect_true("xSum: region boundary" %in% lines$label)
  info <- env$.simDebugLocInfo(region, row)
  value <- function(tbl, nm) tbl$value[tbl$name == nm]
  expect_identical(value(info$gating, "locThresholdMethod"), "region")
  expect_identical(value(info$estimate, "threshold method"), "region")
  expect_identical(
    value(env$.simDebugLocInfo(match, row)$estimate, "threshold method"),
    "match"
  )
})
