root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

.compare_plot_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-plot-style.R", "sim-bandwidth-analysis-plot.R",
    "sim-compare-freq_bs.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  env
}

test_that(
  "analysis/8-sim-compare-freq_bs-batch.qmd does not source benchmarking cyt",
  {
    qmd_path <- file.path(
      root_dir,
      "analysis",
      "8-sim-compare-freq_bs-batch.qmd"
    )
    expect_true(file.exists(qmd_path))

    lines <- readLines(qmd_path, warn = FALSE)
    expect_false(
      any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
      info = "analysis 8 should not source functionsForBenchmarking-Cyt.R"
    )
  }
)

test_that(
  ".simCompareFreqBs forwards shift and sd multiplier to simcyto",
  {
  withr::local_preserve_seed()
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_misc, local = env)
    source(script_bw, local = env)
    source(script_comp, local = env)

    orig_simcyto_experiment <- simcyto::simCytExperiment
    captured_shift <- NULL
    captured_sd_mult <- NULL

    testthat::with_mocked_bindings(
      simCytExperiment = function(...,
                                  stimMeanShift = 0,
                                  stimSdMultiplier = 1) {
        captured_shift <<- stimMeanShift
        captured_sd_mult <<- stimSdMultiplier
        orig_simcyto_experiment(
          ...,
          stimMeanShift = stimMeanShift,
          stimSdMultiplier = stimSdMultiplier
        )
      },
      .package = "simcyto",
      {
        set.seed(42)
        res <- env$.simCompareFreqBs(
          nSample = 1L,
          nMarker = 1L,
          nCondition = 2L,
          nCluster = 2L,
          nIter = 1L,
          biasUns = 0,
          bw = 0.1,
          bwMtd = "hpi1",
          nCellStim = 200L,
          probResponse = 0.1,
          meanPos = 5,
          transformation = "gaussian",
          samplePerturbationSd = 0,
          conditionPerturbationSd = 0,
          clusterPerturbationSd = 0,
          backgroundRelativeToResponse = 0.1,
          ncellUnsRelativeToStim = 1,
          tailgateAutoTol = TRUE,
          stimMeanShift = 0.05,
          stimSdMultiplier = 1.05
        )

        expect_equal(captured_shift, 0.05)
        expect_equal(captured_sd_mult, 1.05)
        expect_s3_class(res, "data.frame")
        expect_true("stimMeanShift" %in% names(res))
        expect_true("stimSdMultiplier" %in% names(res))
        expect_equal(res$stimMeanShift[[1]], 0.05)
        expect_equal(res$stimSdMultiplier[[1]], 1.05)
      }
    )
  }
)

test_that(".simCompareSimCytExperiment applies selective mismatch exactly once", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_comp, local = env)

  make_out <- function() {
    labels <- c("gn", "gp", "gn", "gp")
    values <- c(1, 2, 3, 4)
    expr1 <- matrix(values, ncol = 1L)
    colnames(expr1) <- "F1"
    list(
      flowFrameList = list(
        flowCore::flowFrame(expr = expr1),
        flowCore::flowFrame(expr = expr1)
      ),
      labelsList = list(labels, labels),
      nCondition = 2L
    )
  }

  cluster_aware_seen <- NULL
  testthat::with_mocked_bindings(
    simCytExperiment = function(...,
                                stimMeanShift = 0,
                                stimSdMultiplier = 1,
                                stimMeanShiftClusters = NULL,
                                stimSdMultiplierClusters = NULL) {
      cluster_aware_seen <<- list(
        stimMeanShiftClusters = stimMeanShiftClusters,
        stimSdMultiplierClusters = stimSdMultiplierClusters,
        stimMeanShift = stimMeanShift,
        stimSdMultiplier = stimSdMultiplier
      )
      env$.simCompareApplyClusterMismatch(
        make_out(),
        stimMeanShift = stimMeanShift,
        stimSdMultiplier = stimSdMultiplier,
        stimMeanShiftClusters = stimMeanShiftClusters,
        stimSdMultiplierClusters = stimSdMultiplierClusters
      )
    },
    .package = "simcyto",
    {
      out_cluster_aware <- env$.simCompareSimCytExperiment(
        nSample = 1L,
        nMarker = 1L,
        nCondition = 2L,
        nCluster = 2L,
        nCellByCondition = c(4L, 4L),
        stimMeanShift = 2,
        stimSdMultiplier = 1.5,
        stimMeanShiftClusters = "gn",
        stimSdMultiplierClusters = "gp"
      )

      expect_equal(cluster_aware_seen$stimMeanShiftClusters, "gn")
      expect_equal(cluster_aware_seen$stimSdMultiplierClusters, "gp")
      expect_equal(cluster_aware_seen$stimMeanShift, 2)
      expect_equal(cluster_aware_seen$stimSdMultiplier, 1.5)
      expect_equal(
        as.vector(flowCore::exprs(out_cluster_aware[["flowFrameList"]][[2L]])),
        c(3, 1.5, 5, 4.5),
        tolerance = 1e-8
      )
    }
  )

  legacy_seen <- NULL
  testthat::with_mocked_bindings(
    simCytExperiment = function(...,
                                stimMeanShift = 0,
                                stimSdMultiplier = 1) {
      legacy_seen <<- list(
        stimMeanShift = stimMeanShift,
        stimSdMultiplier = stimSdMultiplier,
        dots = list(...)
      )
      make_out()
    },
    .package = "simcyto",
    {
      out_legacy <- env$.simCompareSimCytExperiment(
        nSample = 1L,
        nMarker = 1L,
        nCondition = 2L,
        nCluster = 2L,
        nCellByCondition = c(4L, 4L),
        stimMeanShift = 2,
        stimSdMultiplier = 1.5,
        stimMeanShiftClusters = "gn",
        stimSdMultiplierClusters = "gp"
      )

      expect_equal(legacy_seen$stimMeanShift, 0)
      expect_equal(legacy_seen$stimSdMultiplier, 1)
      expect_false("stimMeanShiftClusters" %in% names(legacy_seen$dots))
      expect_false("stimSdMultiplierClusters" %in% names(legacy_seen$dots))
      expect_equal(
        as.vector(flowCore::exprs(out_legacy[["flowFrameList"]][[2L]])),
        c(3, 1.5, 5, 4.5),
        tolerance = 1e-8
      )
    }
  )
})

test_that(".simCompareFreqBs with zero mismatch reproduces clean baseline", {
  withr::local_preserve_seed()
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  set.seed(123)
  res_clean <- env$.simCompareFreqBs(
    nSample = 2L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0,
    bw = 0.1,
    bwMtd = "hpi1",
    nCellStim = 300L,
    probResponse = 0.05,
    meanPos = 5,
    transformation = "gaussian",
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.1,
    ncellUnsRelativeToStim = 1,
    tailgateAutoTol = TRUE
  )

  set.seed(123)
  res_zero_mismatch <- env$.simCompareFreqBs(
    nSample = 2L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0,
    bw = 0.1,
    bwMtd = "hpi1",
    nCellStim = 300L,
    probResponse = 0.05,
    meanPos = 5,
    transformation = "gaussian",
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.1,
    ncellUnsRelativeToStim = 1,
    tailgateAutoTol = TRUE,
    stimMeanShift = 0,
    stimSdMultiplier = 1
  )

  # Disregarding the added stimMeanShift and stimSdMultiplier columns
  common_cols <- intersect(names(res_clean), names(res_zero_mismatch))
  expect_equal(res_clean[common_cols], res_zero_mismatch[common_cols])
})

test_that(
  ".simCompareSummariseFreqBs correctly handles mismatch scenarios",
  {
  withr::local_preserve_seed()
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_misc, local = env)
    source(script_bw, local = env)
    source(script_comp, local = env)

    grid <- data.frame(
      transformation = c("gaussian", "gaussian"),
      mean_pos = c(5, 5),
      prob_response = c(0.1, 0.1),
      n_cell = c(200, 200),
      bias_uns = c(0.15, 0.15),
      bw = c(0.1, 0.1),
      sample_perturbation_sd = c(0, 0),
      condition_perturbation_sd = c(0, 0),
      cluster_perturbation_sd = c(0, 0),
      background_relative_to_response = c(0.1, 0.1),
      n_cell_uns_relative_to_stim = c(1, 1),
      stim_mean_shift = c(0, 0.05),
      stim_sd_multiplier = c(1, 1),
      mismatch_type = c("mean_shift", "mean_shift"),
      mismatch_val = c(0, 0.05),
      stringsAsFactors = FALSE
    )

    set.seed(42)
    raw_res <- env$.simCompareFreqBsGrid(
      sim_grid = grid,
      nSample = 1,
      nIter = 1,
      nMarker = 1,
      nCondition = 2,
      nCluster = 2,
      probExact = TRUE,
      tailgateAutoTol = TRUE
    )

    expect_s3_class(raw_res, "data.frame")
    expect_true(nrow(raw_res) > 0)

    summ <- env$.simCompareSummariseFreqBs(raw_res)
    expect_s3_class(summ, "data.frame")
    req_cols <- c(
      "med_abs_rel_error",
      "q90_abs_rel_error",
      "q95_abs_rel_error",
      "max_abs_rel_error"
    )
    expect_true(all(req_cols %in% names(summ)))
    expect_true("mismatch_val" %in% names(summ))
    expect_true(nrow(summ) > 0)
  }
)

test_that(
  ".simCompareFreqBsGrid parallel and serial runs produce equivalent results",
  {
  withr::local_preserve_seed()
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_misc, local = env)
    source(script_bw, local = env)
    source(script_comp, local = env)

    grid <- data.frame(
      sim_id = c(1L, 2L),
      transformation = c("gaussian", "gaussian"),
      mean_pos = c(5, 5),
      prob_response = c(0.1, 0.1),
      n_cell = c(100, 100),
      bias_uns = c(0, 0),
      bw = c(0.1, 0.1),
      sample_perturbation_sd = c(0, 0),
      condition_perturbation_sd = c(0, 0),
      cluster_perturbation_sd = c(0, 0),
      background_relative_to_response = c(0.1, 0.1),
      n_cell_uns_relative_to_stim = c(1, 1),
      stim_mean_shift = c(0, 0.025),
      stim_sd_multiplier = c(1, 1),
      mismatch_type = c("mean_shift", "mean_shift"),
      mismatch_val = c(0, 0.025),
      stringsAsFactors = FALSE
    )

    set.seed(999)
    res_serial <- env$.simCompareFreqBsGrid(
      sim_grid = grid,
      nSample = 1,
      nIter = 1,
      nMarker = 1,
      nCondition = 2,
      nCluster = 2,
      probExact = TRUE,
      tailgateAutoTol = TRUE,
      parallel = FALSE,
      progress = FALSE
    )

    set.seed(999)
    res_parallel <- env$.simCompareFreqBsGrid(
      sim_grid = grid,
      nSample = 1,
      nIter = 1,
      nMarker = 1,
      nCondition = 2,
      nCluster = 2,
      probExact = TRUE,
      tailgateAutoTol = TRUE,
      parallel = TRUE,
      workers = 2L,
      progress = FALSE
    )

    expect_equal(res_serial$sim_id, res_parallel$sim_id)
    expect_equal(res_serial$propRespEst, res_parallel$propRespEst)
    expect_equal(res_serial$threshold, res_parallel$threshold)
    expect_equal(res_serial$propRespTruth, res_parallel$propRespTruth)
  }
)

test_that(".simCompareFreqBsGrid supports per-scenario caching and resume", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  tmp_cache <- file.path(
    tempdir(),
    paste0("test_sim_cache_", Sys.getpid(), "_", sample.int(1e6, 1))
  )
  tmp_log <- file.path(tmp_cache, "progress.txt")
  dir.create(tmp_cache, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_cache, recursive = TRUE, force = TRUE), add = TRUE)

  grid <- data.frame(
    sim_id = c(1L, 2L),
    transformation = c("gaussian", "gaussian"),
    mean_pos = c(5, 5),
    prob_response = c(0.1, 0.1),
    n_cell = c(100, 100),
    bias_uns = c(0, 0),
    bw = c(0.1, 0.1),
    sample_perturbation_sd = c(0, 0),
    condition_perturbation_sd = c(0, 0),
    cluster_perturbation_sd = c(0, 0),
    background_relative_to_response = c(0.1, 0.1),
    n_cell_uns_relative_to_stim = c(1, 1),
    stim_mean_shift = c(0, 0.05),
    stim_sd_multiplier = c(1, 1),
    mismatch_type = c("mean_shift", "mean_shift"),
    mismatch_val = c(0, 0.05),
    stringsAsFactors = FALSE
  )

  res1 <- env$.simCompareFreqBsGrid(
    sim_grid = grid,
    nSample = 1,
    nIter = 1,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    probExact = TRUE,
    tailgateAutoTol = TRUE,
    dirCache = tmp_cache,
    pathProgress = tmp_log,
    resume = TRUE,
    parallel = FALSE,
    progress = FALSE
  )

  cached_files <- env$.simCompareFindScenarioOutputs(tmp_cache)
  expect_equal(length(cached_files), 2L)

  mtimes_before <- file.info(cached_files)$mtime

  # Run again with resume = TRUE: cached files should not be rewritten
  res2 <- env$.simCompareFreqBsGrid(
    sim_grid = grid,
    nSample = 1,
    nIter = 1,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    probExact = TRUE,
    tailgateAutoTol = TRUE,
    dirCache = tmp_cache,
    pathProgress = tmp_log,
    resume = TRUE,
    parallel = FALSE,
    progress = FALSE
  )

  mtimes_after <- file.info(cached_files)$mtime
  expect_equal(mtimes_before, mtimes_after)
  expect_equal(res1$propRespEst, res2$propRespEst)

  # Check progress log records skipped scenarios
  log_lines <- readLines(tmp_log, warn = FALSE)
  expect_true(any(grepl("Skipped", log_lines)))

  # Invalidate scenario 1 by changing mismatch setting
  grid_mod <- grid
  grid_mod$stim_mean_shift[1] <- 0.1
  grid_mod$mismatch_val[1] <- 0.1

  res3 <- env$.simCompareFreqBsGrid(
    sim_grid = grid_mod,
    nSample = 1,
    nIter = 1,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    probExact = TRUE,
    tailgateAutoTol = TRUE,
    dirCache = tmp_cache,
    pathProgress = tmp_log,
    resume = TRUE,
    parallel = FALSE,
    progress = FALSE
  )

  # Scenario 1 was recomputed with new mismatch setting
  expect_equal(res3$stim_mean_shift[res3$sim_id == 1L][[1]], 0.1)
  # Scenario 2 was skipped from cache
  expect_equal(res3$stim_mean_shift[res3$sim_id == 2L][[1]], 0.05)
})

test_that(".simCompareCollateScenarioOutputs collates scenario files", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  tmp_cache <- file.path(
    tempdir(),
    paste0("test_collate_", Sys.getpid(), "_", sample.int(1e6, 1))
  )
  dir.create(tmp_cache, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_cache, recursive = TRUE, force = TRUE), add = TRUE)

  f1 <- file.path(tmp_cache, "compare_raw-sim_id_000001.rds")
  f2 <- file.path(tmp_cache, "compare_raw-sim_id_000002.rds")

  d1 <- data.frame(
    sim_id = 1L,
    propRespEst = 0.05,
    approach = "stimgate",
    method = "stimgate_loc_sample",
    stringsAsFactors = FALSE
  )
  d2 <- data.frame(
    sim_id = 2L,
    propRespEst = 0.10,
    approach = "stimgate",
    method = "stimgate_loc_sample",
    stringsAsFactors = FALSE
  )

  saveRDS(d1, f1)
  saveRDS(d2, f2)

  collated <- env$.simCompareCollateScenarioOutputs(dirCache = tmp_cache)
  expect_s3_class(collated, "data.frame")
  expect_equal(nrow(collated), 2L)
  expect_equal(collated$sim_id, c(1L, 2L))
  expect_equal(collated$propRespEst, c(0.05, 0.10))
})

test_that("Analysis 7 and Analysis 8 namespaces are isolated", {
  qmd7_path <- file.path(root_dir, "analysis", "7-sim-compare-freq_bs.qmd")
  qmd8_path <- file.path(
    root_dir,
    "analysis",
    "8-sim-compare-freq_bs-batch.qmd"
  )

  expect_true(file.exists(qmd7_path))
  expect_true(file.exists(qmd8_path))

  lines7 <- readLines(qmd7_path, warn = FALSE)
  lines8 <- readLines(qmd8_path, warn = FALSE)

  # Analysis 7 should not reference freq_bs_batch
  expect_false(any(grepl("freq_bs_batch", lines7)))

  # Analysis 8 should use dedicated freq_bs_batch namespace
  expect_true(any(grepl('"freq_bs_batch"', lines8)))

  content8 <- paste(lines8, collapse = "\n")
  # Analysis 8 run context should be in freq_bs_batch namespace
  expect_true(grepl('c\\("sim",\\s*"compare",\\s*"freq_bs_batch"\\)', content8))
  # The shared runtime derives the run-state namespace from analysis_key.
  expect_true(
    grepl(
      'analysis_key\\s*<-\\s*c\\([^)]*"freq_bs_batch"',
      content8
    ) &&
      grepl("analysis_key = analysis_key", content8, fixed = TRUE)
  )
})

test_that("Inner nIter execution remains serial without nested parallelism", {
  # Verify that .simCompareFreqBs executes iterations sequentially via
  # purrr::map_df and does not call future::plan() or furrr inside
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)
  fn_body_txt <- paste(deparse(body(env$.simCompareFreqBs)), collapse = "\n")

  expect_true(grepl("purrr::map_df", fn_body_txt))
  expect_false(grepl("future_map", fn_body_txt))
  expect_false(grepl("future::plan", fn_body_txt))
})

test_that(".simCompareRunScenario handles errors and writes log", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  tmp_cache <- file.path(
    tempdir(),
    paste0("test_err_", Sys.getpid(), "_", sample.int(1e6, 1))
  )
  tmp_log <- file.path(tmp_cache, "progress.txt")
  dir.create(tmp_cache, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_cache, recursive = TRUE, force = TRUE), add = TRUE)

  # Invalid scenario that causes an error in .simCompareFreqBs
  bad_row <- data.frame(
    sim_id = 99L,
    transformation = "nonexistent_trans_xyz",
    mean_pos = 5,
    prob_response = 0.1,
    n_cell = 100,
    bias_uns = 0,
    mismatch_type = "none",
    mismatch_val = 0,
    stringsAsFactors = FALSE
  )

  err_res <- env$.simCompareRunScenario(
    row = bad_row,
    nSample = 1,
    nIter = 1,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    probExact = TRUE,
    dirCache = tmp_cache,
    pathProgress = tmp_log
  )

  expect_s3_class(err_res, "data.frame")
  expect_true("error" %in% names(err_res))
  expect_true(!is.na(err_res$error[[1]]) && nzchar(err_res$error[[1]]))

  log_lines <- readLines(tmp_log, warn = FALSE)
  expect_true(any(grepl("Error \\[sim_id = 99", log_lines)))
})

test_that(
  "analysis/8-sim-compare-freq_bs-batch.qmd includes mean_shift_negative",
  {
    qmd_path <- file.path(
      root_dir,
      "analysis",
      "8-sim-compare-freq_bs-batch.qmd"
    )
    lines <- readLines(qmd_path, warn = FALSE)
    expect_true(any(grepl('"mean_shift_negative"', lines)))
    expect_true(any(grepl('"mean_shift_all"', lines)))
    expect_true(any(grepl('"sd_inflation"', lines)))
    expect_true(any(grepl('stim_mean_shift_clusters = "gn"', lines)))
  }
)

test_that(
  ".simCompareFreqBs forwards stimMeanShiftClusters to simcyto",
  {
  withr::local_preserve_seed()
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_misc, local = env)
    source(script_bw, local = env)
    source(script_comp, local = env)

    orig_simcyto_experiment <- simcyto::simCytExperiment
    captured_clusters <- NULL

    testthat::with_mocked_bindings(
      simCytExperiment = function(...,
                                  stimMeanShift = 0,
                                  stimSdMultiplier = 1,
                                  stimMeanShiftClusters = NULL) {
        captured_clusters <<- stimMeanShiftClusters
        call_args <- list(...)
        call_args$stimMeanShift <- stimMeanShift
        call_args$stimSdMultiplier <- stimSdMultiplier
        if ("stimMeanShiftClusters" %in% names(formals(orig_simcyto_experiment))) {
          call_args$stimMeanShiftClusters <- stimMeanShiftClusters
        }
        do.call(orig_simcyto_experiment, call_args)
      },
      .package = "simcyto",
      {
        set.seed(42)
        res <- env$.simCompareFreqBs(
          nSample = 1L,
          nMarker = 1L,
          nCondition = 2L,
          nCluster = 2L,
          nIter = 1L,
          biasUns = 0,
          bw = 0.1,
          bwMtd = "hpi1",
          nCellStim = 200L,
          probResponse = 0.1,
          meanPos = 5,
          transformation = "gaussian",
          samplePerturbationSd = 0,
          conditionPerturbationSd = 0,
          clusterPerturbationSd = 0,
          backgroundRelativeToResponse = 0.1,
          ncellUnsRelativeToStim = 1,
          tailgateAutoTol = TRUE,
          stimMeanShift = 0.05,
          stimMeanShiftClusters = "gn"
        )

        expect_equal(captured_clusters, "gn")
        expect_s3_class(res, "data.frame")
        expect_true("stimMeanShiftClusters" %in% names(res))
        expect_equal(res$stimMeanShiftClusters[[1]], "gn")
      }
    )
  }
)

test_that(
  "cache validation distinguishes all-component and negative-only shift",
  {
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_misc, local = env)
    source(script_bw, local = env)
    source(script_comp, local = env)

    cached_all <- data.frame(
      sim_id = 1L,
      iter = 1L,
      sample = "sample1",
      transformation = "gaussian",
      mismatch_type = "mean_shift_all",
      mismatch_val = 0.05,
      stim_mean_shift = 0.05,
      stim_mean_shift_clusters = NA_character_,
      stringsAsFactors = FALSE
    )

    row_neg <- data.frame(
      sim_id = 1L,
      transformation = "gaussian",
      mismatch_type = "mean_shift_negative",
      mismatch_val = 0.05,
      stim_mean_shift = 0.05,
      stim_mean_shift_clusters = "gn",
      stringsAsFactors = FALSE
    )

    # Cached all-component output should NOT validate for negative-only row
    expect_false(
      env$.simCompareValidateScenarioCache(
        cached = cached_all,
        row = row_neg,
        nSample = 1,
        nIter = 1
      )
    )

    cached_neg <- data.frame(
      sim_id = 1L,
      iter = 1L,
      sample = "sample1",
      transformation = "gaussian",
      mismatch_type = "mean_shift_negative",
      mismatch_val = 0.05,
      stim_mean_shift = 0.05,
      stim_mean_shift_clusters = "gn",
      stringsAsFactors = FALSE
    )

    # Cached negative-only output SHOULD validate for negative-only row
    expect_true(
      env$.simCompareValidateScenarioCache(
        cached = cached_neg,
        row = row_neg,
        nSample = 1,
        nIter = 1
      )
    )
  }
)

test_that(
  "negative-only zero-shift agrees with clean baseline semantics",
  {
  withr::local_preserve_seed()
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_misc, local = env)
    source(script_bw, local = env)
    source(script_comp, local = env)

    set.seed(42)
    res_clean <- env$.simCompareFreqBs(
      nSample = 1L,
      nMarker = 1L,
      nCondition = 2L,
      nCluster = 2L,
      nIter = 1L,
      biasUns = 0,
      bw = 0.1,
      bwMtd = "hpi1",
      nCellStim = 200L,
      probResponse = 0.05,
      meanPos = 5,
      transformation = "gaussian",
      samplePerturbationSd = 0,
      conditionPerturbationSd = 0,
      clusterPerturbationSd = 0,
      backgroundRelativeToResponse = 0.1,
      ncellUnsRelativeToStim = 1,
      tailgateAutoTol = TRUE,
      stimMeanShift = 0,
      stimSdMultiplier = 1
    )

    set.seed(42)
    res_zero_neg <- env$.simCompareFreqBs(
      nSample = 1L,
      nMarker = 1L,
      nCondition = 2L,
      nCluster = 2L,
      nIter = 1L,
      biasUns = 0,
      bw = 0.1,
      bwMtd = "hpi1",
      nCellStim = 200L,
      probResponse = 0.05,
      meanPos = 5,
      transformation = "gaussian",
      samplePerturbationSd = 0,
      conditionPerturbationSd = 0,
      clusterPerturbationSd = 0,
      backgroundRelativeToResponse = 0.1,
      ncellUnsRelativeToStim = 1,
      tailgateAutoTol = TRUE,
      stimMeanShift = 0,
      stimSdMultiplier = 1,
      stimMeanShiftClusters = "gn"
    )

    common_cols <- setdiff(
      intersect(names(res_clean), names(res_zero_neg)),
      "stimMeanShiftClusters"
    )
    expect_equal(res_clean[common_cols], res_zero_neg[common_cols])
  }
)

test_that(
  "analysis/8-sim-compare-freq_bs-batch.qmd targets gn for SD inflation",
  {
    qmd_path <- file.path(
      root_dir,
      "analysis",
      "8-sim-compare-freq_bs-batch.qmd"
    )
    lines <- readLines(qmd_path, warn = FALSE)
    expect_true(any(grepl('stim_sd_multiplier_clusters = "gn"', lines)))
  }
)

test_that(
  ".simCompareFreqBs forwards stimSdMultiplierClusters to simcyto",
  {
  withr::local_preserve_seed()
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_misc, local = env)
    source(script_bw, local = env)
    source(script_comp, local = env)

    orig_simcyto_experiment <- simcyto::simCytExperiment
    captured_sd_clusters <- NULL

    testthat::with_mocked_bindings(
      simCytExperiment = function(...,
                                  stimMeanShift = 0,
                                  stimSdMultiplier = 1,
                                  stimSdMultiplierClusters = NULL) {
        captured_sd_clusters <<- stimSdMultiplierClusters
        call_args <- list(...)
        call_args$stimMeanShift <- stimMeanShift
        call_args$stimSdMultiplier <- stimSdMultiplier
        if ("stimSdMultiplierClusters" %in% names(formals(orig_simcyto_experiment))) {
          call_args$stimSdMultiplierClusters <- stimSdMultiplierClusters
        }
        do.call(orig_simcyto_experiment, call_args)
      },
      .package = "simcyto",
      {
        set.seed(42)
        res <- env$.simCompareFreqBs(
          nSample = 1L,
          nMarker = 1L,
          nCondition = 2L,
          nCluster = 2L,
          nIter = 1L,
          biasUns = 0,
          bw = 0.1,
          bwMtd = "hpi1",
          nCellStim = 200L,
          probResponse = 0.1,
          meanPos = 5,
          transformation = "gaussian",
          samplePerturbationSd = 0,
          conditionPerturbationSd = 0,
          clusterPerturbationSd = 0,
          backgroundRelativeToResponse = 0.1,
          ncellUnsRelativeToStim = 1,
          tailgateAutoTol = TRUE,
          stimSdMultiplier = 1.2,
          stimSdMultiplierClusters = "gn"
        )

        expect_equal(captured_sd_clusters, "gn")
        expect_s3_class(res, "data.frame")
        expect_true("stimSdMultiplierClusters" %in% names(res))
        expect_equal(res$stimSdMultiplierClusters[[1]], "gn")
      }
    )
  }
)

test_that(
  "stimSdMultiplierClusters = 'gn' leaves positive and unstim unchanged",
  {
  withr::local_preserve_seed()
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_comp, local = env)

    set.seed(42)
    trans <- simcyto::simCytTransformIdentity()
    meanExprMat <- matrix(c(0, 5), byrow = TRUE, ncol = 1)
    clusterLabelVec <- c("gn", "gp")
    probVecUns <- c(0.9, 0.1)
    probResponseVecByStimCondition <- list(c(-0.05, 0.05))

    sim_clean <- simcyto::simCytExperiment(
      nSample = 1,
      nMarker = 1,
      nCondition = 2,
      nCluster = 2,
      nCellByCondition = c(1000, 1000),
      transformationFunc = trans,
      mixtureType = "gaussianOnly",
      meanExprMat = meanExprMat,
      clusterLabelVec = clusterLabelVec,
      probVecUns = probVecUns,
      probExact = TRUE,
      probResponseVecByStimCondition = probResponseVecByStimCondition,
      covEvMin = 1,
      covEvMax = 1,
      stimMeanShift = 0,
      stimSdMultiplier = 1
    )

    set.seed(42)
    sim_neg_sd <- env$.simCompareApplyClusterMismatch(
      outListExperiment = sim_clean,
      stimMeanShift = 0,
      stimSdMultiplier = 1.25,
      stimMeanShiftClusters = NULL,
      stimSdMultiplierClusters = "gn"
    )

    unstim_clean <- flowCore::exprs(sim_clean$flowFrameList[[1]])[, "F1"]
    unstim_neg <- flowCore::exprs(sim_neg_sd$flowFrameList[[1]])[, "F1"]
    expect_equal(unstim_clean, unstim_neg)

    stim_clean_expr <- flowCore::exprs(sim_clean$flowFrameList[[2]])[, "F1"]
    stim_neg_expr <- flowCore::exprs(sim_neg_sd$flowFrameList[[2]])[, "F1"]
    labels <- sim_clean$labelsList[[2]]

    expect_equal(
      stim_clean_expr[labels == "gp"],
      stim_neg_expr[labels == "gp"]
    )

    sd_clean_gn <- sd(stim_clean_expr[labels == "gn"])
    sd_neg_gn <- sd(stim_neg_expr[labels == "gn"])
    expect_equal(sd_neg_gn / sd_clean_gn, 1.25, tolerance = 1e-10)
  }
)

test_that(
  "cache validation distinguishes all-component and negative-only SD inflation",
  {
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_misc, local = env)
    source(script_bw, local = env)
    source(script_comp, local = env)

    cached_all_sd <- data.frame(
      sim_id = 2L,
      iter = 1L,
      sample = "sample1",
      transformation = "gaussian",
      mismatch_type = "sd_inflation",
      mismatch_val = 0.10,
      stim_sd_multiplier = 1.10,
      stim_sd_multiplier_clusters = NA_character_,
      stringsAsFactors = FALSE
    )

    row_neg_sd <- data.frame(
      sim_id = 2L,
      transformation = "gaussian",
      mismatch_type = "sd_inflation",
      mismatch_val = 0.10,
      stim_sd_multiplier = 1.10,
      stim_sd_multiplier_clusters = "gn",
      stringsAsFactors = FALSE
    )

    # Cached all-component SD output should NOT validate for negative-only row
    expect_false(
      env$.simCompareValidateScenarioCache(
        cached = cached_all_sd,
        row = row_neg_sd,
        nSample = 1,
        nIter = 1
      )
    )

    # Cached object missing the stim_sd_multiplier_clusters column must NOT validate
    cached_legacy <- data.frame(
      sim_id = 2L,
      iter = 1L,
      sample = "sample1",
      transformation = "gaussian",
      mismatch_type = "sd_inflation",
      mismatch_val = 0.10,
      stim_sd_multiplier = 1.10,
      stringsAsFactors = FALSE
    )
    expect_false(
      env$.simCompareValidateScenarioCache(
        cached = cached_legacy,
        row = row_neg_sd,
        nSample = 1,
        nIter = 1
      )
    )

    cached_neg_sd <- data.frame(
      sim_id = 2L,
      iter = 1L,
      sample = "sample1",
      transformation = "gaussian",
      mismatch_type = "sd_inflation",
      mismatch_val = 0.10,
      stim_sd_multiplier = 1.10,
      stim_sd_multiplier_clusters = "gn",
      stringsAsFactors = FALSE
    )

    # Cached negative-only SD output SHOULD validate for negative-only row
    expect_true(
      env$.simCompareValidateScenarioCache(
        cached = cached_neg_sd,
        row = row_neg_sd,
        nSample = 1,
        nIter = 1
      )
    )
  }
)

test_that(
  "negative-only zero-increase SD inflation agrees with clean baseline",
  {
  withr::local_preserve_seed()
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_misc, local = env)
    source(script_bw, local = env)
    source(script_comp, local = env)

    set.seed(42)
    res_clean <- env$.simCompareFreqBs(
      nSample = 1L,
      nMarker = 1L,
      nCondition = 2L,
      nCluster = 2L,
      nIter = 1L,
      biasUns = 0,
      bw = 0.1,
      bwMtd = "hpi1",
      nCellStim = 200L,
      probResponse = 0.05,
      meanPos = 5,
      transformation = "gaussian",
      samplePerturbationSd = 0,
      conditionPerturbationSd = 0,
      clusterPerturbationSd = 0,
      backgroundRelativeToResponse = 0.1,
      ncellUnsRelativeToStim = 1,
      tailgateAutoTol = TRUE,
      stimMeanShift = 0,
      stimSdMultiplier = 1
    )

    set.seed(42)
    res_zero_sd <- env$.simCompareFreqBs(
      nSample = 1L,
      nMarker = 1L,
      nCondition = 2L,
      nCluster = 2L,
      nIter = 1L,
      biasUns = 0,
      bw = 0.1,
      bwMtd = "hpi1",
      nCellStim = 200L,
      probResponse = 0.05,
      meanPos = 5,
      transformation = "gaussian",
      samplePerturbationSd = 0,
      conditionPerturbationSd = 0,
      clusterPerturbationSd = 0,
      backgroundRelativeToResponse = 0.1,
      ncellUnsRelativeToStim = 1,
      tailgateAutoTol = TRUE,
      stimMeanShift = 0,
      stimSdMultiplier = 1,
      stimSdMultiplierClusters = "gn"
    )

    common_cols <- setdiff(
      intersect(names(res_clean), names(res_zero_sd)),
      c("stimMeanShiftClusters", "stimSdMultiplierClusters")
    )
    expect_equal(res_clean[common_cols], res_zero_sd[common_cols])
  }
)


test_that(".simComparePrimaryOutputComplete requires exact primary coverage", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  primary <- tidyr::expand_grid(
    iter = 1:2,
    sample = as.character(1:2),
    method = c("stimgate", "fbeta", "tailgate")
  ) |>
    dplyr::mutate(
      propRespTruth = 0.1,
      propRespEst = 0.1,
      nCellStim = 100,
      nPosStim = 12L,
      nTruePos = 10L,
      nFalsePos = 2L,
      nFalseNeg = 1L,
      nTrueNeg = 87L,
      unsExprSum = 1.5,
      error = NA_character_
    )

  expect_true(
    env$.simComparePrimaryOutputComplete(
      primary,
      nSample = 2,
      nIter = 2
    )
  )
  expect_false(
    env$.simComparePrimaryOutputComplete(
      primary[-1, , drop = FALSE],
      nSample = 2,
      nIter = 2
    )
  )

  duplicated <- dplyr::bind_rows(primary, primary[1, , drop = FALSE])
  expect_false(
    env$.simComparePrimaryOutputComplete(
      duplicated,
      nSample = 2,
      nIter = 2
    )
  )

  failed <- primary
  failed$error[[1]] <- "competitor failed"
  expect_false(
    env$.simComparePrimaryOutputComplete(
      failed,
      nSample = 2,
      nIter = 2
    )
  )

  nonfinite <- primary
  nonfinite$propRespEst[[1]] <- NA_real_
  expect_false(
    env$.simComparePrimaryOutputComplete(
      nonfinite,
      nSample = 2,
      nIter = 2
    )
  )

  # Outputs saved before the label-based counts existed are incomplete.
  for (col in c("nTruePos", "nFalsePos", "nFalseNeg", "nTrueNeg", "unsExprSum")) {
    legacy <- primary[, setdiff(names(primary), col)]
    expect_false(
      env$.simComparePrimaryOutputComplete(legacy, nSample = 2, nIter = 2),
      info = col
    )
  }
  missing_count <- primary
  missing_count$nTruePos[[1]] <- NA_integer_
  expect_false(
    env$.simComparePrimaryOutputComplete(missing_count, nSample = 2, nIter = 2)
  )
  # Counts must reproduce the method's own gated count and the tube size.
  wrong_gated <- primary
  wrong_gated$nPosStim[[1]] <- 13L
  expect_false(
    env$.simComparePrimaryOutputComplete(wrong_gated, nSample = 2, nIter = 2)
  )
  wrong_total <- primary
  wrong_total$nTrueNeg[[1]] <- 86L
  expect_false(
    env$.simComparePrimaryOutputComplete(wrong_total, nSample = 2, nIter = 2)
  )
})

test_that(".simCompareGridOutputStatus catches missing, failed, and empty chunks", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  make_scenario <- function(sim_id) {
    tidyr::expand_grid(
      sim_id = sim_id,
      iter = 1:2,
      sample = as.character(1:2),
      method = c("stimgate", "fbeta", "tailgate")
    ) |>
      dplyr::mutate(
        propRespTruth = 0.1,
        propRespEst = 0.1,
        nCellStim = 100,
        nPosStim = 12L,
        nTruePos = 10L,
        nFalsePos = 2L,
        nFalseNeg = 1L,
        nTrueNeg = 87L,
        unsExprSum = 1.5,
        error = NA_character_
      )
  }

  grid <- tibble::tibble(sim_id = 1:2)
  complete <- dplyr::bind_rows(make_scenario(1L), make_scenario(2L))
  status <- env$.simCompareGridOutputStatus(
    complete,
    sim_grid = grid,
    nSample = 2,
    nIter = 2
  )
  expect_true(status$collate_ok)
  expect_true(status$validation_ok)
  expect_equal(status$completed_ids, 1:2)
  expect_length(status$failed_ids, 0L)

  missing <- env$.simCompareGridOutputStatus(
    make_scenario(1L),
    sim_grid = grid,
    nSample = 2,
    nIter = 2
  )
  expect_false(missing$collate_ok)
  expect_false(missing$validation_ok)
  expect_equal(missing$missing_ids, 2L)

  failed_data <- complete
  failed_idx <- which(failed_data$sim_id == 2L)[[1]]
  failed_data$error[[failed_idx]] <- "failed"
  failed <- env$.simCompareGridOutputStatus(
    failed_data,
    sim_grid = grid,
    nSample = 2,
    nIter = 2
  )
  expect_true(failed$collate_ok)
  expect_false(failed$validation_ok)
  expect_equal(failed$failed_ids, 2L)

  empty <- env$.simCompareGridOutputStatus(
    tibble::tibble(),
    sim_grid = tibble::tibble(sim_id = integer()),
    nSample = 2,
    nIter = 2
  )
  expect_true(empty$collate_ok)
  expect_true(empty$validation_ok)
  expect_length(empty$completed_ids, 0L)
})

test_that("alternative comparator exceptions remain explicit run errors", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  assign(
    ".simCompareFbetaEnvironment",
    function(...) new.env(parent = emptyenv()),
    envir = env
  )
  assign(
    ".simCompareFbetaThreshold",
    function(...) stop("fbeta boom"),
    envir = env
  )
  assign(
    ".simCompareTailgateThreshold",
    function(...) {
      list(
        threshold = 0,
        thresholdMetric = NA_real_,
        thresholdOrigin = "calculated"
      )
    },
    envir = env
  )

  x_uns <- matrix(c(-2, -1, 0, 1), ncol = 1)
  x_stim <- matrix(c(-1, 0, 1, 2), ncol = 1)
  colnames(x_uns) <- "F1"
  colnames(x_stim) <- "F1"

  flow_frames <- list(
    flowCore::flowFrame(expr = x_uns),
    flowCore::flowFrame(expr = x_stim)
  )
  labels <- list(
    c("gn", "gn", "gp", "gp"),
    c("gn", "gn", "gp", "gp")
  )

  res <- env$.simCompareAlternativeRows(
    flowFrameList = flow_frames,
    labelsList = labels,
    nSample = 1,
    nCondition = 2,
    chnl = "F1",
    fallbackHighValue = TRUE
  )

  fbeta_row <- res[res$method == "fbeta", , drop = FALSE]
  tailgate_row <- res[res$method == "tailgate", , drop = FALSE]

  expect_equal(fbeta_row$error[[1]], "fbeta boom")
  expect_equal(
    fbeta_row$gateReturnPoint[[1]],
    "fbeta_error_fallback_high_value"
  )
  expect_true(isTRUE(fbeta_row$thresholdFallbackUsed[[1]]))
  expect_true(is.finite(fbeta_row$propRespEst[[1]]))

  expect_true(is.na(tailgate_row$error[[1]]))
  expect_equal(tailgate_row$gateReturnPoint[[1]], "tailgate_calculated")
})

test_that("mismatch error summaries average scenario statistics equally", {
  env <- .compare_plot_env()
  # Two scenarios (cell counts) of one transformation and mismatch size. The
  # first has samples with relative errors 0.1 and 0.3 (truth 1); the second
  # has a single sample with relative error 0.5; a zero-truth sample is
  # ignored.
  raw <- tibble::tibble(
    base_scenario_id = c(1L, 1L, 2L, 2L),
    transformation = "gaussian", mean_pos_setting = "high",
    prob_response = 0.1, n_cell = c(100, 100, 1000, 1000),
    mismatch_type = "sd_inflation", mismatch_val = 0.1, method = "stimgate",
    propRespTruth = c(1, 1, 1, 0),
    propRespEst = c(1.1, 0.7, 1.5, 0.2)
  )
  scen <- env$.simCompareUnsignedErrorSummary(
    raw, scenarioCols = c("base_scenario_id", "n_cell", "transformation",
      "mean_pos_setting", "prob_response", "mismatch_type", "mismatch_val",
      "method")
  )
  expect_equal(scen$median, c(0.2, 0.5))
  expect_equal(scen$max, c(0.3, 0.5))
  expect_equal(scen$q95, c(0.1 + 0.95 * 0.2, 0.5))
  # Equal weight per scenario, not per sample (which would give 0.3).
  avg <- env$.simCompareErrorAverage(
    scen, c("transformation", "mismatch_type", "mismatch_val", "method")
  )
  expect_equal(nrow(avg), 1L)
  expect_equal(avg$median, 0.35)
  expect_equal(avg$max, 0.4)
  # Averaging over response probabilities only keeps cell counts apart.
  by_cell <- env$.simCompareErrorAverage(
    scen, c("transformation", "mismatch_type", "mismatch_val", "method",
      "n_cell")
  )
  expect_equal(by_cell$median, c(0.2, 0.5))

  signed <- env$.simCompareSignedErrorSummary(
    raw, scenarioCols = c("base_scenario_id", "n_cell", "transformation",
      "mismatch_type", "mismatch_val", "method")
  )
  signed_avg <- env$.simBandwidthSignedErrorAverage(
    signed, c("transformation", "mismatch_type", "mismatch_val", "method")
  )
  over <- signed_avg[signed_avg$direction == "over", ]
  # Scenario 1: over 0.1 (share 1/2); scenario 2: over 0.5 (share 1).
  expect_equal(over$median, 0.3)
  expect_equal(over$prop, 0.75)

  for (by_prob in c(FALSE, TRUE)) {
    plot_unsigned <- env$.simComparePlotMismatchError(scen, by_prob = by_prob)
    expect_no_error(ggplot2::ggplotGrob(plot_unsigned))
    expect_setequal(
      as.character(plot_unsigned$data$statistic),
      c("Median", "95th percentile", "Maximum")
    )
    plot_signed <- env$.simComparePlotSignedError(
      signed, x = "mismatch_val", x_log = FALSE, by_prob = FALSE
    )
    expect_no_error(ggplot2::ggplotGrob(plot_signed))
  }
})

test_that("figure loop combines several extra columns in one heading", {
  env <- .compare_plot_env()
  dir <- tempfile("compare-fig-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  data <- tidyr::expand_grid(
    method = c("stimgate", "tailgate"), mean_pos_setting = "high",
    mismatch_type = c("sd_inflation", "mean_shift_all"), n_cell = c(100, 1000),
    mismatch_val = c(0, 1)
  ) |>
    dplyr::mutate(value = seq_len(dplyr::n()))
  out <- utils::capture.output(
    env$.simCompareFigureLoop(
      data,
      make_plot = function(d) {
        ggplot2::ggplot(d, ggplot2::aes(mismatch_val, value, colour = method)) +
          ggplot2::geom_line()
      },
      dir = dir,
      file_fn = function(pos, extra) {
        paste0(extra$mismatch_type, "_", extra$n_cell, ".png")
      },
      extra_col = c("mismatch_type", "n_cell"),
      extra_heading = function(extra) {
        paste0(extra$mismatch_type, "; cells: ", extra$n_cell)
      },
      height = 6, level = 4L
    )
  )
  expect_true(any(out == "###### mean_shift_all; cells: 100"))
  expect_true(file.exists(file.path(dir, "no_tailgate", "sd_inflation_1000.png")))
})

test_that(".simCompareFreqBsGrid writes the bandwidth-style progress summary", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(
    file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R"),
    local = env
  )
  source(script_comp, local = env)

  tmp_dir <- tempfile("compare-progress-")
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  dir_cache <- file.path(tmp_dir, "output")
  dir_jobs <- file.path(tmp_dir, "jobs")
  path_progress <- file.path(tmp_dir, "progress.txt")

  grid <- data.frame(
    sim_id = c(1L, 2L),
    transformation = c("gaussian", "gaussian"),
    mean_pos = c(5, 5),
    prob_response = c(0.1, 0.1),
    n_cell = c(100, 100),
    bias_uns = c(0, 0),
    bw = c(0.1, 0.1),
    sample_perturbation_sd = c(0, 0),
    condition_perturbation_sd = c(0, 0),
    cluster_perturbation_sd = c(0, 0),
    background_relative_to_response = c(0.1, 0.1),
    n_cell_uns_relative_to_stim = c(1, 1),
    stringsAsFactors = FALSE
  )
  run_grid <- function() {
    env$.simCompareFreqBsGrid(
      sim_grid = grid, nSample = 1, nIter = 1, nMarker = 1,
      nCondition = 2, nCluster = 2, probExact = TRUE,
      tailgateAutoTol = TRUE, dirCache = dir_cache,
      pathProgress = path_progress, dirJobs = dir_jobs,
      progressHeading = "TEST PROGRESS",
      resume = TRUE, parallel = FALSE, progress = FALSE
    )
  }

  run_grid()
  expect_setequal(list.files(dir_jobs), c("completed-1", "completed-2"))
  summary_lines <- readLines(path_progress, warn = FALSE)
  expect_true(any(grepl("TEST PROGRESS", summary_lines, fixed = TRUE)))
  expect_true(any(grepl("Total Simulations  : 2", summary_lines, fixed = TRUE)))
  expect_true(any(grepl("Completed (Success): 2", summary_lines, fixed = TRUE)))
  expect_false(any(grepl("^Running: ", summary_lines)))

  # Resumed rows keep their completed markers and leave none running.
  run_grid()
  expect_setequal(list.files(dir_jobs), c("completed-1", "completed-2"))
  expect_true(any(grepl(
    "In Progress        : 0", readLines(path_progress, warn = FALSE),
    fixed = TRUE
  )))
})

test_that("signed-error plot draws varying-width lines for methods and directions", {
  env <- .compare_plot_env()
  raw <- tidyr::expand_grid(
    n_cell = c(1000, 10000, 100000),
    method = c("stimgate", "fbeta"),
    rep = 1:6
  ) |>
    dplyr::mutate(
      transformation = "gaussian",
      propRespTruth = 0.01,
      propRespEst = 0.01 * (1 + c(-0.5, -0.2, 0.1, 0.3, 0.8, 1.5)[rep])
    )
  tbl <- env$.simCompareSignedErrorSummary(
    raw, scenarioCols = c("transformation", "n_cell", "method")
  )
  expect_setequal(tbl$direction, c("over", "under"))
  plot <- env$.simComparePlotSignedError(tbl)
  expect_no_error(ggplot2::ggplotGrob(plot))
  expect_setequal(
    unique(as.character(plot$data$statistic)),
    c("Median", "95th percentile", "Maximum")
  )
  segment_layers <- Filter(
    function(l) inherits(l$geom, "GeomSegment"), plot$layers
  )
  expect_length(segment_layers, 1L)
  expect_true("linewidth" %in% names(segment_layers[[1]]$mapping))
  expect_gt(length(unique(segment_layers[[1]]$data$prop_segment)), 1L)
  expect_true(all(c("stimgate", "fbeta") %in% plot$data$method))
})

test_that("figure loop writes each method set to its own folder with headings", {
  env <- .compare_plot_env()
  dir <- tempfile("compare-fig-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  data <- tidyr::expand_grid(
    method = c("stimgate", "tailgate", "fbeta"),
    mean_pos_setting = c("low", "high"),
    mismatch_val = c(0, 1)
  ) |>
    dplyr::mutate(value = seq_len(dplyr::n()))
  make_plot <- function(d) {
    ggplot2::ggplot(d, ggplot2::aes(mismatch_val, value, colour = method)) +
      ggplot2::geom_line() +
      env$.analysis_scale_method()
  }
  out <- utils::capture.output(
    env$.simCompareFigureLoop(
      data, make_plot, dir = dir,
      file_fn = function(pos, extra) paste0("fig_", pos, ".png"),
      height = 6, level = 4L
    )
  )
  expect_true(any(out == "#### All methods"))
  expect_true(any(out == "#### Without Tailgate"))
  expect_true(any(out == "##### Mean position: low"))
  for (set in c("all_methods", "no_tailgate")) {
    expect_true(file.exists(file.path(dir, set, "fig_low.png")))
    expect_true(file.exists(file.path(dir, set, "fig_high.png")))
  }
})

test_that("a missing bias_uns lets StimGate set the bias from its bandwidth", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)
  env$.simCompareEnsureCurrentCheckout <- function(...) invisible(TRUE)
  captured <- new.env()
  env$.simCompareFreqBs <- function(...) {
    args <- list(...)
    captured$has_bias <- "biasUns" %in% names(args)
    captured$bias <- args$biasUns
    captured$factor <- args$biasUnsFactor
    captured$scope <- args$bwScope
    stop("stop after capturing arguments")
  }
  row <- data.frame(
    sim_id = 1L, transformation = "gaussian", mean_pos = 5,
    prob_response = 0.1, n_cell = 100, bias_uns = NA_real_,
    bias_uns_factor = 4, bw_mtd = "nrd0"
  )
  suppressWarnings(env$.simCompareRunScenario(row, nSample = 1, nIter = 1))
  expect_true(captured$has_bias)
  expect_null(captured$bias)
  expect_identical(captured$factor, 4)
  expect_identical(captured$scope, "cytokine")

  row$bw_scope <- "sample"
  suppressWarnings(env$.simCompareRunScenario(row, nSample = 1, nIter = 1))
  expect_identical(captured$scope, "sample")

  row$bias_uns <- 0.15
  suppressWarnings(env$.simCompareRunScenario(row, nSample = 1, nIter = 1))
  expect_identical(captured$bias, 0.15)
})

test_that("comparison error plots cap errors at 16 times the truth", {
  env <- .compare_plot_env()
  scen <- tidyr::expand_grid(
    transformation = "gaussian", method = c("stimgate", "fbeta"),
    mismatch_val = c(0, 0.1)
  ) |>
    dplyr::mutate(median = c(0.1, 30, 0.2, 0.4), q95 = 2 * median, max = 3 * median)
  unsigned <- env$.simComparePlotMismatchError(scen)
  expect_equal(max(unsigned$data$value_shown), 15)
  expect_equal(max(unsigned$data$value), 90)
  expect_no_error(ggplot2::ggplotGrob(unsigned))
  y_labels <- ggplot2::layer_scales(unsigned)$y$get_labels()
  expect_true(any(grepl("≥ 1,500%", y_labels, fixed = TRUE)))

  signed <- scen |>
    dplyr::mutate(direction = "over", prop = 1) |>
    dplyr::bind_rows(
      scen |> dplyr::mutate(
        direction = "under", prop = 0.5,
        dplyr::across(c(median, q95, max), ~ -pmin(.x, 0.9))
      )
    )
  plot_signed <- env$.simComparePlotSignedError(
    signed, x = "mismatch_val", x_log = FALSE
  )
  expect_equal(max(plot_signed$data$value_shown), 15)
  expect_no_error(ggplot2::ggplotGrob(plot_signed))
  expect_equal(env$.simBandwidthAbsErrorLabel(c(0.5, 15), cap = 15), c("50%", "≥ 1,500%"))
})
