root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

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

test_that(".simCompareSimCytExperiment keeps cluster mismatch local and supports old and new simcyto APIs", {
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

  legacy_seen <- NULL
  testthat::with_mocked_bindings(
    simCytExperiment = function(...,
                                stimMeanShift = 0,
                                stimSdMultiplier = 1,
                                stimMeanShiftClusters = NULL,
                                stimSdMultiplierClusters = NULL) {
      legacy_seen <<- list(
        stimMeanShiftClusters = stimMeanShiftClusters,
        stimSdMultiplierClusters = stimSdMultiplierClusters,
        stimMeanShift = stimMeanShift,
        stimSdMultiplier = stimSdMultiplier
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

      expect_equal(legacy_seen$stimMeanShiftClusters, "gn")
      expect_equal(legacy_seen$stimSdMultiplierClusters, "gp")
      expect_equal(
        as.vector(flowCore::exprs(out_legacy[["flowFrameList"]][[2L]])),
        c(3, 1.5, 5, 4.5),
        tolerance = 1e-8
      )
    }
  )

  new_seen <- NULL
  testthat::with_mocked_bindings(
    simCytExperiment = function(...,
                                stimMeanShift = 0,
                                stimSdMultiplier = 1) {
      new_seen <<- list(...)
      make_out()
    },
    .package = "simcyto",
    {
      out_new <- env$.simCompareSimCytExperiment(
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

      expect_false("stimMeanShiftClusters" %in% names(new_seen))
      expect_false("stimSdMultiplierClusters" %in% names(new_seen))
      expect_equal(
        as.vector(flowCore::exprs(out_new[["flowFrameList"]][[2L]])),
        c(3, 1.5, 5, 4.5),
        tolerance = 1e-8
      )
    }
  )
})

test_that(".simCompareFreqBs with zero mismatch reproduces clean baseline", {
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

    sign_check <- tibble::tibble(
      method = rep("stimgate", 6L),
      approach = rep("stimgate", 6L),
      error = rep("", 6L),
      propRespEst = c(0.02, 0.1, 0.16, 0.22, 0.36, 0.44),
      propRespTruth = rep(0.2, 6L),
      propStim = rep(0.2, 6L),
      propUns = rep(0.1, 6L),
      threshold = rep(1, 6L),
      thresholdFallbackUsed = rep(FALSE, 6L),
      mismatch_type = rep("mean_shift", 6L),
      mismatch_val = c(0, 0.01, 0.02, 0.03, 0.04, 0.05)
    )

    sign_summary <- env$.simCompareSummariseFreqBs(
      sign_check,
      scenarioCols = "method"
    )

    expect_lt(sign_summary$med_abs_rel_error[[1]], 0)
    expect_gt(sign_summary$max_abs_rel_error[[1]], 0)
  }
)

test_that(
  ".simCompareFreqBsGrid parallel and serial runs produce equivalent results",
  {
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
  expect_true(grepl('c\\("sim",\\s*"compare",\\s*"freq_bs_batch"\\)', content8) ||
    grepl('"log"[^)]*"freq_bs_batch"', content8))
  # Analysis 8 progress log should be in freq_bs_batch log namespace
  expect_true(
    grepl('"log"[^)]*"freq_bs_batch"', content8) ||
      grepl('analysis_key\\s*=\\s*c\\([^)]*"freq_bs_batch"', content8)
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

test_that("Temporary project directories are unique and collision-safe", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  source(script_comp, local = env)

  paths <- replicate(20, {
    file.path(
      tempdir(),
      "stimgate-sim-compare",
      paste0(
        "pid-",
        Sys.getpid(),
        "-iter-",
        1,
        "-",
        format(Sys.time(), "%Y%m%d%H%M%OS6"),
        "-",
        sample.int(1e9, 1)
      )
    )
  })

  expect_equal(length(unique(paths)), length(paths))
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
