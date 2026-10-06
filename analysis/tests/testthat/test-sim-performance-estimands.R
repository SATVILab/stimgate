.performance_test_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c(
    "analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R"
  )) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

.performance_fixture <- function(n_dataset = 20L, n_sample = 20L) {
  tidyr::expand_grid(iter = seq_len(n_dataset), sample = as.character(seq_len(n_sample))) |>
    dplyr::mutate(
      scenario = "a", sim_id = 1L, sim_seed = 12345L,
      method = "stimgate", transformation = "gaussian",
      mean_pos_setting = "high", prob_response = 0.01, n_cell = 1000,
      mismatch_type = "sd_inflation", mismatch_val = 0.1,
      propRespTruth = 0.01, propRespEst = 0.01,
      propStim = 0.02, propUns = 0.01, threshold = 1,
      thresholdFallbackUsed = FALSE, error = NA_character_
    )
}

.performance_set_errors <- function(data, rel_error) {
  data$propRespEst <- data$propRespTruth * (1 + rel_error)
  data
}

test_that("main performance percentiles pool tubes rather than averaging dataset percentiles", {
  env <- .performance_test_env()
  raw <- .performance_fixture()
  # Ten A datasets: eighteen exact errors and two 100x errors each;
  # ten B datasets: twenty 55x errors each. Type-7 pooled q95 is 57.25,
  # whereas the mean of the twenty within-dataset q95 values is 77.5.
  err <- ifelse(raw$iter <= 10L, ifelse(raw$sample %in% c("19", "20"), 100, 0), 55)
  raw <- .performance_set_errors(raw, err)
  cols <- c("scenario", "method")
  plain <- env$.simCompareUnsignedErrorSummary(raw, cols)
  with_intervals <- env$.simCompareUnsignedErrorSummary(raw, cols, mcse = TRUE)
  expect_equal(plain$median, 55)
  expect_equal(plain$q95, 57.25)
  expect_equal(mean(vapply(split(err, raw$iter),
    function(x) stats::quantile(x, 0.95, names = FALSE), numeric(1))), 77.5)
  expect_false(isTRUE(all.equal(plain$median,
    mean(vapply(split(err, raw$iter), stats::median, numeric(1))))))
  expect_equal(plain$median, with_intervals$median)
  expect_equal(plain$q95, with_intervals$q95)
  signed <- env$.simCompareSignedErrorSummary(raw, cols, mcse = TRUE)
  over <- signed[signed$direction == "over", ]
  under <- signed[signed$direction == "under", ]
  expect_equal(over$median, 55)
  expect_equal(over$prop, 0.55)
  expect_equal(under$prop, 0)
  expect_true(is.na(under$median))
  summary <- env$.simCompareSummariseFreqBs(raw, cols, mcse = TRUE)
  expect_equal(summary$med_abs_rel_error, plain$median)
  expect_equal(summary$q95_abs_rel_error, plain$q95)
})

test_that("dataset maximum occurrence and conditional severity have distinct estimands", {
  env <- .performance_test_env()
  raw <- .performance_fixture()
  err <- ifelse(raw$iter <= 10L, ifelse(raw$sample == "1", 0.4, 0.1), 0)
  raw <- .performance_set_errors(raw, err)
  plain <- env$.simCompareDatasetMaxSummary(raw, c("scenario", "method"))
  bounded <- env$.simCompareDatasetMaxSummary(raw, c("scenario", "method"), mcse = TRUE)
  over <- bounded[bounded$direction == "over", ]
  under <- bounded[bounded$direction == "under", ]
  expect_equal(over$n_dataset_total, 20L)
  expect_equal(over$n_eligible, 20L)
  expect_equal(over$n_incomplete, 0L)
  expect_equal(over$n_affected, 10L)
  expect_equal(over$occurrence, 0.5)
  expect_equal(over$severity, 0.4)
  expect_equal(under$occurrence, 0)
  expect_equal(under$n_affected, 0L)
  expect_true(is.na(under$severity))
  expect_equal(plain$occurrence, bounded$occurrence)
  expect_equal(plain$severity, bounded$severity)
  expect_true(over$interval_available_occurrence)
  expect_true(over$interval_available_severity)
  expect_true(is.na(under$severity_lower) && is.na(under$severity_upper))
  expect_equal(under$n_bootstrap_finite_severity, 0L)
  # Occurrence resamples all twenty datasets, not only the ten affected ones.
  expect_lt(over$occurrence_lower, 0.5)
  expect_gt(over$occurrence_upper, 0.5)
})

test_that("exact agreement belongs to neither error direction", {
  env <- .performance_test_env()
  raw <- .performance_fixture(5L)
  raw <- .performance_set_errors(raw, rep(c(-0.2, 0, 0.4, 0), 25L))
  signed <- env$.simCompareSignedErrorSummary(raw, c("scenario", "method"))
  expect_equal(signed$prop, rep(0.25, 2L))
  expect_equal(signed$median[signed$direction == "over"], 0.4)
  expect_equal(signed$median[signed$direction == "under"], -0.2)
  empty <- env$.simCompareDatasetMaxSummary(.performance_fixture(),
    c("scenario", "method"), mcse = TRUE)
  expect_true(all(empty$occurrence == 0))
  expect_true(all(empty$n_affected == 0))
  expect_true(all(is.na(empty$severity)))
  expect_true(all(empty$n_eligible == 20L))
})

test_that("fixed-size maxima exclude incomplete outcomes instead of calling them unaffected", {
  env <- .performance_test_env()
  raw <- .performance_set_errors(.performance_fixture(), 0.4)
  raw$propRespEst[raw$iter == 1L & raw$sample == "1"] <- NA_real_
  raw$error[raw$iter == 2L & raw$sample == "1"] <- "runtime error with finite diagnostic"
  raw <- raw[!(raw$iter == 3L & raw$sample == "1"), ]
  # Dataset four has twenty rows, but only nineteen distinct sample IDs.
  raw$sample[raw$iter == 4L & raw$sample == "1"] <- "2"
  raw$sample[raw$iter == 5L & raw$sample == "1"] <- NA_character_
  out <- env$.simCompareDatasetMaxSummary(raw, c("scenario", "method"))
  expect_true(all(out$n_dataset_total == 20L))
  expect_true(all(out$n_eligible == 15L))
  expect_true(all(out$n_incomplete == 5L))
  expect_equal(out$n_affected[out$direction == "over"], 15L)
  expect_equal(out$occurrence[out$direction == "over"], 1)
  expect_equal(out$severity[out$direction == "over"], 0.4)
  # Entirely absent datasets remain incomplete when the intended count is supplied.
  absent <- env$.simCompareDatasetMaxSummary(raw[raw$iter != 20L, ],
    c("scenario", "method"), expected_datasets = 20L)
  expect_true(all(absent$n_dataset_total == 20L))
  expect_true(all(absent$n_incomplete == 6L))
  expect_true(all(absent$n_eligible == 14L))
  # Show differing eligible cohorts for methods; don't force a shared complete-case pool.
  comparison <- dplyr::bind_rows(raw,
    dplyr::mutate(.performance_set_errors(.performance_fixture(), 0.2), method = "fbeta"))
  methods <- env$.simCompareDatasetMaxSummary(comparison, c("scenario", "method"))
  expect_true(all(methods$n_eligible[methods$method == "stimgate"] == 15L))
  expect_true(all(methods$n_eligible[methods$method == "fbeta"] == 20L))
})

test_that("rare-direction bootstrap coverage suppresses unreliable conditional intervals", {
  env <- .performance_test_env()
  raw <- .performance_fixture()
  raw <- .performance_set_errors(raw, ifelse(raw$iter == 1L, 0.4, 0))
  out <- env$.simCompareDatasetMaxSummary(raw, c("scenario", "method"), mcse = TRUE)
  over <- out[out$direction == "over", ]
  expect_equal(over$n_affected, 1L)
  expect_equal(over$severity, 0.4)
  expect_equal(over$n_eligible, 20L)
  expect_lt(over$bootstrap_coverage_severity, 0.95)
  expect_gt(over$n_bootstrap_finite_severity, 0L)
  expect_equal(over$bootstrap_coverage_severity,
    over$n_bootstrap_finite_severity / over$n_bootstrap_finite_occurrence)
  expect_false(over$interval_available_severity)
  expect_true(is.na(over$severity_lower) && is.na(over$severity_upper))
  expect_true(over$interval_available_occurrence)
  # Eligible datasets <5 suppress even otherwise finite degenerate bounds.
  few <- env$.simCompareDatasetMaxSummary(
    .performance_set_errors(.performance_fixture(4L), 0.4),
    c("scenario", "method"), mcse = TRUE
  )
  expect_true(all(is.na(few$occurrence_lower) & is.na(few$severity_lower)))
})

test_that("classification remains a distribution of tube proportions", {
  env <- .performance_test_env()
  raw <- .performance_fixture(5L, 2L) |>
    dplyr::mutate(
      nTruePos = ifelse(sample == "1", 1L, 0L),
      nFalsePos = ifelse(sample == "1", 0L, 99L),
      nFalseNeg = ifelse(sample == "1", 0L, 1L),
      nTrueNeg = ifelse(sample == "1", 99L, 0L)
    )
  off <- env$.simCompareClassificationSummary(raw, c("scenario", "method"))
  on <- env$.simCompareClassificationSummary(raw, c("scenario", "method"), mcse = TRUE)
  expect_equal(off$fdp_median, 0.5)
  expect_equal(off$sensitivity_median, 0.5)
  expect_equal(off$fpr_median, 0.5)
  expect_equal(off$fdp_q90, 1)
  expect_equal(off$sensitivity_q10, 0)
  expect_equal(off$fpr_q90, 1)
  expect_equal(on$fdp_median, off$fdp_median)
  expect_equal(on$sensitivity_median, off$sensitivity_median)
  expect_equal(on$fpr_median, off$fpr_median)
  # Cell-count weighting would give 99/100 instead of the median tube FDP.
  expect_false(isTRUE(all.equal(off$fdp_median,
    sum(raw$nFalsePos) / sum(raw$nTruePos + raw$nFalsePos))))
})


test_that("whole-dataset bootstrap retains every selected sample with multiplicity", {
  env <- .performance_test_env()
  withr::local_preserve_seed()
  set.seed(312)
  before <- .Random.seed
  x <- c(0, 0, 1, 10)
  unit <- c(1L, 1L, 2L, 3L)
  indices <- env$.analysis_mcse_bootstrap_indices(unit, "paired-biological", reps = 41L)
  expect_equal(dim(indices), c(3L, 41L))
  blocks <- split(x, unit)
  manual <- vapply(seq_len(ncol(indices)), function(b) {
    values <- unlist(blocks[indices[, b]], use.names = FALSE)
    stats::median(values)
  }, numeric(1))
  draws <- env$.analysis_mcse_block_draws(
    x, unit, stats::median, "paired-biological", reps = 41L
  )
  expect_equal(draws, manual)
  lengths <- env$.analysis_mcse_block_draws(
    x, unit, length, "paired-biological", reps = 41L
  )
  expect_equal(lengths, colSums(matrix(c(2, 1, 1)[indices], nrow = 3L)))
  expect_true(any(lengths > length(x)))
  expect_true(any(lengths < length(x)))
  # Repeated draws of one dataset must repeat its samples, not drop duplicates.
  repeated <- which(apply(indices, 2L, function(ids) anyDuplicated(ids) > 0L))
  expect_gt(length(repeated), 0L)
  expect_identical(.Random.seed, before)
  expect_identical(indices,
    env$.analysis_mcse_bootstrap_indices(unit, "paired-biological", reps = 41L))
})

test_that("bootstrap indices pair methods and deterministic mismatch settings", {
  env <- .performance_test_env()
  x <- rep(c(0.1, 0.2, 0.4, 0.8, 1.6), each = 20L)
  unit <- rep(1:5, each = 20L)
  same <- env$.analysis_mcse_bootstrap_indices(unit, 12345L, reps = 41L)
  expect_identical(same,
    env$.analysis_mcse_bootstrap_indices(unit, 12345L, reps = 41L))
  independent <- env$.analysis_mcse_bootstrap_indices(unit, 12346L, reps = 41L)
  expect_false(identical(same, independent))
  base <- env$.analysis_mcse_block_draws(x, unit, mean, 12345L, reps = 41L)
  shifted_method <- env$.analysis_mcse_block_draws(
    2 * x + 3, unit, mean, 12345L, reps = 41L
  )
  expect_equal(shifted_method, 2 * base + 3)
})

test_that("scenario averages recompute the full statistic with paired dataset draws", {
  env <- .performance_test_env()
  raw <- .performance_fixture(5L)
  raw <- .performance_set_errors(raw, 0.1 * raw$iter^2)
  pair <- dplyr::bind_rows(raw, dplyr::mutate(raw,
    scenario = "b", sim_id = 2L, mismatch_val = 0.2))
  cols <- c("scenario", "method", "mismatch_val", "sim_seed")
  components <- env$.simCompareUnsignedErrorSummary(pair, cols, mcse = TRUE)
  averaged <- env$.simComparePerformanceAverage(pair, cols, "method", mcse = TRUE)
  single <- components[components$scenario == "a", ]
  expect_equal(averaged$median, single$median)
  expect_equal(averaged$q95, single$q95)
  expect_equal(averaged$median_mcse, single$median_mcse)
  expect_equal(averaged$median_lower, single$median_lower)
  expect_equal(averaged$median_upper, single$median_upper)
  expect_false(isTRUE(all.equal(averaged$median_mcse,
    single$median_mcse / sqrt(2))))
  expect_equal(components$.boot_median[[1]], components$.boot_median[[2]])
  expect_equal(averaged$.boot_median[[1]], components$.boot_median[[1]])
  # This point estimate is a mean of scenario pooled percentiles, not one
  # percentile taken after pooling both scenarios' tubes together.
  unequal <- pair
  unequal <- unequal[unequal$scenario != "b" | unequal$sample %in% c("1", "2"), ]
  unequal$propRespEst[unequal$scenario == "b"] <- 0.01 * 11
  out <- env$.simComparePerformanceAverage(unequal, cols, "method", mcse = TRUE)
  component_stats <- env$.simCompareUnsignedErrorSummary(unequal, cols, mcse = TRUE)
  expect_equal(out$median, mean(component_stats$median))
  expect_equal(out$.boot_median[[1]],
    (component_stats$.boot_median[[1]] + component_stats$.boot_median[[2]]) / 2)
  all_errors <- abs((unequal$propRespEst - unequal$propRespTruth) / unequal$propRespTruth)
  expect_false(isTRUE(all.equal(out$median, stats::median(all_errors))))
  plain <- env$.simComparePerformanceAverage(unequal, cols, "method")
  expect_equal(out$median, plain$median)
  expect_equal(out$q95, plain$q95)
})

test_that("twenty-sample scientific caches reject ten-sample and old-estimand outputs", {
  env <- .performance_test_env()
  withr::local_envvar(ANALYSIS_EXPECTED_RUN_ID = NA)
  current <- withr::local_tempdir()
  file.create(file.path(current, "COMPLETE"))
  saveRDS(1L, file.path(current, "compare_raw.rds"))
  ctx <- list(analysis_key = c("sim", "compare"), current_dir = current,
    qmd_path = "analysis/7-sim-compare-freq_bs.qmd")
  required <- list(comparison_semantics_version = "corrected-comparison-v13",
    n_sample_sim = 20L, n_iter_sim = 20L, sim_size = "final")
  write_manifest <- function(settings) {
    saveRDS(list(analysis_key = ctx$analysis_key, params = settings),
      file.path(current, "manifest.rds"))
  }
  old_size <- required
  old_size$n_sample_sim <- 10L
  write_manifest(old_size)
  expect_error(env$.analysis_current_file(ctx, "compare_raw.rds", required),
    "n_sample_sim", fixed = TRUE)
  old_estimand <- required
  old_estimand$comparison_semantics_version <- "corrected-comparison-v10"
  write_manifest(old_estimand)
  expect_error(env$.analysis_current_file(ctx, "compare_raw.rds", required),
    "comparison_semantics_version", fixed = TRUE)
  write_manifest(required)
  expect_true(file.exists(env$.analysis_current_file(ctx, "compare_raw.rds", required)))
  raw <- dplyr::bind_rows(lapply(c("stimgate", "fbeta", "tailgate"), function(method) {
    dplyr::mutate(.performance_fixture(1L, 10L),
      method = method, nCellStim = 100L, nPosStim = 1L,
      nTruePos = 1L, nFalsePos = 0L, nFalseNeg = 0L, nTrueNeg = 99L,
      unsExprSum = 1)
  }))
  status <- env$.simCompareGridOutputStatus(raw,
    sim_grid = tibble::tibble(sim_id = 1L), nSample = 20L, nIter = 1L)
  expect_false(status$validation_ok)
  expect_equal(status$failed_ids, 1L)
})

test_that("dataset maxima retain repeated dataset draws and expose differing cohorts", {
  env <- .performance_test_env()
  raw <- .performance_fixture()
  raw <- .performance_set_errors(raw,
    ifelse(raw$iter <= 10 & raw$sample == "1", raw$iter / 10, 0))
  raw$propRespEst[raw$iter == 20 & raw$sample == "1"] <- NA_real_
  out <- env$.simCompareDatasetMaxSummary(raw, c("scenario", "method"),
    expected_datasets = 20L, mcse = TRUE)
  over <- out[out$direction == "over", ]
  # Dataset IDs are sorted as characters by the public resampling helper.
  ids <- sort(as.character(1:20))
  indices <- env$.analysis_mcse_bootstrap_indices(ids, "sim_seed:12345")
  eligible <- as.integer(ids) != 20L
  affected <- as.integer(ids) <= 10L
  magnitude <- as.integer(ids) / 10
  occurrence <- apply(indices, 2, function(i) mean(affected[i][eligible[i]]))
  severity <- apply(indices, 2, function(i) {
    selected <- affected[i] & eligible[i]
    if (any(selected)) mean(magnitude[i][selected]) else NA_real_
  })
  expect_true(any(apply(indices, 2, anyDuplicated) > 0L))
  expect_equal(over$.boot_occurrence[[1]], occurrence)
  expect_equal(over$.boot_severity[[1]], severity)
  # Equal counts need not imply equal scientific cohorts.
  other <- dplyr::mutate(.performance_fixture(), method = "fbeta")
  other <- .performance_set_errors(other, 0.2)
  other$propRespEst[other$iter == 19 & other$sample == "1"] <- NA_real_
  cohorts <- env$.simCompareDatasetMaxSummary(dplyr::bind_rows(raw, other),
    c("scenario", "method"), expected_datasets = 20L)
  expect_true(all(cohorts$n_eligible == 19L))
  expect_true(all(cohorts$cohort_matches_stimgate[cohorts$method == "stimgate"]))
  expect_false(any(cohorts$cohort_matches_stimgate[cohorts$method == "fbeta"]))
  invalid <- .performance_fixture()
  invalid$sample[invalid$iter == 1 & invalid$sample == "20"] <- "21"
  invalid_summary <- env$.simCompareDatasetMaxSummary(invalid,
    c("scenario", "method"), expected_datasets = 20L)
  expect_true(all(invalid_summary$n_incomplete == 1L))
})

test_that("independent scenario families use independent full-statistic draws", {
  env <- .performance_test_env()
  a <- .performance_fixture(5L)
  a <- .performance_set_errors(a, a$iter^2 / 10)
  b <- dplyr::mutate(a, scenario = "b", sim_id = 2L, sim_seed = 98765L)
  raw <- dplyr::bind_rows(a, b)
  cols <- c("scenario", "method", "sim_seed")
  components <- env$.simCompareUnsignedErrorSummary(raw, cols, mcse = TRUE)
  expect_false(identical(components$.boot_median[[1]], components$.boot_median[[2]]))
  avg <- env$.simComparePerformanceAverage(raw, cols, "method", mcse = TRUE)
  expected <- (components$.boot_median[[1]] + components$.boot_median[[2]]) / 2
  expect_equal(avg$.boot_median[[1]], expected)
  expect_equal(avg$median_mcse, stats::sd(expected))
})

test_that("scenario-average draws never substitute a changing finite subset", {
  env <- .performance_test_env()
  a <- rep(0.2, 999L)
  b <- rep(0.8, 999L)
  b[1:100] <- NA_real_
  components <- tibble::tibble(method = "stimgate", median = c(0.2, 0.8),
    n_dataset_median = 20L, .boot_median = list(a, b))
  out <- env$.analysis_mcse_bootstrap_average(components, "method", "median")
  expect_equal(out$median, 0.5)
  expect_true(all(is.na(out$.boot_median[[1]][1:100])))
  expect_equal(out$.boot_median[[1]][101:999], rep(0.5, 899L))
  expect_equal(out$n_bootstrap_finite_median, 899L)
  expect_equal(out$bootstrap_coverage_median, 899 / 999)
  expect_false(out$interval_available_median)
  expect_true(is.na(out$median_lower) && is.na(out$median_upper))
})
