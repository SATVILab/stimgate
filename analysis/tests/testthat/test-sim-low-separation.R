root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

.lowSepEnv <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (f in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "coexpression-gates.R", "sim-low-separation.R")) {
    source(file.path(root_dir, "scripts", "r", f), local = env)
  }
  env
}

# One dataset's scenario and row columns (main grid, first scenario).
.lowSepToyScen <- function(env, iter = 1L) {
  scen <- env$.simLowSepGrid(55700L)[1, ]
  scen$n_cell <- 100
  dplyr::mutate(scen, iter = iter, sim_id = iter, sim_seed = 1L)
}

# Confusion counts of one sample; the co-expression rules copy the
# antimode calls.
.lowSepToyCounts <- function(env, stim, uns, iter = 1L, sample = 1L) {
  copy_coex <- function(x) {
    anti <- x[x$gate_type == "antimode", ]
    dplyr::bind_rows(x, dplyr::mutate(anti, gate_type = "coex_screened"),
      dplyr::mutate(anti, gate_type = "coex_unscreened"))
  }
  rows <- function(tube, x) {
    x <- copy_coex(x)
    tibble::tibble(
      sample = sample, tube = tube, gate_type = x$gate_type,
      truth = x$truth, call = x$call, n = as.integer(x$n)
    )
  }
  body <- dplyr::bind_rows(rows("stim", stim), rows("uns", uns))
  dplyr::bind_cols(.lowSepToyScen(env, iter)[rep(1L, nrow(body)), ], body)
}

.lowSepToyStim <- function() {
  tibble::tribble(
    ~gate_type, ~truth, ~call, ~n,
    "ordinary", "TNF+IFNg+", "TNF+IFNg+", 4,
    "ordinary", "TNF+IFNg+", "TNF+IFNg-", 6,
    "ordinary", "TNF-IFNg-", "TNF-IFNg-", 90,
    "antimode", "TNF+IFNg+", "TNF+IFNg+", 10,
    "antimode", "TNF-IFNg-", "TNF-IFNg-", 88,
    "antimode", "TNF-IFNg-", "TNF+IFNg+", 2
  )
}

.lowSepToyUns <- function() {
  tibble::tribble(
    ~gate_type, ~truth, ~call, ~n,
    "ordinary", "TNF-IFNg-", "TNF-IFNg-", 100,
    "antimode", "TNF-IFNg-", "TNF-IFNg-", 100
  )
}

test_that("the four rules share one set of labels, colours and short labels", {
  env <- .lowSepEnv()
  rules <- c("ordinary", "antimode", "coex_screened", "coex_unscreened")
  expect_identical(names(env$.simLowSepGateTypes), rules)
  expect_identical(names(env$.simLowSepGateColours), rules)
  expect_identical(names(env$.simLowSepGateShort), rules)
  expect_identical(env$.simLowSepCoexScreen, c(coex_screened = TRUE, coex_unscreened = FALSE))
  expect_equal(anyDuplicated(env$.simLowSepGateColours), 0L)
  expect_false(any(grepl("refine|cytpos", c(names(env$.simLowSepGateTypes), env$.simLowSepGateTypes))))
})

test_that("the scenario grid keeps earlier IDs and seeds and appends new families", {
  env <- .lowSepEnv()
  grid <- env$.simLowSepGrid(55700L)
  expect_equal(nrow(grid), 40L)
  expect_equal(grid$scenario_id, 1:40)
  # The main grid, the first families and the appended families each keep
  # their own seed stream, so earlier scenarios keep their data.
  stream <- function(seed, n) as.integer(env$.analysis_with_seed(seed, sample.int(.Machine$integer.max, n)))
  expect_equal(grid$scenario_seed[1:8], stream(55700L, 8L))
  expect_equal(grid$scenario_seed[9:32], stream(55701L, 24L))
  expect_equal(grid$scenario_seed[33:40], stream(55702L, 8L))
  expect_true(all(grid$family[1:8] == "main"))
  expect_setequal(unique(grid$family[33:40]), c("weak_double", "null_background_coexpression"))
  expect_setequal(unique(grid$family), names(env$.simLowSepFamilies))
  expect_true(all(grid$mismatch_type == "none"))
  expect_equal(grid$base_scenario_id, grid$scenario_id)
  expect_identical(grid, env$.simLowSepGrid(55700L))
  expect_equal(anyDuplicated(grid$scenario_seed), 0L)
  expect_true(all(env$.simLowSepScenCols %in% names(grid)))
  # Every split sums to one.
  for (s in env$.simLowSepSplits) expect_equal(sum(s), 1)
})

test_that("dataset rows get stable IDs and seeds drawn from the scenario seed", {
  env <- .lowSepEnv()
  grid <- env$.simLowSepGrid(55700L)
  rows <- env$.simLowSepRowGrid(grid)
  expect_equal(nrow(rows), 40L * 20L)
  expect_equal(rows$sim_id, seq_len(800L))
  expect_equal(rows$sim_id, (rows$scenario_id - 1L) * 20L + rows$iter)
  first <- rows[rows$iter == 1L, ]
  expect_equal(first$sim_seed, first$scenario_seed)
  # Seeds depend only on the scenario, not on which scenarios are run.
  sub <- env$.simLowSepRowGrid(grid[c(5, 33), ])
  expect_equal(sub$sim_seed, rows$sim_seed[rows$scenario_id %in% c(5, 33)])
  expect_equal(sub$sim_id, rows$sim_id[rows$scenario_id %in% c(5, 33)])
  # Fewer datasets keep the same IDs and seeds for the datasets they run.
  few <- env$.simLowSepRowGrid(grid[c(5, 33), ], n_iter = 5L)
  keep <- rows$scenario_id %in% c(5, 33) & rows$iter <= 5L
  expect_equal(few$sim_id, rows$sim_id[keep])
  expect_equal(few$sim_seed, rows$sim_seed[keep])
  expect_error(env$.simLowSepRowGrid(grid[1, ], n_iter = 21L), "At most")
  expect_equal(anyDuplicated(rows$sim_seed), 0L)
  expect_equal(env$.simLowSepSizes("final"), list(n_dataset = 20L, n_sample = 12L))
  expect_equal(env$.simLowSepSizes("draft"), list(n_dataset = 5L, n_sample = 12L))
  expect_equal(env$.simLowSepSizes("final", quick = TRUE), list(n_dataset = 2L, n_sample = 3L))
})

test_that("the batch grid shares each baseline's data across mismatch settings", {
  env <- .lowSepEnv()
  grid_a <- env$.simLowSepGrid(55700L)
  grid_b <- env$.simLowSepBatchGrid(55700L)
  mismatch <- env$.simLowSepBatchMismatch()
  expect_equal(nrow(mismatch), 5L)
  expect_equal(sort(mismatch$stim_mean_shift[mismatch$mismatch_type == "mean_shift"]), c(0, 0.025, 0.05, 0.2))
  expect_equal(mismatch$stim_sd_multiplier[mismatch$mismatch_type == "sd_inflation"], 1.1)
  base <- env$.simLowSepBatchBaseline(grid_a)
  expect_equal(nrow(grid_b), nrow(base) * 5L)
  expect_equal(grid_b$scenario_id, seq_len(nrow(grid_b)))
  expect_true(all(grid_b$separation == "low"))
  expect_setequal(unique(grid_b$n_cell), c(1e5, 5e3))
  # Baseline seeds are Analysis 11a's.
  expect_equal(grid_b$scenario_seed, grid_a$scenario_seed[match(grid_b$base_scenario_id, grid_a$scenario_id)])
  rows <- env$.simLowSepRowGrid(grid_b)
  seeds <- dplyr::summarise(dplyr::group_by(rows, .data$base_scenario_id, .data$iter),
    n_seed = dplyr::n_distinct(.data$sim_seed), .groups = "drop")
  expect_true(all(seeds$n_seed == 1L))
  expect_equal(anyDuplicated(rows$sim_id), 0L)
  # No-mismatch rows reproduce Analysis 11a's data.
  rows_a <- env$.simLowSepRowGrid(grid_a)
  none <- rows[rows$stim_mean_shift == 0 & rows$stim_sd_multiplier == 1, ]
  expect_equal(none$sim_seed,
    rows_a$sim_seed[match(paste(none$base_scenario_id, none$iter), paste(rows_a$scenario_id, rows_a$iter))])
  expect_length(env$.simLowSepMismatchArgs(none[1, ]), 0L)
})

test_that("the mismatch changes only the stimulated tube's all-negative component", {
  testthat::skip_if_not_installed("simcyto")
  testthat::skip_if_not(all(c("stimMeanShiftClusters", "stimSdMultiplierClusters") %in%
    names(formals(simcyto::simCytExperiment))))
  env <- .lowSepEnv()
  withr::local_preserve_seed()
  settings <- env$.simLowSepMainSettings()
  grid <- env$.simLowSepBatchGrid(55700L)
  pick <- function(shift, sd) {
    row <- grid[grid$family == "main" & grid$n_cell == 5e3 & grid$response_level == "higher" &
      grid$stim_mean_shift == shift & grid$stim_sd_multiplier == sd, ]
    row$n_cell <- 400
    row
  }
  sim <- function(row) env$.analysis_with_seed(1L, env$.simLowSepSimulate(row, 1L, settings))
  ex <- function(res, i) flowCore::exprs(res$flowFrameList[[i]])
  base <- sim(pick(0, 1))
  shifted <- sim(pick(0.2, 1))
  wider <- sim(pick(0, 1.1))
  expect_identical(env$.simLowSepMismatchArgs(pick(0.2, 1))$stimMeanShiftClusters, "TNF-IFNg-")
  for (res in list(shifted, wider)) {
    expect_identical(res$labelsList, base$labelsList)
    # The unstimulated tube is unchanged.
    expect_equal(ex(res, 1L), ex(base, 1L))
  }
  neg <- base$labelsList[[2L]] == "TNF-IFNg-"
  expect_true(any(neg) && any(!neg))
  # Positive stimulated cells are unchanged; negatives move by 0.2 on both markers.
  expect_equal(ex(shifted, 2L)[!neg, ], ex(base, 2L)[!neg, ])
  expect_equal(unname(ex(shifted, 2L)[neg, ] - ex(base, 2L)[neg, ]), matrix(0.2, sum(neg), 2L),
    tolerance = 1e-5)
  expect_equal(ex(wider, 2L)[!neg, ], ex(base, 2L)[!neg, ])
  centred <- function(m) sweep(m, 2L, colMeans(m))
  expect_equal(centred(ex(wider, 2L)[neg, ]), 1.1 * centred(ex(base, 2L)[neg, ]), tolerance = 1e-5)
})

test_that("only the Gaussian transformation is accepted", {
  env <- .lowSepEnv()
  settings <- env$.simLowSepMainSettings()
  settings$transformation <- "gamma"
  row <- env$.simLowSepGrid(1L)[1, ]
  expect_error(env$.simLowSepSimulate(row, 1L, settings), "Gaussian")
})

test_that("antimode and co-expression calls follow their rules", {
  env <- .lowSepEnv()
  gate <- c(tnf = 4, ifng = 4, tnf_cyt = 2, ifng_cyt = 2)
  x_tnf <- c(5, 1, 3, 5, 3)
  x_ifng <- c(3, 3, 3, 5, 5)
  expect_equal(
    env$.simLowSepCall(x_tnf, x_ifng, gate, "ordinary"),
    c("TNF+IFNg-", "TNF-IFNg-", "TNF-IFNg-", "TNF+IFNg+", "TNF-IFNg+")
  )
  # IFNg between the gates counts only alongside ordinary TNF positivity, and
  # two cytokines that are both only above their lower gates stay negative.
  expect_equal(
    env$.simLowSepCall(x_tnf, x_ifng, gate, "antimode"),
    c("TNF+IFNg+", "TNF-IFNg-", "TNF-IFNg-", "TNF+IFNg+", "TNF+IFNg+")
  )
  # The package's positivity helper gives the same calls.
  ex <- data.frame(F1 = x_tnf, F2 = x_ifng)
  gate_tbl <- data.frame(chnl = c("F1", "F2"), gate = c(4, 4), gateCyt = c(2, 2))
  pos <- stimgate:::.getPosIndByChnl(ex, gate_tbl, c("F1", "F2"), "cyt")
  expect_equal(
    env$.simLowSepCombn[1L + pos$F1 + 2L * pos$F2],
    env$.simLowSepCall(x_tnf, x_ifng, gate, "antimode")
  )
  # A lowered IFNg gate (2.5) for cells above a raised TNF cut (4.5).
  low <- data.frame(a = "TNF", b = "IFNg", cut = 2.5, condCut = 4.5, lowered = TRUE)
  expect_equal(
    env$.simLowSepCall(x_tnf, x_ifng, gate, "coex_unscreened", low),
    c("TNF+IFNg+", "TNF-IFNg-", "TNF-IFNg-", "TNF+IFNg+", "TNF-IFNg+")
  )
  expect_equal(
    env$.simLowSepCall(x_tnf, x_ifng, gate, "coex_screened", NULL),
    env$.simLowSepCall(x_tnf, x_ifng, gate, "ordinary")
  )
  expect_error(env$.simLowSepCall(1, 1, gate, "other"), "gate_type")
  expect_error(env$.simLowSepCall(1, 1, gate, "cytpos"), "gate_type")
})

test_that("recovering IFNg in TNF+ cells corrects the exclusive single-positive overestimate", {
  env <- .lowSepEnv()
  # Stimulated tube: 10 true TNF+IFNg+ cells; the ordinary gates call 6 of them
  # TNF+IFNg- only, the antimode gates recover them all.
  counts <- .lowSepToyCounts(env, .lowSepToyStim(), .lowSepToyUns())
  freq <- env$.simLowSepFrequencies(counts)
  get <- function(gt, q, col) freq[[col]][freq$gate_type == gt & freq$quantity == q]
  expect_equal(get("ordinary", "TNF+IFNg-", "prop_bs_est"), 0.06)
  expect_equal(get("antimode", "TNF+IFNg-", "prop_bs_est"), 0)
  expect_equal(get("ordinary", "TNF+IFNg-", "prop_bs_truth"), 0)
  expect_equal(get("ordinary", "TNF+", "prop_bs_est"), get("ordinary", "TNF+", "prop_bs_truth"))
  expect_equal(get("antimode", "IFNg+", "prop_bs_est"), 0.12)
  expect_equal(get("antimode", "TNF+IFNg+", "error_bs"), 0.02)
  expect_equal(get("antimode", "IFNg+", "n_est_stim"), 12L)
  expect_equal(get("antimode", "IFNg+", "n_tube_uns"), 100L)
  expect_setequal(unique(freq$quantity_type), c("marginal", "exclusive"))

  rec <- env$.simLowSepRecovery(counts)
  ifng <- rec[rec$marker == "IFNg" & rec$tube == "stim", ]
  expect_equal(ifng$sensitivity[ifng$gate_type == "ordinary"], 0.4)
  expect_equal(ifng$sensitivity[ifng$gate_type == "antimode"], 1)
  expect_equal(ifng$n_contaminating[ifng$gate_type == "antimode"], 2L)
  expect_equal(ifng$fdp[ifng$gate_type == "antimode"], 2 / 12)

  added <- env$.simLowSepAdded(rec)
  a <- added[added$marker == "IFNg" & added$tube == "stim" & added$gate_type == "antimode", ]
  expect_equal(c(a$n_added, a$n_added_true, a$n_added_negative), c(8L, 6L, 2L))
  expect_false("ordinary" %in% added$gate_type)

  summ <- env$.simLowSepFrequencySummary(freq)
  row <- summ[summ$quantity == "TNF+IFNg-" & summ$gate_type == "antimode", ]
  expect_equal(row$n_closer, 1L)
  expect_equal(row$n_further + row$n_unchanged, 0L)
  expect_true(is.na(summ$n_closer[summ$quantity == "TNF+IFNg-" & summ$gate_type == "ordinary"]))
})

test_that("cell-level F1 counts true positives within each quantity", {
  env <- .lowSepEnv()
  f1 <- env$.simLowSepF1(.lowSepToyCounts(env, .lowSepToyStim(), .lowSepToyUns()))
  get <- function(gt, q, col) f1[[col]][f1$gate_type == gt & f1$quantity == q]
  expect_equal(get("ordinary", "TNF+IFNg+", "f1"), 2 * 4 / (4 + 10))
  expect_equal(get("antimode", "TNF+IFNg+", "f1"), 2 * 10 / (12 + 10))
  expect_equal(get("antimode", "TNF+IFNg+", "precision"), 10 / 12)
  # No true and no called cells: undefined, not zero.
  expect_true(is.na(get("antimode", "TNF+IFNg-", "f1")))
  pct <- env$.simLowSepPercentiles(f1, "f1")
  expect_equal(pct$p50[pct$gate_type == "antimode" & pct$quantity == "TNF+IFNg+"], 10 / 11)
})

test_that("the one-sided Fisher test equals fisher.test", {
  env <- .lowSepEnv()
  tabs <- list(c(12, 100, 0, 100), c(5, 1000, 3, 900), c(0, 50, 4, 60), c(30, 5000, 10, 5000))
  for (t in tabs) {
    expected <- stats::fisher.test(matrix(c(t[1], t[3], t[2] - t[1], t[4] - t[3]), 2),
      alternative = "greater")$p.value
    expect_equal(env$.simLowSepFisherGreater(t[1], t[2], t[3], t[4]), expected)
  }
  # Vectorised.
  expect_length(env$.simLowSepFisherGreater(c(1, 2), c(10, 10), c(0, 0), c(10, 10)), 2L)
})

test_that("a missing gate calls every cell negative under every rule", {
  env <- .lowSepEnv()
  gate <- c(tnf = NA_real_, ifng = 1, tnf_cyt = NA_real_, ifng_cyt = 0.5)
  for (gt in names(env$.simLowSepGateTypes)) {
    call <- env$.simLowSepCall(c(2, 0, 3), c(2, 0.7, 0), gate, gt)
    expect_false(anyNA(call), info = gt)
    expect_false(any(grepl("TNF\\+", call)), info = gt)
  }
})

test_that("paired tests treat no difference as no evidence", {
  env <- .lowSepEnv()
  expect_equal(env$.simLowSepPairedP(c(0, 0, 0, 0), "wilcoxon"), 1)
  expect_equal(env$.simLowSepPairedP(c(0, 0, 0, 0), "t"), 1)
  expect_equal(env$.simLowSepPairedP(c(0.1, 0.1, 0.1), "t"), 0)
  # Constant up to rounding: t.test() would stop.
  expect_equal(env$.simLowSepPairedP(0.1 + c(0, 1e-17, 2e-17, -1e-17), "t"), 0)
  expect_true(is.na(env$.simLowSepPairedP(0.1, "wilcoxon")))
  d <- c(0.01, 0.02, -0.005, 0.03, 0.015, 0.04)
  expect_equal(env$.simLowSepPairedP(d, "wilcoxon"),
    stats::wilcox.test(d, alternative = "greater")$p.value)
  expect_equal(env$.simLowSepPairedP(d, "t"), stats::t.test(d, alternative = "greater")$p.value)
})

test_that("rejection rates count samples and datasets with their intervals", {
  env <- .lowSepEnv()
  # Two datasets of two samples: the stimulated tube has many more positives
  # in dataset 1 only.
  make <- function(iter, sample, stim_pos) {
    stim <- tibble::tribble(
      ~gate_type, ~truth, ~call, ~n,
      "ordinary", "TNF-IFNg-", "TNF+IFNg+", stim_pos,
      "ordinary", "TNF-IFNg-", "TNF-IFNg-", 1000 - stim_pos,
      "antimode", "TNF-IFNg-", "TNF+IFNg+", stim_pos,
      "antimode", "TNF-IFNg-", "TNF-IFNg-", 1000 - stim_pos
    )
    uns <- tibble::tribble(
      ~gate_type, ~truth, ~call, ~n,
      "ordinary", "TNF-IFNg-", "TNF+IFNg+", 1,
      "ordinary", "TNF-IFNg-", "TNF-IFNg-", 999,
      "antimode", "TNF-IFNg-", "TNF+IFNg+", 1,
      "antimode", "TNF-IFNg-", "TNF-IFNg-", 999
    )
    .lowSepToyCounts(env, stim, uns, iter = iter, sample = sample)
  }
  counts <- dplyr::bind_rows(make(1L, 1L, 40), make(1L, 2L, 30), make(2L, 1L, 1), make(2L, 2L, 0))
  freq <- env$.simLowSepFrequencies(counts)
  sample_tests <- env$.simLowSepSampleTests(freq)
  dataset_tests <- env$.simLowSepDatasetTests(freq)
  expected <- env$.simLowSepExpectedNet(.lowSepToyScen(env))
  rates <- env$.simLowSepRejectionRates(sample_tests, dataset_tests, expected, mcse = TRUE)
  dp <- rates[rates$quantity == "TNF+IFNg+" & rates$gate_type == "ordinary", ]
  fisher <- dp[dp$test == "fisher", ]
  expect_equal(fisher$rate, 0.5)
  expect_equal(fisher$n_unit, 4L)
  # Fewer than five datasets: no bootstrap interval.
  expect_false(fisher$interval_available)
  expect_true(is.na(fisher$rate_lower))
  t_row <- dp[dp$test == "t", ]
  expect_equal(t_row$n_unit, 2L)
  ci <- stats::binom.test(t_row$n_reject, 2L)$conf.int
  expect_equal(c(t_row$rate_lower, t_row$rate_upper), as.numeric(ci))
  expect_setequal(rates$test, c("fisher", "wilcoxon", "t"))
  expect_equal(unique(dp$response_truth), "response")
  off <- env$.simLowSepRejectionRates(sample_tests, dataset_tests, expected, mcse = FALSE)
  expect_equal(off$rate, rates$rate)
  expect_true(all(is.na(off$rate_lower)))
})

test_that("expected net frequencies say where a rejection is a false positive", {
  env <- .lowSepEnv()
  grid <- env$.simLowSepGrid(55700L)
  net <- env$.simLowSepExpectedNet(grid)
  get <- function(fam, q, level) {
    ids <- grid$scenario_id[grid$family == fam & grid$response_level == level]
    unique(round(net$net_expected[net$scenario_id %in% ids & net$quantity == q], 12))
  }
  expect_equal(get("main", "TNF+IFNg+", "lower"), 0.01 * 0.6 * 0.8)
  expect_equal(get("no_coexpression", "TNF+IFNg+", "higher"), 0)
  expect_equal(get("weak_double", "TNF+IFNg+", "lower"), round(0.01 * 0.06 * 0.8, 12))
  expect_equal(get("weak_double", "TNF+", "higher"), round(0.05 * 0.53 * 0.8, 12))
  for (fam in c("null", "null_background_coexpression")) {
    ids <- grid$scenario_id[grid$family == fam]
    expect_true(all(net$response_truth[net$scenario_id %in% ids] == "none"))
  }
  expect_setequal(unique(net$quantity), env$.simLowSepQuantities$quantity)
})

test_that("validation and collation check and unnest each dataset's tables", {
  env <- .lowSepEnv()
  scen <- .lowSepToyScen(env)
  n_sample <- 1L
  counts <- .lowSepToyCounts(env, .lowSepToyStim(), .lowSepToyUns())
  counts <- tidyr::complete(dplyr::select(counts, "sample", "tube", "gate_type", "truth", "call", "n"),
    sample, tube, gate_type, truth = env$.simLowSepCombn, call = env$.simLowSepCombn,
    fill = list(n = 0L))
  gates <- tibble::tibble(sample = 1L, marker = c("TNF", "IFNg"), gate = 4, gate_cyt = c(3, 4),
    locThresholdMethod = "cap")
  coex <- tibble::tibble(gate_type = rep(c("coex_screened", "coex_unscreened"), each = 2), sample = 1L,
    a = c("TNF", "IFNg", "TNF", "IFNg"), b = c("IFNg", "TNF", "IFNg", "TNF"),
    gateB = 4, cut = c(4, 4, 3, 4), condCut = 4, lowered = c(FALSE, FALSE, TRUE, FALSE))
  cells <- tibble::tibble(sample = 1L, tnf = c(1, 5), ifng = c(1, 5), truth = c("TNF-IFNg-", "TNF+IFNg+"))
  row <- dplyr::bind_cols(scen, tibble::tibble(n_sample_run = n_sample, gates = list(gates),
    coex_gates = list(coex), counts = list(counts), cells = list(cells)))
  validate <- env$.simLowSepValidator(n_sample, "cap")
  expect_identical(validate(row), character())
  expect_match(env$.simLowSepValidator(2L, "cap")(row), "unexpected result sizes")
  expect_match(env$.simLowSepValidator(n_sample, "region")(row), "locThresholdMethod")
  short <- row
  short$coex_gates <- list(coex[1:2, ])
  expect_match(validate(short), "sim_id: 1")

  out <- env$.simLowSepCollate(row)
  expect_named(out, env$.simLowSepCollatedNames)
  expect_equal(nrow(out$counts), nrow(counts))
  expect_true(all(c("scenario_id", "iter", "family", "mismatch_type") %in% names(out$gates)))
  expect_named(out$cells, c("scenario_id", "sim_id", "iter", "sample", "tnf", "ifng", "truth"))
  expect_equal(nrow(out$scenarios), 1L)

  gating <- env$.simLowSepGatingSummary(out$gates, out$coex_gates, out$counts)
  ifng <- gating[gating$marker == "IFNg", ]
  expect_equal(ifng$share_lowered[ifng$gate_type == "coex_unscreened"], 1)
  expect_equal(ifng$share_lowered[ifng$gate_type == "coex_screened"], 0)
  expect_equal(ifng$median_shift_lowered[ifng$gate_type == "coex_unscreened"], 1)
  tnf <- gating[gating$marker == "TNF" & gating$gate_type == "antimode", ]
  expect_equal(tnf$share_lowered, 1)
  expect_equal(tnf$median_shift_lowered, 1)
})

test_that("count validation rejects disagreement with the package statistics", {
  env <- .lowSepEnv()
  counts <- tibble::tibble(
    sample = 1L, tube = rep(c("stim", "uns"), each = 4), gate_type = "antimode",
    truth = "TNF-IFNg-", call = rep(env$.simLowSepCombn, 2), n = c(90L, 4L, 2L, 4L, 98L, 1L, 1L, 0L)
  )
  stats_tbl <- tibble::tibble(
    ind = "2",
    cytCombn = c("F1~-~F2~-~", "F1~+~F2~-~", "F1~-~F2~+~", "F1~+~F2~+~"),
    countStim = c(90L, 4L, 2L, 4L), countUns = c(98L, 1L, 1L, 0L)
  )
  expect_true(env$.simLowSepValidateCounts(counts, stats_tbl))
  stats_tbl$countStim[[4]] <- 5L
  expect_error(env$.simLowSepValidateCounts(counts, stats_tbl), "differ")
  expect_error(env$.simLowSepValidateCounts(counts, stats_tbl[-1, ]), "cover")
})

test_that("the mismatch sections add the mismatch axis", {
  env <- .lowSepEnv()
  tbl <- dplyr::mutate(env$.simLowSepBatchGrid(55700L), value = 1)
  ms <- env$.simLowSepMismatchSection(tbl, "mean_shift")
  expect_identical(levels(ms$mismatch_x), c("None", "0.025", "0.05", "0.2"))
  sd <- env$.simLowSepMismatchSection(tbl, "sd_inflation")
  expect_identical(levels(sd$mismatch_x), c("None", "SD +10%"))
  expect_equal(nrow(sd), 2L * nrow(tbl) / 5L)
  expect_equal(nrow(env$.simLowSepMismatchSection(tbl[tbl$stim_mean_shift == 0.2, ], "sd_inflation")), 0L)
  expect_error(env$.simLowSepMismatchSection(tbl, "other"), "type")
})

test_that("one small dataset gates, validates against the package and plots", {
  testthat::skip_if_not_installed("simcyto")
  env <- .lowSepEnv()
  withr::local_preserve_seed()
  row <- env$.simLowSepRowGrid(env$.simLowSepGrid(55700L)[1, ])[1, ]
  row$n_cell <- 2000
  row$prob_response <- 0.05
  settings <- env$.simLowSepMainSettings(n_sample = 2L)
  run <- function() env$.analysis_with_seed(row$sim_seed, env$.simLowSepScenario(row, settings))
  res <- run()
  expect_equal(nrow(res), 1L)
  expect_identical(env$.simLowSepValidator(2L, "region")(dplyr::bind_cols(row, res)), character())
  out <- env$.simLowSepCollate(dplyr::bind_cols(row, res))
  # Both ordered pairs and both co-expression rules for every stimulated
  # sample, never above the ordinary gate.
  expect_equal(nrow(out$coex_gates), 8L)
  expect_true(all(out$coex_gates$cut <= out$coex_gates$gateB))
  expect_equal(nrow(out$gates), 4L)
  expect_true(all(out$gates$gate_cyt <= out$gates$gate))
  # Every cell of every tube is counted once per rule.
  tube_n <- out$counts |>
    dplyr::group_by(.data$sample, .data$tube, .data$gate_type) |>
    dplyr::summarise(n = sum(.data$n), .groups = "drop")
  expect_true(all(tube_n$n == 2000L))
  expect_setequal(unique(tube_n$gate_type), names(env$.simLowSepGateTypes))
  expect_equal(nrow(out$cells), 2000L)
  # Rerunning the same dataset reproduces it exactly.
  expect_identical(run(), res)

  summ <- env$.simLowSepSummaries(out, settings, mcse = FALSE)
  cells <- dplyr::inner_join(out$cells, out$scenarios, by = "scenario_id")
  plot_dir <- withr::local_tempdir()
  withr::local_dir(plot_dir)
  plots <- list(
    env$.simLowSepPlotHex(cells, out$gates, out$coex_gates),
    env$.simLowSepPlotConditionalDensity(cells, out$gates, out$coex_gates),
    env$.simLowSepPlotGateShift(out$gates, out$coex_gates),
    env$.simLowSepPlotPercentiles(summ$f1_pct, "F1"),
    env$.simLowSepPlotPercentiles(summ$freq_pct, "Frequency", truth = TRUE, scale = 100),
    env$.simLowSepPlotRecovery(summ$recovery),
    env$.simLowSepPlotFrequencyError(summ$freq),
    env$.simLowSepPlotRates(summ$rates),
    env$.simLowSepPlotAdded(summ$gating)
  )
  # The 11b figures, with the dataset relabelled as a mismatch row.
  mm <- function(tbl) dplyr::mutate(tbl, mismatch_type = "mean_shift", stim_mean_shift = 0,
    mismatch_x = factor("None"))
  plots <- c(plots, list(
    env$.simLowSepPlotMismatchRates(dplyr::filter(mm(summ$rates), .data$test == "fisher"), "mean_shift"),
    env$.simLowSepPlotMismatchError(mm(summ$error_pct), "mean_shift"),
    env$.simLowSepPlotMismatchGating(dplyr::filter(mm(summ$gating), .data$gate_type != "ordinary"), "sd_inflation")
  ))
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
    expect_null(p$labels$title)
  }
  expect_length(list.files(plot_dir, recursive = TRUE), 0L)
})

test_that("QMDs 11a and 11b follow the analysis execution contract", {
  for (name in c("11a-sim-low-separation-cyt-pos.qmd", "11b-sim-low-separation-batch.qmd")) {
    qmd <- readLines(file.path(root_dir, "analysis", name), warn = FALSE)
    text <- paste(qmd, collapse = "\n")
    expect_match(text, "run_simulations: true", info = name)
    expect_match(text, "run_plots: false", info = name)
    expect_match(text, "sim-low-separation.R", fixed = TRUE, info = name)
    expect_false(grepl("stimgate:::", text, fixed = TRUE), info = name)
    expect_false(grepl("functionsForBenchmarking-Cyt", text, fixed = TRUE), info = name)
    expect_match(text, "if (isTRUE(run_simulations))", fixed = TRUE, info = name)
    expect_match(text, "if (isTRUE(run_plots))", fixed = TRUE, info = name)
    expect_match(text, "#| label: rerun-one-simulation\n#| include: false\n#| eval: false", fixed = TRUE, info = name)
    expect_match(text, ".analysis_sim_size()", fixed = TRUE, info = name)
    expect_match(text, ".analysis_run_context(", fixed = TRUE, info = name)
    expect_match(text, ".analysis_results_context(", fixed = TRUE, info = name)
    expect_match(text, ".simBandwidthFinishChunk(", fixed = TRUE, info = name)
    expect_match(text, "sim_grid_spec = analysis_grid_spec", fixed = TRUE, info = name)
    expect_match(text, "scenario_fn = .simLowSepScenario", fixed = TRUE, info = name)
    expect_false(grepl("ggtitle|labs\\(title", text), info = name)
    expect_false(grepl("error: true", text, fixed = TRUE), info = name)
    # The final empty chunk is kept for "run all above".
    expect_identical(utils::tail(qmd[nzchar(qmd)], 2L), c("```{r}", "```"), info = name)
  }
  text_b <- paste(readLines(file.path(root_dir, "analysis", "11b-sim-low-separation-batch.qmd"), warn = FALSE),
    collapse = "\n")
  expect_match(text_b, ".simLowSepBatchGrid(simulation_seed)", fixed = TRUE)
  expect_match(text_b, '"mean_shift"', fixed = TRUE)
  expect_match(text_b, "sd_inflation", fixed = TRUE)
  expect_false(file.exists(file.path(root_dir, "analysis", "11-sim-low-separation-cyt-pos.qmd")))
})
