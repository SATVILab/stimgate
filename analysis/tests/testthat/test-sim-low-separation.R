root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

.lowSepEnv <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (f in c("analysis-runtime.R", "analysis-plot-style.R", "sim-low-separation.R")) {
    source(file.path(root_dir, "scripts", "r", f), local = env)
  }
  env
}

.lowSepToyCounts <- function(stim, uns) {
  scen <- tibble::tibble(
    sim_id = 1L, separation = "low", mean_pos = 4.5,
    response_level = "higher", prob_response = 0.05, n_cell = 100
  )
  rows <- function(tube, x) {
    tibble::tibble(
      sample = 1L, tube = tube, gate_type = x$gate_type,
      truth = x$truth, call = x$call, n = as.integer(x$n)
    )
  }
  dplyr::bind_cols(scen, dplyr::bind_rows(rows("stim", stim), rows("uns", uns)))
}

test_that("scenario grid assigns IDs and seeds on the full grid", {
  env <- .lowSepEnv()
  grid <- env$.simLowSepGrid(55700L)
  expect_equal(nrow(grid), 8L)
  expect_equal(grid$sim_id, 1:8)
  expect_setequal(unique(grid$mean_pos), c(4.5, 3.5))
  expect_setequal(unique(grid$prob_response), c(0.01, 0.05))
  expect_setequal(unique(grid$n_cell), c(1e5, 5e3))
  expect_identical(grid, env$.simLowSepGrid(55700L))
  expect_false(identical(grid$sim_seed, env$.simLowSepGrid(1L)$sim_seed))
  expect_equal(anyDuplicated(grid$sim_seed), 0L)
})

test_that("only the Gaussian transformation is accepted", {
  env <- .lowSepEnv()
  settings <- env$.simLowSepMainSettings()
  settings$transformation <- "gamma"
  row <- env$.simLowSepGrid(1L)[1, ]
  expect_error(env$.simLowSepSimulate(row, 1L, settings), "Gaussian")
})

test_that("cytokine-positive calls follow the package rule", {
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
    env$.simLowSepCall(x_tnf, x_ifng, gate, "cytpos"),
    c("TNF+IFNg+", "TNF-IFNg-", "TNF-IFNg-", "TNF+IFNg+", "TNF+IFNg+")
  )
  # The package's positivity helper gives the same calls.
  ex <- data.frame(F1 = x_tnf, F2 = x_ifng)
  gate_tbl <- data.frame(chnl = c("F1", "F2"), gate = c(4, 4), gateCyt = c(2, 2))
  pos <- stimgate:::.getPosIndByChnl(ex, gate_tbl, c("F1", "F2"), "cyt")
  expect_equal(
    env$.simLowSepCombn[1L + pos$F1 + 2L * pos$F2],
    env$.simLowSepCall(x_tnf, x_ifng, gate, "cytpos")
  )
  expect_error(env$.simLowSepCall(1, 1, gate, "other"), "gate_type")
})

test_that("recovering IFNg in TNF+ cells corrects the exclusive single-positive overestimate", {
  env <- .lowSepEnv()
  # Stimulated tube: 10 true TNF+IFNg+ cells; the ordinary gates call 6 of them
  # TNF+IFNg- only, the cytokine-positive gates recover them all.
  stim <- tibble::tribble(
    ~gate_type, ~truth, ~call, ~n,
    "ordinary", "TNF+IFNg+", "TNF+IFNg+", 4,
    "ordinary", "TNF+IFNg+", "TNF+IFNg-", 6,
    "ordinary", "TNF-IFNg-", "TNF-IFNg-", 90,
    "cytpos", "TNF+IFNg+", "TNF+IFNg+", 10,
    "cytpos", "TNF-IFNg-", "TNF-IFNg-", 88,
    "cytpos", "TNF-IFNg-", "TNF+IFNg+", 2
  )
  uns <- tibble::tribble(
    ~gate_type, ~truth, ~call, ~n,
    "ordinary", "TNF-IFNg-", "TNF-IFNg-", 100,
    "cytpos", "TNF-IFNg-", "TNF-IFNg-", 100
  )
  counts <- .lowSepToyCounts(stim, uns)
  freq <- env$.simLowSepFrequencies(counts)
  get <- function(gt, q, col) freq[[col]][freq$gate_type == gt & freq$quantity == q]
  expect_equal(get("ordinary", "TNF+IFNg-", "prop_bs_est"), 0.06)
  expect_equal(get("cytpos", "TNF+IFNg-", "prop_bs_est"), 0)
  expect_equal(get("ordinary", "TNF+IFNg-", "prop_bs_truth"), 0)
  expect_equal(get("ordinary", "TNF+", "prop_bs_est"), get("ordinary", "TNF+", "prop_bs_truth"))
  expect_equal(get("cytpos", "IFNg+", "prop_bs_est"), 0.12)
  expect_equal(get("cytpos", "TNF+IFNg+", "error_bs"), 0.02)
  expect_setequal(unique(freq$quantity_type), c("marginal", "exclusive"))

  rec <- env$.simLowSepRecovery(counts)
  ifng <- rec[rec$marker == "IFNg", ]
  expect_equal(ifng$sensitivity[ifng$gate_type == "ordinary"], 0.4)
  expect_equal(ifng$sensitivity[ifng$gate_type == "cytpos"], 1)
  expect_equal(ifng$n_contaminating[ifng$gate_type == "cytpos"], 2L)
  expect_equal(ifng$fdp[ifng$gate_type == "cytpos"], 2 / 12)

  summ <- env$.simLowSepFrequencySummary(freq)
  row <- summ[summ$quantity == "TNF+IFNg-", ]
  expect_equal(row$n_closer, 1L)
  expect_equal(row$n_further + row$n_unchanged, 0L)
})

test_that("count validation rejects disagreement with the package statistics", {
  env <- .lowSepEnv()
  counts <- tibble::tibble(
    sample = 1L, tube = rep(c("stim", "uns"), each = 4), gate_type = "cytpos",
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

test_that("cache reads reject results made with other settings", {
  env <- .lowSepEnv()
  path <- withr::local_tempfile(fileext = ".rds")
  grid <- env$.simLowSepGrid(1L)
  settings <- env$.simLowSepCacheSettings(grid, 2L, env$.simLowSepMainSettings(), 1L, list(sim_size = "final"))
  env$.simLowSepWriteCache(list(x = 1), settings, path)
  withr::local_envvar(ANALYSIS_EXPECTED_RUN_ID = NA)
  expect_equal(env$.simLowSepReadCache(path, settings, c("sim", "test")), list(x = 1))
  other <- env$.simLowSepCacheSettings(grid, 3L, env$.simLowSepMainSettings(), 1L, list(sim_size = "final"))
  expect_error(env$.simLowSepReadCache(path, other, c("sim", "test")), "different settings")
})

test_that("one small scenario gates, validates against the package and plots", {
  env <- .lowSepEnv()
  withr::local_preserve_seed()
  row <- env$.simLowSepGrid(55700L)[1, ]
  row$n_cell <- 2000
  row$prob_response <- 0.05
  res <- env$.simLowSepRunScenario(row, n_sample = 2L)
  expect_named(res, c("gates", "counts", "cells"))
  expect_equal(nrow(res$gates), 4L)
  expect_true(all(res$gates$gate_cyt <= res$gates$gate))
  # Every cell of every tube is counted once per gate type.
  tube_n <- res$counts |>
    dplyr::group_by(.data$sample, .data$tube, .data$gate_type) |>
    dplyr::summarise(n = sum(.data$n), .groups = "drop")
  expect_true(all(tube_n$n == 2000L))
  expect_equal(nrow(res$cells), 2000L)

  # Rerunning the same scenario reproduces it exactly.
  expect_identical(env$.simLowSepRunScenario(row, n_sample = 2L), res)

  recovery <- env$.simLowSepRecovery(res$counts)
  freq <- env$.simLowSepFrequencies(res$counts)
  plot_dir <- withr::local_tempdir()
  withr::local_dir(plot_dir)
  plots <- list(
    env$.simLowSepPlotHex(res$cells, res$gates),
    env$.simLowSepPlotConditionalDensity(res$cells, res$gates),
    env$.simLowSepPlotGateShift(res$gates),
    env$.simLowSepPlotRecovery(recovery),
    env$.simLowSepPlotFrequencyError(freq)
  )
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
    expect_null(p$labels$title)
  }
  expect_length(list.files(plot_dir, recursive = TRUE), 0L)
})

test_that("QMD 11 follows the analysis execution contract", {
  qmd <- readLines(file.path(root_dir, "analysis", "11-sim-low-separation-cyt-pos.qmd"), warn = FALSE)
  text <- paste(qmd, collapse = "\n")
  expect_match(text, "run_simulations: true")
  expect_match(text, "run_plots: false")
  expect_match(text, "sim-low-separation.R", fixed = TRUE)
  expect_false(grepl("stimgate:::", text, fixed = TRUE))
  expect_false(grepl("functionsForBenchmarking-Cyt", text, fixed = TRUE))
  expect_match(text, "if (isTRUE(run_simulations))", fixed = TRUE)
  expect_match(text, "if (isTRUE(run_plots))", fixed = TRUE)
  expect_match(text, "#| label: rerun-one-simulation\n#| include: false\n#| eval: false", fixed = TRUE)
  expect_match(text, ".analysis_sim_size()", fixed = TRUE)
  expect_false(grepl("ggtitle|labs\\(title", text))
})
