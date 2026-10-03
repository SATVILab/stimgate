.load_est_base_run_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "sim-misc.R", "sim-bandwidth.R",
    "sim-bandwidth-analysis-io.R", "sim-bandwidth-analysis-run.R"
  )) {
    source(file.path(root, "scripts", "r", fn), local = env)
  }
  env
}

.est_base_test_row <- function() {
  tibble::tibble(
    transformation = "gaussian", prob_response = 0.002, n_cell = 1000,
    mean_pos_setting = "high", mean_pos = 8, bw_mtd = "hpi1",
    bias_uns_setting = "low", bias_uns = 0.05,
    sim_id = 7L, sim_seed = 12351L
  )
}

test_that("base scenario forwards settings and reruns the stored RNG stream", {
  env <- .load_est_base_run_env()
  row <- .est_base_test_row()
  settings <- list(
    nSample = 2L, nIter = 1L, bwFallback = NA_real_, bwMin = -Inf,
    bwMax = Inf, capStimRange = FALSE, excMin = TRUE, summarise = FALSE
  )
  # A cheap estimator stub checks forwarding without a cytometry simulation.
  env$.simBandwidthEstBwDirect <- function(...) {
    args <- list(...)
    expect_identical(args, c(settings, list(
      biasUns = 0.05, bwMtd = "hpi1", nCellStim = 1000,
      probResponse = 0.002, meanPos = 8, transformation = "gaussian"
    )))
    tibble::tibble(iter = 1L, sample = c("1", "2"), bw = stats::runif(2))
  }
  withr::local_seed(10L)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  before <- .Random.seed
  project <- withr::local_tempdir()
  ctx <- env$.analysis_run_context(
    c("sim", "bw", "est", "base"), run_id = "rerun",
    path_root = project
  )
  stored <- env$.simBandwidthRunGrid(
    row, scenario_fn = env$.simBandwidthEstBaseScenario,
    settings = settings, run_ctx = ctx, error_col = "error", workers = 1L
  )[[1]]
  RNGkind("Mersenne-Twister")
  set.seed(99L)
  rerun <- env$.simBandwidthRunRow(
    row, env$.simBandwidthEstBaseScenario, settings, error_col = "error"
  )
  expect_identical(stored, rerun)
  expect_true(all(is.na(rerun$error)))
  # Direct row execution also restores the caller's stream.
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", before, envir = .GlobalEnv)
  env$.simBandwidthRunRow(row, env$.simBandwidthEstBaseScenario, settings)
  expect_identical(.Random.seed, before)
})

test_that("base collation preserves failed estimates and validates row counts", {
  env <- .load_est_base_run_env()
  row <- .est_base_test_row()
  tbl <- dplyr::bind_cols(row[rep(1L, 3L), ], tibble::tibble(
    bw_stim = c(0.2, NA_real_, Inf), bw_uns = c(0.1, 0.3, Inf),
    bw = c(0.1, 0.3, Inf), error = NA_character_
  ))
  settings <- list(nSample = 3L, nIter = 1L)
  expect_identical(env$.simBandwidthEstBaseValidate(tbl, settings), character())
  expect_match(
    env$.simBandwidthEstBaseValidate(tbl[1:2, ], settings),
    "unexpected sample-row counts for sim_id: 7"
  )
  res <- env$.simBandwidthEstBaseCollate(tbl, names(row))
  expect_named(res, c("bw_list_raw_mtd", "bw_tbl_results"))
  expect_identical(res$bw_list_raw_mtd, tbl)
  summary <- res$bw_tbl_results
  expect_equal(summary$n_bw_total, 3L)
  expect_equal(summary$n_bw_stim_finite, 1L)
  expect_equal(summary$n_bw_uns_finite, 2L)
  expect_equal(summary$n_bw_finite, 1L)
  expect_equal(summary$prop_bw_finite, 1 / 3)
  expect_equal(summary$mean_bw, 0.1)
  tbl$bw_stim <- NA_real_
  tbl$bw_uns <- NA_real_
  all_failed <- env$.simBandwidthEstBaseCollate(tbl, names(row))$bw_tbl_results
  expect_equal(all_failed$prop_bw_finite, 0)
  expect_true(is.na(all_failed$mean_bw))
  expect_identical(env$.simBandwidthEstBaseValidate(tbl, settings), character())
})

test_that("base full-grid IDs survive dev filtering and external chunking", {
  env <- .load_est_base_run_env()
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  lines <- readLines(file.path(root, "analysis", "3-sim-bw-est-base.qmd"))
  chunk <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1]
    lines[(start + 1L):(end - 1L)]
  }
  build <- function(dev, index, chunks) {
    env$analysis_dev <- dev
    env$simulation_seed <- 12345L
    env$sim_grid_shuffle_seed <- 8L
    env$sim_grid_chunk_index <- index
    env$sim_grid_n_chunks <- chunks
    eval(parse(text = chunk("actual-settings")), env)
    invisible(utils::capture.output(
      eval(parse(text = chunk("bw-estimate-grid")), env)
    ))
    list(full = env$sim_grid_full, selected = env$sim_grid_all,
         chunk = env$sim_grid, spec = env$analysis_grid_spec)
  }
  full <- build(FALSE, 1L, 1L)
  one <- build(TRUE, 1L, 2L)
  two <- build(TRUE, 2L, 2L)
  expect_identical(one$full, full$full)
  expect_identical(one$spec, two$spec)
  expect_gt(nrow(one$selected), 0L)
  combined <- dplyr::bind_rows(one$chunk, two$chunk)
  expect_setequal(combined$sim_id, one$selected$sim_id)
  expect_identical(anyDuplicated(combined$sim_id), 0L)
  expect_equal(combined$sim_seed, 12345L + combined$sim_id - 1L)
  expect_equal(sort(unique(full$selected$n_cell)), c(1e3, 5e3, 2e4, 1e5))
  expect_identical(unique(full$selected$bias_uns_setting), "low")
})

test_that("base failures use typed shared error rows and retry on resume", {
  env <- .load_est_base_run_env()
  row <- .est_base_test_row()
  project <- withr::local_tempdir()
  ctx <- env$.analysis_run_context(
    c("sim", "bw", "est", "base"), run_id = "retry",
    path_root = project
  )
  env$.simBandwidthEstBwDirect <- function(...) stop("estimator crashed")
  run <- function() {
    env$.simBandwidthRunRowResumable(
      row, env$.simBandwidthEstBaseScenario, list(), ctx, 1L,
      error_col = "error"
    )
  }
  failed <- run()
  expect_identical(failed$error, "estimator crashed")
  expect_named(failed, c(names(row), "error"))
  env$.simBandwidthEstBwDirect <- function(...) {
    tibble::tibble(sample = "1", ind = 1L, bw = 0.2)
  }
  retried <- run()
  expect_true(is.na(retried$error))
  expect_identical(env$.simBandwidthChunkMarkerCounts(ctx), list(
    completed = 1L, failed = 0L
  ))
  combined <- dplyr::bind_rows(failed, retried)
  expect_type(combined$ind, "integer")
  expect_type(combined$bw, "double")
  expect_type(combined$sample, "character")
  expect_true(is.na(combined$bw[[1]]))
})
