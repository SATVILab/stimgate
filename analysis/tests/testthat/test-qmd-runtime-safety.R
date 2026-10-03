root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

test_that("analysis QMDs do not overwrite sourced helper functions", {
  script_paths <- list.files(
    file.path(root_dir, "scripts", "r"),
    pattern = "[.]R$",
    full.names = TRUE
  )
  helper_names <- unique(unlist(lapply(script_paths, function(path) {
    lines <- readLines(path, warn = FALSE)
    matches <- regmatches(
      lines,
      regexec(
        "^\\s*([.A-Za-z][A-Za-z0-9._]*)\\s*<-\\s*function\\s*\\(",
        lines,
        perl = TRUE
      )
    )
    vapply(
      matches[lengths(matches) > 0L],
      function(match) match[[2]],
      character(1)
    )
  })))

  violations <- character()
  qmd_paths <- list.files(
    file.path(root_dir, "analysis"),
    pattern = "[.]qmd$",
    full.names = TRUE
  )

  for (qmd_path in qmd_paths) {
    lines <- readLines(qmd_path, warn = FALSE)
    for (i in seq_along(lines)) {
      match <- regmatches(
        lines[[i]],
        regexec(
          "^\\s*([.A-Za-z][A-Za-z0-9._]*)\\s*<-\\s*(.*)$",
          lines[[i]],
          perl = TRUE
        )
      )[[1]]
      if (length(match) == 0L) {
        next
      }

      lhs <- match[[2]]
      rhs <- trimws(match[[3]])
      if (
        lhs %in% helper_names &&
          !grepl("^function\\s*\\(", rhs, perl = TRUE)
      ) {
        violations <- c(
          violations,
          paste0(basename(qmd_path), ":", i, ": ", trimws(lines[[i]]))
        )
      }
    }
  }

  expect_identical(violations, character())
})

test_that("analysis 8 is paired, transactional, and read-only for plots", {
  qmd_path <- file.path(
    root_dir,
    "analysis",
    "8-sim-compare-freq_bs-batch.qmd"
  )
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("simulation_seed:\\s*1", content))
  expect_true(grepl(
    'comparison_semantics_version <- "batch-mismatch-comparison-v2"',
    content,
    fixed = TRUE
  ))
  expect_true(grepl("analysis_dev <- isTRUE(.isDev())", content, fixed = TRUE))
  expect_true(grepl("base_scenario_id = dplyr::row_number()", content, fixed = TRUE))
  expect_true(grepl(
    "sim_seed = as.integer(simulation_seed + base_scenario_id - 1L)",
    content,
    fixed = TRUE
  ))
  expect_false(grepl(".simCompareRunScenarioUnseeded", content, fixed = TRUE))

  expect_true(grepl(
    "run_ctx <- .analysis_results_context(",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("paired_mismatch_rng = TRUE", content, fixed = TRUE))
  expect_true(grepl("analysis_grid_spec = analysis_grid_spec", content, fixed = TRUE))
  expect_true(grepl("retryErrors = TRUE", content, fixed = TRUE))
  expect_true(grepl(".simCompareGridOutputStatus(", content, fixed = TRUE))
  expect_true(grepl(
    "F-beta comparator preflight failed",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "cytoUtils' is required for the tailgate comparison",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Refusing to promote analysis 8",
    content,
    fixed = TRUE
  ))

  expect_true(grepl(
    "run_plots is false, so stopping after simulation/collation.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Skipping plots during a multi-chunk simulation render.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("knitr::knit_exit()", content, fixed = TRUE))

  expect_true(grepl(".analysis_current_file(", content, fixed = TRUE))
  expect_true(grepl("simulation_seed = simulation_seed", content, fixed = TRUE))
  expect_true(grepl("analysis_dev = analysis_dev", content, fixed = TRUE))
  expect_true(grepl("n_sample_sim = n_sample_sim", content, fixed = TRUE))
  expect_true(grepl("n_iter_sim = n_iter_sim", content, fixed = TRUE))

  expect_true(grepl(
    "dir.create(dirname(path_p), recursive = TRUE, showWarnings = FALSE)",
    content,
    fixed = TRUE
  ))
})


.qmd_r_chunks <- function(lines) {
  starts <- grep("^```\\{r", lines)
  fences <- grep("^```\\s*$", lines)
  lapply(starts, function(i) {
    end <- fences[fences > i][[1]]
    lines[seq.int(i + 1L, length.out = end - i - 1L)]
  })
}

test_that("analysis 2 chunks are balanced and plot chunks are guarded", {
  lines <- readLines(
    file.path(root_dir, "analysis", "2-sim-bw-freq_bs-global.qmd"),
    warn = FALSE
  )
  expect_identical(
    length(grep("^```\\{r", lines)),
    length(grep("^```\\s*$", lines))
  )
  chunks <- .qmd_r_chunks(lines)

  is_eval_false <- vapply(chunks, function(x) {
    any(grepl("^#\\|\\s*eval:\\s*false", x))
  }, logical(1))
  expect_identical(sum(is_eval_false), 1L)
  rerun <- chunks[is_eval_false][[1]]
  expect_true(any(grepl("label: rerun-one-simulation", rerun, fixed = TRUE)))
  expect_true(any(grepl("sim_id_target <-", rerun, fixed = TRUE)))
  expect_true(any(grepl("sim_grid_full", rerun, fixed = TRUE)))
  expect_true(any(grepl(".simBandwidthRunRow(", rerun, fixed = TRUE)))
  expect_true(any(grepl(".simBandwidthFindSimOutput(", rerun, fixed = TRUE)))
  expect_false(any(grepl("^#\\|\\s*#", unlist(chunks))))

  uses_results <- vapply(chunks, function(x) {
    code <- x[!grepl("^#\\|", x)]
    any(grepl("ggsave|bw_tbl_results_|knitr::kable", code)) &&
      !any(grepl("read_current", code, fixed = TRUE))
  }, logical(1))
  expect_gt(sum(uses_results), 0L)
  for (x in chunks[uses_results]) {
    code <- x[!grepl("^#\\|", x) & nzchar(trimws(x))]
    expect_identical(code[[1]], "if (isTRUE(run_plots)) {")
  }
})

test_that("analysis 2 uses shared seeded runners and canonical reads", {
  qmd_path <- file.path(
    root_dir,
    "analysis",
    "2-sim-bw-freq_bs-global.qmd"
  )
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
  has <- function(x) grepl(x, content, fixed = TRUE)
  pos <- function(x) regexpr(x, content, fixed = TRUE)[[1]]

  expect_true(grepl("simulation_seed:\\s*12345", content))
  expect_true(grepl("sim_retry_errors:\\s*true", content))
  expect_true(has('analysis_semantics_version <- "global-bw-freq-v4"'))
  expect_true(has("sim_seed = as.integer(simulation_seed + sim_id - 1L)"))
  # sim_id/sim_seed are fixed on the full grid before filtering/shuffling.
  expect_lt(pos("sim_seed = as.integer("), pos("sim_grid_full <- sim_grid"))
  expect_lt(pos("sim_grid_full <- sim_grid"), pos("if (analysis_dev)"))
  expect_lt(pos("if (analysis_dev)"), pos("dplyr::slice_sample(prop = 1)"))
  expect_true(has("n_cell == max(n_cell)"))
  expect_false(has("n_cell == 1e5"))

  expect_true(has('"sim-bandwidth-analysis-run.R"'))
  expect_true(has(".simBandwidthRunGrid("))
  expect_true(has("scenario_fn = .simBandwidthFreqBsGlobalScenario"))
  expect_true(has("retry_errors = sim_retry_errors"))
  expect_true(has(".simBandwidthFinishChunk("))
  expect_true(has(".simBandwidthFreqBsGlobalCollate("))
  expect_true(has(".analysis_current_file("))
  expect_true(has("required_params = analysis_required_params"))
  expect_true(has("sim_grid_spec = analysis_grid_spec"))
  expect_true(has("scenario_settings = scenario_settings"))
  expect_true(has(
    "Skipping plots during a multi-chunk simulation render."
  ))
  expect_true(has(
    "run_plots is false, so stopping after simulation/collation."
  ))

  expect_false(has("set.seed(as.integer(sim_seed))"))
  expect_false(has("promote_analysis2_if_ready"))
  expect_false(has("knitr::knit_exit()"))
  expect_false(has("bw_list_raw_chunk"))
  expect_false(has("collate_suffix"))
  expect_false(has("tibble::tibble(\n                transformation"))
  expect_false(grepl("dens_tbl|rug_tbl|mean_line_tbl|main_settings", content))
  expect_false(has("current_manifest$params"))
  expect_identical(
    lengths(regmatches(content, gregexpr("projr_path_get(", content,
      fixed = TRUE
    ))),
    1L
  )
  expect_false(has("make_bw_colour_values <- function"))
  expect_false(has("make_bw_linetype_scale <- function"))
  expect_true(has("format_bw_lab(.data$bw)"))
  expect_true(has("abs(propRespEst - propRespTruth) / propRespTruth"))
})


test_that("analysis 3 is chunk-stable, read-only, and retains estimator failure coverage", {
  qmd_path <- file.path(
    root_dir,
    "analysis",
    "3-sim-bw-est-base.qmd"
  )
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("simulation_seed:\\s*12345", content))
  expect_true(grepl(
    'analysis_semantics_version <- "bandwidth-est-base-v3"',
    content,
    fixed = TRUE
  ))
  expect_true(grepl("analysis_grid_spec", content, fixed = TRUE))
  expect_true(grepl("sim_grid_spec = analysis_grid_spec", content, fixed = TRUE))
  expect_true(grepl(
    "sim_seed = as.integer(simulation_seed + sim_id - 1L)",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "set.seed(as.integer(sim_seed))",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("bw_fallback <- NA_real_", content, fixed = TRUE))
  expect_false(grepl("0.23482348792138919129198282389", content, fixed = TRUE))
  expect_true(grepl("cap_stim_range <- FALSE", content, fixed = TRUE))
  expect_equal(
    lengths(regmatches(
      content,
      gregexpr("capStimRange = cap_stim_range", content, fixed = TRUE)
    )),
    2L
  )

  expect_true(grepl(
    "run_ctx <- .analysis_results_context(",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("run_ctx$chunk_dir", content, fixed = TRUE))
  expect_true(grepl(".analysis_current_file(", content, fixed = TRUE))
  expect_true(grepl("analysis_required_params", content, fixed = TRUE))

  expect_true(grepl("expected_rows_per_sim", content, fixed = TRUE))
  expect_true(grepl("row_count_bad_ids", content, fixed = TRUE))
  expect_true(grepl("seed_ok", content, fixed = TRUE))
  expect_true(grepl("expected_full_ids", content, fixed = TRUE))
  expect_true(grepl(
    "Refusing to promote analysis 3",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "No simulations were assigned to this chunk; marked it complete.",
    content,
    fixed = TRUE
  ))

  expect_true(grepl("n_bw_total", content, fixed = TRUE))
  expect_true(grepl("n_bw_stim_finite", content, fixed = TRUE))
  expect_true(grepl("n_bw_uns_finite", content, fixed = TRUE))
  expect_true(grepl("n_bw_finite", content, fixed = TRUE))
  expect_true(grepl("prop_bw_finite", content, fixed = TRUE))
  expect_true(grepl(
    "pmin(.data$bw_stim, .data$bw_uns)",
    content,
    fixed = TRUE
  ))

  expect_false(grepl("#\\| error:\\s*true", content))
  expect_true(grepl("old_plan <- future::plan()", content, fixed = TRUE))
  expect_true(grepl(
    "finally = future::plan(old_plan)",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "run_plots is false, so stopping after simulation/collation.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Skipping plots during a multi-chunk simulation render.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("promotion_done <-", content, fixed = TRUE))
  expect_true(grepl(
    "Refusing to plot analysis 3: this run was not promoted",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "dir.create(dirname(path_plot)",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("grid::unit(0.9", content, fixed = TRUE))
})


test_that("analysis 4 uses paired estimator seeds and transactional chunk promotion", {
  qmd_path <- file.path(
    root_dir,
    "analysis",
    "4-sim-bw-est-norm.qmd"
  )
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("simulation_seed:\\s*12345", content))
  expect_true(grepl(
    "sim_seed = as.integer(simulation_seed + dplyr::cur_group_id() - 1L)",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    'bw_mtd %in% c("hpi1", "hpi1Norm")',
    content,
    fixed = TRUE
  ))
  expect_false(grepl(
    'grepl("hpi1", bw_mtd),',
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "set.seed(as.integer(sim_seed))",
    content,
    fixed = TRUE
  ))

  expect_true(grepl(
    "run_ctx <- .analysis_results_context(",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "output_dir = run_ctx$chunk_dir",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("row_count_bad_ids", content, fixed = TRUE))
  expect_true(grepl("seed_ok", content, fixed = TRUE))
  expect_true(grepl("expected_sim_ids", content, fixed = TRUE))
  expect_true(grepl("validation$error_ids", content, fixed = TRUE))
  expect_true(grepl("promote_analysis4_if_ready", content, fixed = TRUE))
  expect_true(grepl("nrow(sim_grid) == 0L", content, fixed = TRUE))
  expect_true(grepl(
    "No simulations were assigned to this chunk; marked it complete.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Refusing to promote analysis 4",
    content,
    fixed = TRUE
  ))

  expect_true(grepl(
    "run_plots is false, so stopping after simulation/collation.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Skipping plots during a multi-chunk simulation render.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("finally = future::plan(old_plan)", content, fixed = TRUE))

  expect_true(grepl("n_total = dplyr::n()", content, fixed = TRUE))
  expect_true(grepl("estimate_rate = .data$n_est / .data$n_total", content, fixed = TRUE))
  expect_true(grepl("analysis_grid_spec", content, fixed = TRUE))
  expect_true(grepl(
    "sim_grid_spec = analysis_grid_spec",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("bw_fallback <- NA_real_", content, fixed = TRUE))
  expect_true(grepl("capStimRange = FALSE", content, fixed = TRUE))
  expect_false(grepl("#\\| error:\\s*true", content))
  expect_true(grepl("analysis4-collation", content, fixed = TRUE))
  expect_true(grepl("Requested normalisation", content, fixed = TRUE))
  expect_true(grepl("n_norm_fallback", content, fixed = TRUE))
  expect_true(grepl("norm_fallback_rate", content, fixed = TRUE))
})


test_that("analysis 5 matches adaptive estimator semantics and is chunk-stable", {
  qmd_path <- file.path(
    root_dir,
    "analysis",
    "5-sim-bw-est-adaptive.qmd"
  )
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("simulation_seed:\\s*12345", content))
  expect_true(grepl(
    'analysis_semantics_version <- "adaptive-bw-est-v2"',
    content,
    fixed = TRUE
  ))
  expect_true(grepl("analysis_grid_spec", content, fixed = TRUE))
  expect_true(grepl("sim_grid_spec = analysis_grid_spec", content, fixed = TRUE))
  expect_true(grepl("bw_fallback <- NA_real_", content, fixed = TRUE))

  expect_true(grepl("norm_adaptive_ncell <- 2500L", content, fixed = TRUE))
  expect_true(grepl(
    "normAdaptiveNcell = norm_adaptive_ncell",
    content,
    fixed = TRUE
  ))
  expect_false(grepl("bw_ncell_upper", content, fixed = TRUE))
  expect_false(grepl("bwNcellMax =", content, fixed = TRUE))

  expect_true(grepl("data_scenario_id", content, fixed = TRUE))
  expect_true(grepl(
    "sim_seed = as.integer(simulation_seed + .data$data_scenario_id - 1L)",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "set.seed(as.integer(sim_seed))",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    ".simBandwidthEnsureCurrentCheckout(root_dir)",
    content,
    fixed = TRUE
  ))

  expect_true(grepl(
    "run_ctx <- .analysis_results_context(",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(".analysis_current_file(", content, fixed = TRUE))
  expect_true(grepl("analysis_required_params", content, fixed = TRUE))
  expect_true(grepl("run_ctx$chunk_dir", content, fixed = TRUE))

  expect_true(grepl("expected_rows_per_sim", content, fixed = TRUE))
  expect_true(grepl("row_count_bad_ids", content, fixed = TRUE))
  expect_true(grepl("seed_ok", content, fixed = TRUE))
  expect_true(grepl("expected_full_ids", content, fixed = TRUE))
  expect_true(grepl("promote_analysis5_if_ready", content, fixed = TRUE))
  expect_true(grepl("nrow(sim_grid) == 0L", content, fixed = TRUE))
  expect_true(grepl(
    "Refusing to promote analysis 5",
    content,
    fixed = TRUE
  ))
  expect_false(grepl("#\\| error:\\s*true", content))

  expect_true(grepl("n_total = dplyr::n()", content, fixed = TRUE))
  expect_true(grepl("prop_est = n_est / n_total", content, fixed = TRUE))
  expect_true(grepl(
    "means are conditional on finite estimates",
    content,
    fixed = TRUE
  ))

  expect_true(grepl("old_plan <- future::plan()", content, fixed = TRUE))
  expect_true(grepl(
    "finally = future::plan(old_plan)",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "run_plots is false, so stopping after simulation/collation.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Render again with run_simulations = FALSE and run_plots = TRUE.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(".write_rds_atomic(", content, fixed = TRUE))
  expect_false(grepl("saveRDS(", content, fixed = TRUE))
})
