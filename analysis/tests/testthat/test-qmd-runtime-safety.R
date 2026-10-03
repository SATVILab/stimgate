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

  expect_true(grepl("simulation_seed:\\s*12345", content))
  expect_true(grepl(
    'comparison_semantics_version <- "batch-mismatch-comparison-v4"',
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
    ".simComparePromoteIfReady(",
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

test_that("analysis 2a chunks are balanced and plot chunks are guarded", {
  lines <- readLines(
    file.path(root_dir, "analysis", "2a-sim-bw-freq_bs-global.qmd"),
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

test_that("analysis 2a uses shared seeded runners and canonical reads", {
  qmd_path <- file.path(
    root_dir,
    "analysis",
    "2a-sim-bw-freq_bs-global.qmd"
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
  expect_lt(pos("sim_grid_full <- sim_grid"), pos("if (analysis_quick)"))
  expect_lt(pos("if (analysis_quick)"), pos("if (analysis_dev) {"))
  expect_lt(
    pos("if (analysis_dev) {"),
    pos("dplyr::slice_sample(sim_grid_all, prop = 1)")
  )
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


test_that("analysis 3 uses the shared runner and matching canonical results", {
  content <- paste(readLines(file.path(
    root_dir, "analysis", "3-sim-bw-est-base.qmd"
  ), warn = FALSE), collapse = "\n")
  has <- function(x) grepl(x, content, fixed = TRUE)
  expect_true(has(
    'analysis_semantics_version <- "bandwidth-est-base-v5"'
  ))
  for (contract in c(
    "sim_grid_full <- sim_grid",
    "sim_seed = as.integer(simulation_seed + sim_id - 1L)",
    "sim_grid_spec = analysis_grid_spec",
    "scenario_settings = scenario_settings",
    "sim-bandwidth-analysis-run.R",
    ".simBandwidthRunGrid(", ".simBandwidthFinishChunk(",
    ".simBandwidthEstBaseScenario", ".simBandwidthEstBaseCollate",
    ".simBandwidthEstBaseValidate", ".simBandwidthRunRow(",
    "sim_grid_full$sim_id == sim_id_target",
    "retry_errors = sim_retry_errors", "SIM_RETRY_ERRORS",
    "if (is.na(sim_grid_shuffle_seed))",
    "bw_fallback <- NA_real_", "bw_min <- -Inf", "bw_max <- Inf",
    "cap_stim_range <- FALSE", "capStimRange = cap_stim_range",
    ".analysis_results_context(", ".analysis_current_file(",
    "required_params = analysis_required_params",
    "if (isTRUE(run_plots))", "prop_bw_finite",
    "Skipping plots during a multi-chunk simulation render."
  )) {
    expect_true(has(contract), info = contract)
  }
  for (obsolete in c(
    "promote_analysis3_if_ready", "read_base_bw_outputs",
    "validate_base_bw_outputs", "bw-estimate-loop", "set.seed(",
    "bw_list_raw_mtd_chunk_", 'projr::projr_path_get(\n        "output"'
  )) {
    expect_false(has(obsolete), info = obsolete)
  }
  expect_false(grepl("#\\| error:\\s*true", content))
  expect_true(grepl("execute:\\s*warning: false\\s*message: false", content))
  expect_equal(sum(grepl("#| eval: false", strsplit(
    content, "\n", fixed = TRUE
  )[[1]], fixed = TRUE)), 1L)
})


test_that("analysis 4 uses shared seeded runners and canonical reads", {
  qmd_path <- file.path(root_dir, "analysis", "4-sim-bw-est-norm.qmd")
  lines <- readLines(qmd_path, warn = FALSE)
  content <- paste(lines, collapse = "\n")
  has <- function(x) grepl(x, content, fixed = TRUE)
  pos <- function(x) regexpr(x, content, fixed = TRUE)[[1]]

  expect_true(grepl("simulation_seed:\\s*12345", content))
  expect_true(grepl("sim_retry_errors:\\s*true", content))
  expect_true(grepl("warning:\\s*false", content))
  expect_true(has('analysis_semantics_version <- "bandwidth-est-norm-v4"'))
  expect_true(has(
    "sim_seed = as.integer(simulation_seed + dplyr::cur_group_id() - 1L)"
  ))
  # IDs and paired estimator seeds are fixed before the dev filter/shuffle.
  expect_lt(pos("sim_seed = as.integer("), pos("sim_grid_full <- sim_grid"))
  expect_lt(pos("sim_grid_full <- sim_grid"), pos("if (analysis_dev) {"))
  expect_lt(
    pos("if (analysis_dev) {"),
    pos("dplyr::slice_sample(sim_grid_all, prop = 1)")
  )
  expect_true(has('bw_mtd %in% c("hpi1", "hpi1Norm")'))
  expect_false(has('"none", "gaussian"'))
  expect_false(has('"high", "gaussian"'))

  expect_true(has('"sim-bandwidth-analysis-run.R"'))
  expect_true(has(".simBandwidthRunGrid("))
  expect_true(has("scenario_fn = .simBandwidthEstNormScenario"))
  expect_true(has("retry_errors = sim_retry_errors"))
  expect_true(has(".simBandwidthFinishChunk("))
  expect_true(has(".simBandwidthEstNormCollate("))
  expect_true(has(".simBandwidthEstNormValidator("))
  expect_true(has(".analysis_current_file("))
  expect_true(has("required_params = analysis_required_params"))
  expect_true(has("sim_grid_spec = analysis_grid_spec"))
  expect_true(has("scenario_settings = scenario_settings"))
  expect_true(has("capStimRange = FALSE"))
  expect_true(has("Requested normalisation"))
  expect_true(has("pmin(bw_stim, bw_uns)"))
  expect_true(has(
    "Skipping plots during a multi-chunk simulation render."
  ))
  expect_true(has(
    "run_plots is false, so stopping after simulation/collation."
  ))
  expect_true(has("Results were not promoted"))

  expect_false(has("set.seed(as.integer(sim_seed))"))
  expect_false(has("promote_analysis4_if_ready"))
  expect_false(has("read_norm_bw_outputs"))
  expect_false(has("validate_norm_bw_outputs"))
  expect_false(has("collate_suffix"))
  expect_false(has("dir.create"))
  expect_false(has('"cache"'))
  expect_false(has("NA_character_"))
  expect_false(grepl("#\\| error:\\s*true", content))
})

test_that("analysis 4 has one rerun chunk and guarded plotting", {
  lines <- readLines(
    file.path(root_dir, "analysis", "4-sim-bw-est-norm.qmd"),
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

  uses_plots <- vapply(chunks, function(x) {
    code <- x[!grepl("^#\\|", x)]
    any(grepl("ggsave|ggplot\\(", code))
  }, logical(1))
  expect_gt(sum(uses_plots), 0L)
  for (x in chunks[uses_plots]) {
    code <- x[!grepl("^#\\|", x) & nzchar(trimws(x))]
    expect_identical(code[[1]], "if (isTRUE(run_plots)) {")
  }
})


test_that("analysis 5 uses shared seeded runners and canonical reads", {
  qmd_path <- file.path(
    root_dir,
    "analysis",
    "5-sim-bw-est-adaptive.qmd"
  )
  lines <- readLines(qmd_path, warn = FALSE)
  content <- paste(lines, collapse = "\n")
  has <- function(x) grepl(x, content, fixed = TRUE)
  pos <- function(x) regexpr(x, content, fixed = TRUE)[[1]]

  expect_true(grepl("simulation_seed:\\s*12345", content))
  expect_true(grepl("sim_retry_errors:\\s*true", content))
  expect_true(grepl("warning:\\s*false", content))
  expect_true(grepl("message:\\s*false", content))
  expect_true(has('analysis_semantics_version <- "adaptive-bw-est-v3"'))
  expect_true(has("norm_adaptive_ncell <- 2500L"))
  expect_true(has("normAdaptiveNcell = norm_adaptive_ncell"))
  expect_false(has("bw_ncell_upper"))
  expect_false(has("bwNcellMax ="))
  expect_true(has("n_cell_stim_vec <- c(1e3, 5e3, 2e4, 1e5)"))

  # Scenario seeds belong to the data scenario and are fixed on the full grid
  # before the dev filter, shuffling and chunking.
  expect_true(has(
    "sim_seed = as.integer(simulation_seed + .data$data_scenario_id - 1L)"
  ))
  expect_lt(pos("sim_seed = as.integer("), pos("sim_grid_full <- sim_grid"))
  expect_lt(pos("sim_grid_full <- sim_grid"), pos("if (analysis_dev) {"))
  expect_lt(
    pos("if (analysis_dev) {"),
    pos("dplyr::slice_sample(sim_grid_all, prop = 1)")
  )
  expect_false(has("1e4"))
  expect_false(has("sample_n("))

  expect_true(has('"sim-bandwidth-analysis-run.R"'))
  expect_true(has(".simBandwidthRunGrid("))
  expect_true(has("scenario_fn = .simBandwidthEstAdaptiveScenario"))
  expect_true(has("retry_errors = sim_retry_errors"))
  expect_true(has(".simBandwidthFinishChunk("))
  expect_true(has(".simBandwidthEstAdaptiveValidate("))
  expect_true(has(".simBandwidthEstAdaptiveCollate("))
  expect_true(has("run_ctx <- .analysis_results_context("))
  expect_true(has(".analysis_current_file("))
  expect_true(has("required_params = analysis_required_params"))
  expect_true(has("sim_grid_spec = analysis_grid_spec"))
  expect_true(has("scenario_settings = scenario_settings"))
  expect_true(has("Skipping plots during a simulation render."))
  expect_true(has(
    "run_plots is false, so stopping after simulation/collation."
  ))
  expect_true(has("means are conditional on finite"))

  expect_false(has("set.seed(as.integer(sim_seed))"))
  expect_false(has("promote_analysis5_if_ready"))
  expect_false(has("knitr::knit_exit()"))
  expect_false(has("old_plan <- future::plan()"))
  expect_false(has("bw_list_raw_chunk"))
  expect_false(has("collate_suffix"))
  expect_false(has("bw_mtd_norm"))
  expect_false(has("bw_tbl_long"))
  expect_false(has("saveRDS("))
  expect_false(grepl("#\\| error:\\s*true", content))
  expect_identical(
    lengths(regmatches(content, gregexpr("projr_path_get(", content,
      fixed = TRUE
    ))),
    1L
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
  expect_true(any(grepl("sim_grid_full", rerun, fixed = TRUE)))
  expect_true(any(grepl(".simBandwidthRunRow(", rerun, fixed = TRUE)))

  # Every chunk that plots or saves figures is guarded by run_plots.
  plots <- vapply(chunks, function(x) {
    any(grepl("ggsave|ggplot\\(", x[!grepl("^#\\|", x)]))
  }, logical(1))
  expect_gt(sum(plots), 0L)
  for (x in chunks[plots]) {
    code <- x[!grepl("^#\\|", x) & nzchar(trimws(x))]
    expect_identical(code[[1]], "if (isTRUE(run_plots)) {")
  }
})

test_that("analysis 6 presentation chunks are guarded and rerun is singular", {
  lines <- readLines(file.path(
    root_dir, "analysis", "6-sim-bw-freq_bs-adaptive.qmd"
  ), warn = FALSE)
  expect_equal(length(grep("^```\\{r", lines)),
               length(grep("^```\\s*$", lines)))
  chunks <- .qmd_r_chunks(lines)
  disabled <- vapply(chunks, function(x) {
    any(grepl("^#\\|\\s*eval:\\s*false", x))
  }, logical(1))
  expect_equal(sum(disabled), 1L)
  expect_true(any(grepl("label: rerun-one-simulation",
                       chunks[disabled][[1]], fixed = TRUE)))
  for (chunk in chunks[!disabled]) {
    code <- chunk[!grepl("^#\\|", chunk) & nzchar(trimws(chunk))]
    expect_no_error(parse(text = code))
    if (any(grepl("ggsave|bw_tbl_results_|knitr::kable", code))) {
      expect_identical(code[[1]], "if (isTRUE(run_plots)) {")
    }
  }
  expect_false(any(grepl("projr::projr_path_get", lines[-seq_len(40L)],
                        fixed = TRUE)))
})
