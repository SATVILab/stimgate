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

test_that("analysis 8 uses deterministic scenario seeds and full-grid promotion", {
  qmd_path <- file.path(
    root_dir,
    "analysis",
    "8-sim-compare-freq_bs-batch.qmd"
  )
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("simulation_seed:\\s*1", content))
  expect_true(grepl("comparison_semantics_version", content, fixed = TRUE))
  expect_true(grepl(
    "sim_seed = as.integer(simulation_seed + sim_id - 1L)",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "set.seed(as.integer(row$sim_seed[[1]]))",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "path_progress_file <- run_ctx$progress_file",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("recursive\\s*=\\s*TRUE", content))
  expect_true(grepl("expected_sim_ids", content, fixed = TRUE))
  expect_true(grepl("Refusing to promote analysis 8", content, fixed = TRUE))
  expect_true(grepl(".analysis_current_file", content, fixed = TRUE))
  expect_true(grepl("results_available", content, fixed = TRUE))
  expect_true(grepl(
    "skipping summary and plots for this chunk",
    content,
    fixed = TRUE,
    ignore.case = TRUE
  ))
})


test_that("analysis 2 is chunk-stable, read-only when not simulating, and validates promotion", {
  qmd_path <- file.path(
    root_dir,
    "analysis",
    "2-sim-bw-freq_bs-global.qmd"
  )
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("simulation_seed:\\s*12345", content))
  expect_true(grepl("analysis_semantics_version", content, fixed = TRUE))
  expect_true(grepl("analysis_quick", content, fixed = TRUE))
  expect_true(grepl("analysis_dev", content, fixed = TRUE))
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
  expect_false(grepl("furrr:::make_seeds", content, fixed = TRUE))
  expect_false(grepl(
    "12345 + as.integer(sim_grid_chunk_index)",
    content,
    fixed = TRUE
  ))

  expect_true(grepl(
    "run_ctx <- .analysis_results_context(",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "collate_output_dir <- if (results_read_only)",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("run_ctx$chunk_dir", content, fixed = TRUE))

  expect_true(grepl("expected_chunk_ids", content, fixed = TRUE))
  expect_true(grepl("output_error_ids", content, fixed = TRUE))
  expect_true(grepl("expected_full_ids", content, fixed = TRUE))
  expect_true(grepl("promote_analysis2_if_ready", content, fixed = TRUE))
  expect_true(grepl("nrow(sim_grid) == 0L", content, fixed = TRUE))
  expect_true(grepl(
    "No simulations were assigned to this chunk; marked it complete.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "Refusing to promote analysis 2",
    content,
    fixed = TRUE
  ))

  expect_true(grepl(
    "run_plots is false, so stopping after simulation/collation.",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("knitr::knit_exit()", content, fixed = TRUE))
  expect_true(grepl(
    "Skipping plots during a multi-chunk simulation render.",
    content,
    fixed = TRUE
  ))

  expect_false(grepl(
    "make_bw_colour_values <- function",
    content,
    fixed = TRUE
  ))
  expect_false(grepl(
    "make_bw_linetype_scale <- function",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("format_bw_lab(.data$bw)", content, fixed = TRUE))
  expect_true(grepl(".write_rds_atomic(", content, fixed = TRUE))

  expect_true(grepl('.data$method == "loc_sample"', content, fixed = TRUE))
  expect_true(grepl("is.finite(.data$propRespTruth)", content, fixed = TRUE))
  expect_true(grepl("is.finite(.data$propRespEst)", content, fixed = TRUE))
  expect_true(grepl(
    "abs(propRespEst - propRespTruth) / propRespTruth",
    content,
    fixed = TRUE
  ))
  expect_false(grepl(
    "filter(is.finite(threshold) & is.finite(propBsEst))",
    content,
    fixed = TRUE
  ))
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
  expect_true(grepl("expected_chunk_ids", content, fixed = TRUE))
  expect_true(grepl("expected_sim_ids", content, fixed = TRUE))
  expect_true(grepl("output_error_ids", content, fixed = TRUE))
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
  expect_true(grepl("sim_grid_definition", content, fixed = TRUE))
  expect_true(grepl("Requested normalisation", content, fixed = TRUE))
  expect_true(grepl(
    "fallback is not yet exposed in row-level",
    content,
    fixed = TRUE
  ))
})
