root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
script_runtime <- file.path(root_dir, "scripts", "r", "analysis-runtime.R")
analysis_dir <- file.path(root_dir, "analysis")

.load_runtime_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_runtime, local = env)
  env
}

primary_qmds <- c(
  "1-sim-trans.qmd",
  "2a-sim-bw-freq_bs-global.qmd",
  "2b-sim-bias_uns-freq_bs.qmd",
  "3-sim-bw-est-base.qmd",
  "4-sim-bw-est-norm.qmd",
  "_archive/5-sim-bw-est-adaptive.qmd",
  "_archive/6-sim-bw-freq_bs-adaptive.qmd",
  "7-sim-compare-freq_bs.qmd",
  "8-sim-compare-freq_bs-batch.qmd"
)

test_that("all primary analysis QMDs exist and declare standard YAML params", {
  for (qmd_name in primary_qmds) {
    qmd_path <- file.path(analysis_dir, qmd_name)
    expect_true(file.exists(qmd_path), info = paste0("Missing QMD: ", qmd_name))

    lines <- readLines(qmd_path, warn = FALSE)
    yaml_delims <- which(lines == "---")
    expect_gte(length(yaml_delims), 2L)

    yaml_block <- lines[(yaml_delims[1] + 1L):(yaml_delims[2] - 1L)]
    yaml_text <- paste(yaml_block, collapse = "\n")

    expect_true(
      grepl("run_simulations:\\s*true", yaml_text, ignore.case = TRUE),
      info = paste0("Expected 'run_simulations: true' in YAML of ", qmd_name)
    )
    expect_true(
      grepl("run_plots:\\s*false", yaml_text, ignore.case = TRUE),
      info = paste0("Expected 'run_plots: false' in YAML of ", qmd_name)
    )
  }
})

test_that("all primary analysis QMDs source analysis-runtime.R and initialize execution flags", {
  for (qmd_name in primary_qmds) {
    qmd_path <- file.path(analysis_dir, qmd_name)
    lines <- readLines(qmd_path, warn = FALSE)
    content <- paste(lines, collapse = "\n")

    expect_true(
      grepl("analysis-runtime\\.R", content),
      info = paste0("Expected analysis-runtime.R sourced in ", qmd_name)
    )
    expect_true(
      grepl("run_simulations\\s*<-\\s*\\.as_flag", content),
      info = paste0("Expected run_simulations initialization in ", qmd_name)
    )
    expect_true(
      grepl("run_plots\\s*<-\\s*\\.as_flag", content),
      info = paste0("Expected run_plots initialization in ", qmd_name)
    )
  }
})

test_that("execution flag contract respects render defaults, interactive fallbacks, and env overrides", {
  env <- .load_runtime_env()

  # 1. Normal Quarto render simulation: params provided (simulations = TRUE, plots = FALSE)
  env$params <- list(run_simulations = TRUE, run_plots = FALSE)
  render_sims <- env$.as_flag(env$.get_qmd_param_env("run_simulations", "RUN_SIMULATIONS", FALSE))
  render_plots <- env$.as_flag(env$.get_qmd_param_env("run_plots", "RUN_PLOTS", TRUE))
  expect_true(render_sims)
  expect_false(render_plots)

  # 2. Interactive execution simulation: no params present -> fallbacks apply
  env$params <- NULL
  inter_sims <- env$.as_flag(env$.get_qmd_param_env("run_simulations", "RUN_SIMULATIONS", FALSE))
  inter_plots <- env$.as_flag(env$.get_qmd_param_env("run_plots", "RUN_PLOTS", TRUE))
  expect_false(inter_sims)
  expect_true(inter_plots)

  # 3. Environment overrides
  old_run_sims <- Sys.getenv("RUN_SIMULATIONS", unset = NA_character_)
  old_run_plots <- Sys.getenv("RUN_PLOTS", unset = NA_character_)
  on.exit(
    {
      if (is.na(old_run_sims)) {
        Sys.unsetenv("RUN_SIMULATIONS")
      } else {
        Sys.setenv(RUN_SIMULATIONS = old_run_sims)
      }
      if (is.na(old_run_plots)) {
        Sys.unsetenv("RUN_PLOTS")
      } else {
        Sys.setenv(RUN_PLOTS = old_run_plots)
      }
    },
    add = TRUE
  )

  Sys.setenv(RUN_SIMULATIONS = "false", RUN_PLOTS = "true")
  env$params <- list(run_simulations = TRUE, run_plots = FALSE)
  override_sims <- env$.as_flag(env$.get_qmd_param_env("run_simulations", "RUN_SIMULATIONS", FALSE))
  override_plots <- env$.as_flag(env$.get_qmd_param_env("run_plots", "RUN_PLOTS", TRUE))
  expect_false(override_sims)
  expect_true(override_plots)
})

test_that("QMDs with run contexts can read canonical results without run_ctx", {
  for (qmd_name in primary_qmds) {
    lines <- readLines(file.path(analysis_dir, qmd_name), warn = FALSE)
    text <- paste(lines, collapse = "\n")
    if (!grepl("run_ctx <- .analysis_run_context(", text, fixed = TRUE)) {
      next
    }

    starts <- grep("^```\\{r", lines)
    ends <- grep("^```$", lines)
    chunks <- vapply(starts, function(s) {
      paste(lines[(s + 1L):(min(ends[ends > s]) - 1L)], collapse = "\n")
    }, character(1))
    pattern_sim <- "run_ctx <- .analysis_run_context("
    ind_sim <- which(grepl(pattern_sim, chunks, fixed = TRUE))[[1]]
    ind_use <- which(
      seq_along(chunks) > ind_sim & grepl("run_ctx", chunks, fixed = TRUE)
    )
    expect_gt(length(ind_use), 0L, label = qmd_name)
    first_use <- chunks[[ind_use[[1]]]]

    # analysis_key is shared by the run context and the read-only fallback.
    expect_true(
      grepl("analysis_key = analysis_key", chunks[[ind_sim]], fixed = TRUE),
      info = qmd_name
    )
    inline_read_only_branch <- grepl(
      "if (isTRUE(run_simulations))",
      chunks[[ind_sim]],
      fixed = TRUE
    ) && grepl(
      "run_ctx <- .analysis_results_context(",
      chunks[[ind_sim]],
      fixed = TRUE
    )
    later_read_only_fallback <- grepl('if (!exists("run_ctx"))', first_use, fixed = TRUE) &&
      grepl("run_ctx <- .analysis_results_context(", first_use, fixed = TRUE)

    expect_true(
      inline_read_only_branch || later_read_only_fallback,
      info = paste0("QMD lacks a read-only results path: ", qmd_name)
    )

    # Older QMDs select read-only mode in the later collation chunk. For that
    # layout, writes and promotion must sit behind the read-only guard. QMDs
    # that branch before creating run_ctx already prevent the read-only path
    # from entering active-run orchestration.
    if (
      !inline_read_only_branch &&
        grepl(".analysis_promote_run(run_ctx)", first_use, fixed = TRUE)
    ) {
      guard <- regexpr("if (!results_read_only) {", first_use, fixed = TRUE)
      expect_gt(guard, 0L, label = qmd_name)
      writes <- gregexpr(
        paste0(
          "saveRDS\\(|\\.write_rds_atomic\\(|",
          "\\.analysis_mark_chunk\\(|\\.analysis_promote_run\\("
        ),
        first_use
      )[[1]]
      expect_true(
        all(writes > guard),
        info = paste0("Collation writes outside read-only guard: ", qmd_name)
      )
    }
  }
})
