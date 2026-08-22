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
  "2-sim-bw-freq_bs-global.qmd",
  "3-sim-bw-est-base.qmd",
  "4-sim-bw-est-norm.qmd",
  "5-sim-bw-est-adaptive.qmd",
  "6-sim-bw-freq_bs-adaptive.qmd",
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
  on.exit({
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
  }, add = TRUE)

  Sys.setenv(RUN_SIMULATIONS = "false", RUN_PLOTS = "true")
  env$params <- list(run_simulations = TRUE, run_plots = FALSE)
  override_sims <- env$.as_flag(env$.get_qmd_param_env("run_simulations", "RUN_SIMULATIONS", FALSE))
  override_plots <- env$.as_flag(env$.get_qmd_param_env("run_plots", "RUN_PLOTS", TRUE))
  expect_false(override_sims)
  expect_true(override_plots)
})
