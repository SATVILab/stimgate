test_that("bandwidth and bias renders explain how to create missing cached results", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  documents <- c(
    "2a-sim-bw-freq_bs-global.qmd", "2b-sim-bias_uns-freq_bs.qmd"
  )
  prefixes <- c("bw-manual", "bias-uns")
  for (i in seq_along(documents)) {
    lines <- readLines(file.path(root, "analysis", documents[[i]]))
    chunk <- function(label) {
      start <- which(lines == paste0("#| label: ", label))
      end <- which(lines == "```" & seq_along(lines) > start)[1L]
      parse(text = lines[seq.int(start + 1L, end - 1L)])
    }
    env <- new.env(parent = baseenv())
    env$interactive <- function() FALSE
    env$run_simulations <- FALSE
    env$run_plots <- TRUE
    env$analysis_key <- "test"
    env$root_dir <- root
    env$analysis_required_params <- list()
    env$.analysis_results_context <- function(...) stop("No canonical result.")
    command <- paste0(
      "RUN_SIMULATIONS=true RUN_PLOTS=false quarto render analysis/", documents[[i]]
    )
    expect_error(
      eval(chunk(paste0(prefixes[[i]], "-parallel")), env),
      command, fixed = TRUE
    )
    expect_false(exists("run_ctx", env, inherits = FALSE))

    env$run_ctx <- list(read_only = TRUE)
    env$.analysis_current_file <- function(...) stop("Required cache file missing.")
    expect_error(
      eval(chunk(paste0(prefixes[[i]], "-collate")), env),
      command, fixed = TRUE
    )
  }
})

test_that("bandwidth and bias set-up resolves Quarto's document directory", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  for (document in c(
    "2a-sim-bw-freq_bs-global.qmd", "2b-sim-bias_uns-freq_bs.qmd"
  )) {
    lines <- readLines(file.path(root, "analysis", document))
    start <- which(lines == 'if (!file.exists(file.path(root_dir, "DESCRIPTION"))) {')
    end <- which(lines == 'scripts_r_dir <- file.path(root_dir, "scripts", "r")')
    env <- new.env(parent = baseenv())
    for (directory in c(root, file.path(root, "analysis"))) {
      env$root_dir <- directory
      eval(parse(text = lines[seq.int(start, end)]), env)
      expect_identical(env$root_dir, root)
      expect_true(file.exists(file.path(env$scripts_r_dir, "analysis-runtime.R")))
    }
  }
})
