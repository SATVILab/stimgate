test_that("simulation QMDs share actionable cache errors without creating run state", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  documents <- c(
    "2a-sim-bw-freq_bs-global.qmd", "2b-sim-bias_uns-freq_bs.qmd",
    "3-sim-bw-est-base.qmd", "4-sim-bw-est-norm.qmd",
    "_archive/5-sim-bw-est-adaptive.qmd", "_archive/6-sim-bw-freq_bs-adaptive.qmd",
    "7-sim-compare-freq_bs.qmd", "8-sim-compare-freq_bs-batch.qmd"
  )
  for (i in seq_along(documents)) {
    lines <- readLines(file.path(root, "analysis", documents[[i]]))
    # Find the simulation chunk by its shared runtime call, not its display label.
    starts <- grep("^```\\{r", lines)
    ends <- which(lines == "```")
    chunks <- lapply(starts, function(start) {
      lines[seq.int(start + 1L, ends[ends > start][[1L]] - 1L)]
    })
    code <- Filter(function(x) {
      any(grepl("run_ctx <- .analysis_run_context(", x, fixed = TRUE))
    }, chunks)[[1L]]
    env <- new.env(parent = baseenv())
    source(file.path(root, "scripts", "r", "analysis-runtime.R"), local = env)
    env$interactive <- function() FALSE
    env$run_simulations <- FALSE
    env$run_plots <- TRUE
    env$analysis_key <- c("sim", "test")
    env$analysis_qmd <- file.path("analysis", documents[[i]])
    env$root_dir <- root
    missing_dir <- tempfile("missing-analysis-cache-")
    env$.analysis_cache_dir <- function(..., create) {
      expect_false(create)
      missing_dir
    }
    branch <- Filter(function(expr) {
      is.call(expr) && identical(expr[[1L]], as.name("if")) &&
        grepl("run_ctx <- .analysis_run_context(", paste(deparse(expr), collapse = "\n"), fixed = TRUE)
    }, as.list(parse(text = code)))[[1L]]
    expect_error(eval(branch, env), paste0(
      "RUN_SIMULATIONS=true RUN_PLOTS=false SIM_SIZE=final quarto render analysis/", documents[[i]]
    ), fixed = TRUE)
    expect_false(dir.exists(missing_dir))
    expect_false(exists("run_ctx", env, inherits = FALSE))
  }
})

test_that("all QMDs resolve the checkout from its root and analysis directory", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  documents <- list.files(file.path(root, "analysis"), pattern = "[.]qmd$", full.names = TRUE)
  for (document in documents) {
    lines <- readLines(document)
    start <- grep('^root_dir <- normalizePath', lines)
    expect_length(start, 1L)
    expect_identical(lines[[start + 2L]], "knitr::opts_knit$set(root.dir = root_dir)")
    env <- new.env(parent = baseenv())
    for (directory in c(root, file.path(root, "analysis"))) {
      withr::with_dir(directory, eval(parse(text = lines[start:(start + 1L)]), env))
      expect_identical(env$root_dir, root)
    }
  }
})
