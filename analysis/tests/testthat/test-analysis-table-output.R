.table_output_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  source(file.path(testthat::test_path(), "../../../scripts/r/analysis-runtime.R"), local = env)
  env
}

test_that("table directories are output siblings and read-only resolution creates nothing", {
  env <- .table_output_env()
  # Tests use a temporary projr project, never the checkout's output.
  root <- .local_projr_root()
  parts <- c("analysis-name", "draft", "coverage")
  path <- env$.analysis_table_dir(parts, root, create = FALSE)
  expect_equal(normalizePath(path, winslash = "/", mustWork = FALSE),
    .projr_output_path(root, "table", "analysis-name", "draft", "coverage"))
  expect_false(dir.exists(file.path(root, "_tmp")))
  expect_equal(env$.analysis_table_dir(parts, root), path)
  expect_true(dir.exists(path))
  expect_false(dir.exists(file.path(root, "output")))
})

test_that("CSV companions retain all rows and keep method-specific coverage in prose", {
  env <- .table_output_env()
  root <- .local_projr_root()
  tbl <- data.frame(method = rep(c("stimgate", "fbeta"), each = 12),
    setting = rep(1:12, 2), n_total = 20L,
    n_failed = c(rep(0L, 12), rep(2L, 12)),
    n_fallback = c(rep(1L, 12), rep(0L, 12)),
    bootstrap_coverage_median = c(rep(1, 12), rep(NA_real_, 12)),
    interval_available_median = c(rep(TRUE, 12), rep(FALSE, 12)),
    paired_datasets = "5 / 20", estimate = seq_len(24) / 7)
  tbl$.boot_median <- rep(list(c(1, NA_real_, 3)), nrow(tbl))
  parts <- c("analysis-name", "draft", "mcse_on", "coverage.csv")
  output <- paste(utils::capture.output(result <- env$.analysis_report_table(
    tbl, parts, "Coverage and estimates by setting.", root)), collapse = "\n")
  path <- .projr_output_path(root, "table", "analysis-name", "draft", "mcse_on", "coverage.csv")
  saved <- readr::read_csv(path, show_col_types = FALSE)
  expect_equal(nrow(saved), 24L)
  expect_equal(saved$estimate, tbl$estimate)
  expect_equal(saved$n_failed, tbl$n_failed)
  expect_true(all(is.na(saved$bootstrap_coverage_median[13:24])))
  expect_equal(saved$.boot_median, rep("1;NA;3", 24))
  expect_identical(result, tbl)
  expect_match(output, "output/table/analysis-name/draft/mcse_on/coverage.csv", fixed = TRUE)
  expect_match(output, "Full table (24 rows)", fixed = TRUE)
  expect_match(output, "Coverage for stimgate", fixed = TRUE)
  expect_match(output, "Coverage for fbeta", fixed = TRUE)
  expect_match(output, "n_failed: 2 to 2", fixed = TRUE)
  expect_match(output, "12 unavailable", fixed = TRUE)
  expect_match(output, "paired_datasets: 5 to 5 / 20 to 20", fixed = TRUE)
  expect_false(grepl("\n\\s*\\|", output, perl = TRUE))
  expect_false(grepl("# A tibble", output, fixed = TRUE))
  expect_false(dir.exists(file.path(root, "_tmp", "sim")))
  expect_false(dir.exists(file.path(root, "output")))

  utils::capture.output(env$.analysis_report_table(tbl[0, ],
    c("analysis-name", "empty.csv"), "Empty coverage.", root))
  expect_equal(nrow(readr::read_csv(.projr_output_path(root, "table", "analysis-name", "empty.csv"),
    show_col_types = FALSE)), 0L)
})

test_that("QMD table exports are guarded and large inline table printers stay absent", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  files <- list.files(file.path(root, "analysis"), pattern = "\\.qmd$", full.names = TRUE)
  n_exports <- 0L
  for (file in files) {
    lines <- readLines(file, warn = FALSE)
    for (start in which(startsWith(lines, "```{r"))) {
      end <- which(lines == "```" & seq_along(lines) > start)[1L]
      # Empty chunks (e.g. a QMD's trailing "run all above" chunk) hold no code.
      if (end <= start + 1L) next
      code <- lines[seq.int(start + 1L, end - 1L)]
      parsed <- parse(text = code)
      if (any(grepl("#\\| (include|eval): false", code))) next
      visit <- function(expr, guarded = FALSE) {
        if (!is.call(expr) && !is.expression(expr)) return(invisible(NULL))
        if (is.call(expr)) {
          name <- paste(deparse(expr[[1]]), collapse = "")
          # Subscripts contain empty argument symbols; they are data lookup,
          # rather than report orchestration.
          if (name %in% c("[", "[[")) return(invisible(NULL))
          if (name == "if") {
            condition <- paste(deparse(expr[[2]]), collapse = "")
            enabled <- grepl("isTRUE(run_plots)", condition, fixed = TRUE) &&
              !grepl("!isTRUE(run_plots)", condition, fixed = TRUE)
            visit(expr[[3]], guarded || enabled)
            if (length(expr) == 4L) visit(expr[[4]], guarded)
            return(invisible(NULL))
          }
          if (name %in% c(".analysis_report_table", ".simBandwidthPrintCoverage")) {
            n_exports <<- n_exports + 1L
            expect_true(guarded, info = paste(basename(file), start, name))
            if (name == ".simBandwidthPrintCoverage") {
              expect_true("table_parts" %in% names(as.list(expr)))
            }
          }
          if (name == "knitr::kable" && !grepl("real-compare-omip", basename(file), fixed = TRUE)) {
            # Simulation reports keep one inline result table: three methods by
            # six outcome columns. The OMIP real-data reports (one donor or a
            # few mice) print small inline tables of their single dataset.
            # Analysis 6 prints its five best settings per method (ten rows).
            allowed <- c(
              "8-sim-compare-freq_bs-batch.qmd" = ".simCompareMethodOutcomeCounts",
              "6-sim-tune-comparators.qmd" = "knitr::kable(rank_display",
              # One row per population and stimulation (four rows).
              "17-explore-acs-cytof-coexpression.qmd" = "knitr::kable(just_summary"
            )
            expect_true(basename(file) %in% names(allowed), info = basename(file))
            expect_match(paste(deparse(expr), collapse = ""), allowed[[basename(file)]], fixed = TRUE)
          }
        }
        for (child in as.list(expr)[if (is.call(expr)) -1L else seq_along(expr)]) visit(child, guarded)
      }
      for (expr in parsed) visit(expr)
    }
  }
  expect_gt(n_exports, 40L)
})

test_that("disabled table chunks write nothing while retaining correlation computation", {
  env <- .table_output_env()
  env$run_plots <- FALSE
  env$.analysis_report_table <- function(...) stop("Disabled reporting wrote a table")
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  chunk <- function(file, label) {
    lines <- readLines(file.path(root, "analysis", file), warn = FALSE)
    start <- which(lines == paste0("#| label: ", label))
    end <- which(lines == "```" & seq_along(lines) > start)[1L]
    parse(text = lines[seq.int(start + 1L, end - 1L)])
  }
  for (entry in list(
    c("1-sim-trans.qmd", "univariate-counts-by-mean-pos"),
    c("2c-sim-test.qmd", "test-summary"),
    c("9-real-compare-acs-cytof.qmd", "show-manual-summary"),
    c("9-real-compare-acs-cytof.qmd", "manual-error-coverage"),
    c("9-real-compare-acs-cytof.qmd", "manual-donor-uncertainty"),
    c("10-real-compare-acs-cytof-validation.qmd", "validation-error-coverage"),
    c("10-real-compare-acs-cytof-validation.qmd", "validation-donor-uncertainty")
  )) {
    expect_no_error(eval(chunk(entry[1], entry[2]), env))
  }
  env$manual_comparison_tbl <- data.frame(method = "stimgate")
  env$.acsCytofValidationCorrelationTable <- function(x) data.frame(pcc = 0.9)
  expect_no_error(eval(chunk("10-real-compare-acs-cytof-validation.qmd", "correlation-table"), env))
  expect_equal(env$correlation_tbl, data.frame(pcc = 0.9))
})
