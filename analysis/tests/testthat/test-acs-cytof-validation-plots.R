root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)
script_plot <- file.path(root_dir, "scripts", "r", "acs_cytof-plot_cyt.R")

env <- new.env(parent = getNamespace("stimgate"))
source(script_plot, local = env)

.acs_validation_fixture <- function() {
  tidyr::expand_grid(
    method = c("stimgate", "fbeta", "tailgate"),
    pop = c("CD4 T cells", "B cells"),
    cyt = c("IFNg", "IL2"),
    stim = c("p1", "mtb", "ebv", "p4"),
    SampleID = paste0("sample", 1:4)
  ) |>
    dplyr::group_by(method, pop, cyt, stim) |>
    dplyr::mutate(
      freq_bs_man = seq(0.1, 0.4, length.out = dplyr::n()),
      freq_bs_auto = .data$freq_bs_man * 1.1,
      freq_stim_man = .data$freq_bs_man + 0.02,
      freq_uns_man = 0.02
    ) |>
    dplyr::ungroup()
}

test_that("CCC is one for identical vectors", {
  expect_equal(env$.acsCytofValidationCcc(1:5, 1:5), 1)
  expect_true(is.na(env$.acsCytofValidationCcc(rep(1, 5), 1:5)))
})

test_that("correlation table contains PCC and CCC by method and stratum", {
  comparison_tbl <- .acs_validation_fixture()
  correlation_tbl <- env$.acsCytofValidationCorrelationTable(comparison_tbl)

  expect_true(all(c("method", "pop", "cyt", "stim", "pcc", "ccc") %in%
    names(correlation_tbl)))
  expect_equal(nrow(correlation_tbl), 3 * 2 * 2 * 4)
  expect_true(all(abs(correlation_tbl$pcc - 1) < 1e-10))
})

test_that("validation plotting helpers return ggplot objects", {
  comparison_tbl <- .acs_validation_fixture()
  correlation_tbl <- env$.acsCytofValidationCorrelationTable(comparison_tbl)

  expect_s3_class(
    env$.acsCytofValidationPlotScatter(comparison_tbl, "stimgate"),
    "ggplot"
  )
  expect_s3_class(
    env$.acsCytofValidationPlotCorrelation(
      correlation_tbl,
      method = "stimgate",
      metric = "ccc",
      realPopulationsOnly = TRUE
    ),
    "ggplot"
  )
})

test_that("analysis 10 correlation and heatmap chunks run with comparison fixtures", {
  qmd_path <- file.path(
    root_dir, "analysis", "10-real-compare-acs-cytof-validation.qmd"
  )
  lines <- readLines(qmd_path, warn = FALSE)
  chunk_code <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    expect_length(start, 1L)
    end <- which(lines == "```" & seq_along(lines) > start)[1L]
    expect_false(is.na(end))
    parse(text = lines[seq.int(start + 1L, end - 1L)])
  }

  setup <- paste(
    deparse(chunk_code("setup"), width.cutoff = 500L), collapse = "\n"
  )
  expect_true(grepl(
    'source(file.path(scripts_r_dir, "acs_cytof-plot_cyt.R"))',
    setup,
    fixed = TRUE
  ))

  chunk_env <- new.env(parent = getNamespace("stimgate"))
  source(script_plot, local = chunk_env)
  chunk_env$manual_comparison_tbl <- .acs_validation_fixture()
  eval(chunk_code("correlation-table"), envir = chunk_env)
  expect_equal(
    chunk_env$correlation_tbl,
    env$.acsCytofValidationCorrelationTable(chunk_env$manual_comparison_tbl)
  )

  # Capture the document's printed plots without opening graphics devices.
  plots <- list()
  chunk_env$print <- function(x, ...) {
    plots[[length(plots) + 1L]] <<- x
    invisible(x)
  }
  heatmap_labels <- c("show-real-pop-heatmaps", "show-all-pop-heatmaps")
  chunk_env$run_plots <- TRUE
  for (label in heatmap_labels) {
    eval(chunk_code(label), envir = chunk_env)
  }
  expect_length(plots, 12L)
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))

  plots <- list()
  chunk_env$run_plots <- FALSE
  for (label in heatmap_labels) {
    eval(chunk_code(label), envir = chunk_env)
  }
  expect_length(plots, 0L)
})
