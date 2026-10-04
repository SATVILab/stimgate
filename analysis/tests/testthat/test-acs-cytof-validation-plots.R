root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)
script_runtime <- file.path(root_dir, "scripts", "r", "analysis-runtime.R")
script_style <- file.path(root_dir, "scripts", "r", "analysis-plot-style.R")
script_plot <- file.path(root_dir, "scripts", "r", "acs_cytof-plot_cyt.R")
qmd_path <- file.path(
  root_dir,
  "analysis",
  "10-real-compare-acs-cytof-validation.qmd"
)

env <- new.env(parent = getNamespace("stimgate"))
source(script_runtime, local = env)
source(script_style, local = env)
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

test_that("CCC matches the manuscript cccrm estimator on shifted data", {
  expect_equal(env$.acsCytofValidationCcc(1:5, 1:5), 1)
  expect_equal(
    env$.acsCytofValidationCcc(1:5, 2:6),
    0.8,
    tolerance = 1e-12
  )
  expect_equal(
    env$.acsCytofValidationCcc(1:5, 2 * (1:5)),
    8 / 19,
    tolerance = 1e-12
  )
  expect_equal(env$.acsCytofValidationCcc(rep(1, 5), 1:5), 0)
  expect_true(is.na(env$.acsCytofValidationCcc(rep(1, 5), rep(1, 5))))
})

test_that("correlation table contains PCC and CCC by method and stratum", {
  comparison_tbl <- .acs_validation_fixture()
  correlation_tbl <- env$.acsCytofValidationCorrelationTable(comparison_tbl)

  expect_true(all(c("method", "pop", "cyt", "stim", "pcc", "ccc") %in%
    names(correlation_tbl)))
  expect_equal(nrow(correlation_tbl), 3 * 2 * 2 * 4)
  expect_true(all(abs(correlation_tbl$pcc - 1) < 1e-10))
})

test_that("correlation table keeps only manuscript-eligible signal strata", {
  comparison_tbl <- tibble::tibble(
    method = "stimgate",
    pop = "CD4 T cells",
    cyt = "IFNg",
    stim = rep(c("p1", "mtb", "ebv"), each = 4L),
    SampleID = rep(paste0("sample", 1:4), 3L),
    freq_bs_man = c(
      0.10, 0.20, 0.30, 0.40,
      0.10, 0.20, 0.30, 0.40,
      0.001, 0.002, 0.003, 0.004
    ),
    freq_bs_auto = c(
      0.11, 0.22, 0.33, 0.44,
      0.11, 0.22, 0.33, 0.44,
      0.0011, 0.0022, 0.0033, 0.0044
    ),
    freq_stim_man = c(
      0.12, 0.22, 0.32, 0.42,
      0.03, 0.04, 0.05, 0.06,
      0.12, 0.22, 0.32, 0.42
    ),
    freq_uns_man = 0.02
  )

  correlation_tbl <- env$.acsCytofValidationCorrelationTable(comparison_tbl)

  expect_equal(correlation_tbl$stim, "p1")
  expect_equal(correlation_tbl$n, 4L)
  expect_equal(correlation_tbl$pcc, 1, tolerance = 1e-12)
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

test_that("analysis 10 correlation and plot chunks run with comparison fixtures", {
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
  source(script_style, local = chunk_env)
  source(script_plot, local = chunk_env)
  chunk_env$manual_comparison_tbl <- .acs_validation_fixture()
  chunk_env$validation_methods <- chunk_env$.acsCytofValidationMethods()
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
  plot_labels <- c(
    "acs-validation-scatter", "acs-validation-t-cell-correlations",
    "acs-validation-all-correlations"
  )
  chunk_env$run_plots <- TRUE
  # Two method sets: all three methods, then StimGate and F-beta only.
  for (label in plot_labels) {
    out <- capture.output(eval(chunk_code(label), envir = chunk_env))
    expect_true(any(grepl("#### All methods", out, fixed = TRUE)))
    expect_true(any(grepl("#### Without Tailgate", out, fixed = TRUE)))
  }
  expect_length(plots, 25L)
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  expect_equal(
    vapply(
      plots[1:3],
      function(plot) as.character(unique(plot$data$method)),
      character(1)
    ),
    chunk_env$validation_methods
  )
  expect_equal(
    vapply(
      plots[4:5],
      function(plot) as.character(unique(plot$data$method)),
      character(1)
    ),
    c("stimgate", "fbeta")
  )

  plots <- list()
  chunk_env$run_plots <- FALSE
  for (label in plot_labels) {
    out <- capture.output(eval(chunk_code(label), envir = chunk_env))
    expect_false(any(grepl("^#", out)))
  }
  expect_length(plots, 0L)
})

test_that("validation input checks schema, methods, and duplicate keys", {
  comparison_tbl <- .acs_validation_fixture()

  expect_no_error(env$.acsCytofValidationValidateComparisonTable(
    comparison_tbl,
    requiredMethods = c("stimgate", "fbeta", "tailgate")
  ))

  expect_error(
    env$.acsCytofValidationValidateComparisonTable(
      dplyr::select(comparison_tbl, -freq_stim_man)
    ),
    "missing required column"
  )

  expect_error(
    env$.acsCytofValidationValidateComparisonTable(
      dplyr::filter(comparison_tbl, .data$method != "tailgate"),
      requiredMethods = c("stimgate", "fbeta", "tailgate")
    ),
    "missing required method"
  )

  expect_error(
    env$.acsCytofValidationValidateComparisonTable(
      dplyr::bind_rows(comparison_tbl, comparison_tbl[1, , drop = FALSE])
    ),
    "duplicate"
  )
})

test_that("validation directory replacement swaps complete output sets", {
  parent_dir <- tempfile("acs-validation-parent-")
  dir.create(parent_dir)
  withr::defer(unlink(parent_dir, recursive = TRUE))

  target_dir <- file.path(parent_dir, "validation-figures")
  staged_dir <- file.path(parent_dir, "staged")
  dir.create(target_dir)
  dir.create(staged_dir)
  writeLines("old", file.path(target_dir, "old.txt"))
  writeLines("new", file.path(staged_dir, "new.txt"))

  expect_no_error(env$.acsCytofValidationReplaceDirectory(
    stagedDir = staged_dir,
    targetDir = target_dir
  ))
  expect_false(file.exists(file.path(target_dir, "old.txt")))
  expect_equal(readLines(file.path(target_dir, "new.txt")), "new")
})

test_that("failed validation rendering preserves the last good output directory", {
  comparison_tbl <- .acs_validation_fixture()
  parent_dir <- tempfile("acs-validation-save-")
  dir.create(parent_dir)
  withr::defer(unlink(parent_dir, recursive = TRUE))

  target_dir <- file.path(parent_dir, "validation-figures")
  dir.create(target_dir)
  marker_path <- file.path(target_dir, "last-good.txt")
  writeLines("keep me", marker_path)

  old_scatter <- env$.acsCytofValidationPlotScatter
  withr::defer(assign(
    ".acsCytofValidationPlotScatter",
    old_scatter,
    envir = env
  ))
  assign(
    ".acsCytofValidationPlotScatter",
    function(...) stop("plot boom"),
    envir = env
  )

  expect_error(
    env$.acsCytofValidationSavePlots(
      comparisonTbl = comparison_tbl,
      pathDirSave = target_dir
    ),
    "plot boom"
  )
  expect_equal(readLines(marker_path), "keep me")
})

test_that("analysis 10 validates its input and does not delete last good figures", {
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl(
    ".acsCytofValidationValidateComparisonTable(",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "requiredMethods = validation_methods",
    content,
    fixed = TRUE
  ))
  expect_false(grepl(
    "unlink(path_dir_save",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "non-secreted-protein (`p4`) stimulation",
    content,
    fixed = TRUE
  ))

  save_body <- paste(
    deparse(body(env$.acsCytofValidationSavePlots)),
    collapse = "\n"
  )
  expect_true(grepl(".write_rds_atomic(", save_body, fixed = TRUE))
  expect_true(grepl(
    ".acsCytofValidationReplaceDirectory(",
    save_body,
    fixed = TRUE
  ))
})

test_that("validation figures are saved for both method sets in one directory", {
  comparison_tbl <- .acs_validation_fixture()
  parent_dir <- tempfile("acs-validation-sets-")
  dir.create(parent_dir)
  withr::defer(unlink(parent_dir, recursive = TRUE))
  target_dir <- file.path(parent_dir, "validation-figures")

  env$.acsCytofValidationSavePlots(comparison_tbl, target_dir)

  scatter_all <- list.files(file.path(target_dir, "all_methods", "scatter-plots"))
  scatter_no_tg <- list.files(file.path(target_dir, "no_tailgate", "scatter-plots"))
  expect_setequal(scatter_all, paste0(c("stimgate", "fbeta", "tailgate"), ".pdf"))
  expect_setequal(scatter_no_tg, paste0(c("stimgate", "fbeta"), ".pdf"))
  expect_length(
    list.files(file.path(target_dir, "no_tailgate", "heatmaps")),
    2L * 2L * 2L
  )
})
