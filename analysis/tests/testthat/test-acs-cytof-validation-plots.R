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

  # Deparse each setup expression on its own so line wrapping cannot split it.
  setup <- vapply(
    chunk_code("setup"),
    function(expr) paste(deparse(expr, width.cutoff = 500L), collapse = " "),
    character(1)
  )
  expect_true(any(grepl(
    'source(file.path(scripts_r_dir, "acs_cytof-plot_cyt.R"))',
    setup,
    fixed = TRUE
  )))

  chunk_env <- new.env(parent = getNamespace("stimgate"))
  source(script_style, local = chunk_env)
  source(script_plot, local = chunk_env)
  chunk_env$manual_comparison_tbl <- .acs_validation_fixture()
  chunk_env$validation_methods <- chunk_env$.acsCytofValidationMethods()
  chunk_env$run_plots <- FALSE
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
  # Each figure shows one method, so there is no method-set duplication.
  for (label in plot_labels) {
    out <- capture.output(eval(chunk_code(label), envir = chunk_env))
    expect_true(any(grepl("#### Method: StimGate", out, fixed = TRUE)))
    expect_false(any(grepl("Without Tailgate", out, fixed = TRUE)))
  }
  expect_length(plots, 18L)
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  expect_equal(
    vapply(
      plots[1:6],
      function(plot) as.character(unique(plot$data$method)),
      character(1)
    ),
    rep(chunk_env$validation_methods, each = 2L)
  )
  expect_true(all(vapply(
    plots[1:6], function(plot) dplyr::n_distinct(plot$data$pop) == 1L, logical(1)
  )))

  plots <- list()
  chunk_env$run_plots <- FALSE
  for (label in plot_labels) {
    out <- capture.output(eval(chunk_code(label), envir = chunk_env))
    expect_false(any(grepl("^#", out)))
  }
  expect_length(plots, 0L)
})

test_that("analysis 10 keeps default-settings F-beta for its appendix only", {
  lines <- readLines(qmd_path, warn = FALSE)
  chunk_code <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    expect_length(start, 1L)
    end <- which(lines == "```" & seq_along(lines) > start)[1L]
    parse(text = lines[seq.int(start + 1L, end - 1L)])
  }
  appendix_methods <- env$.acsCytofValidationAppendixMethods()
  expect_identical(appendix_methods, "fbeta_default")

  main <- .acs_validation_fixture()
  all_tbl <- dplyr::bind_rows(
    main,
    dplyr::mutate(dplyr::filter(main, .data$method == "fbeta"), method = "fbeta_default")
  )
  attr(all_tbl, "manifest") <- list(x = 1)

  # The main results drop the default-settings rows but keep the attributes.
  main_tbl <- env$.acsCytofMainComparisonTable(all_tbl)
  expect_setequal(unique(main_tbl$method), env$.acsCytofValidationMethods())
  expect_identical(attr(main_tbl, "manifest"), list(x = 1))

  chunk_env <- new.env(parent = getNamespace("stimgate"))
  source(script_runtime, local = chunk_env)
  source(script_style, local = chunk_env)
  source(script_plot, local = chunk_env)
  chunk_env$manual_comparison_all_tbl <- all_tbl
  chunk_env$appendix_methods <- appendix_methods
  chunk_env$run_plots <- FALSE
  eval(chunk_code("appendix-correlation-table"), envir = chunk_env)
  expect_setequal(
    unique(chunk_env$appendix_correlation_tbl$method), appendix_methods
  )

  plots <- list()
  chunk_env$print <- function(x, ...) {
    plots[[length(plots) + 1L]] <<- x
    invisible(x)
  }
  chunk_env$run_plots <- TRUE
  for (label in c("appendix-t-cell-correlations", "appendix-all-correlations")) {
    out <- capture.output(eval(chunk_code(label), envir = chunk_env))
    expect_true(any(grepl("#### Method: F-beta", out, fixed = TRUE)))
    expect_true(any(grepl("Concordance correlation", out, fixed = TRUE)))
  }
  # One method, two metrics, two population sets.
  expect_length(plots, 4L)
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
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
    "requiredMethods = c(validation_methods, appendix_methods)",
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

test_that("validation scatter figures are saved per method and population", {
  comparison_tbl <- .acs_validation_fixture()
  parent_dir <- tempfile("acs-validation-sets-")
  dir.create(parent_dir)
  withr::defer(unlink(parent_dir, recursive = TRUE))
  target_dir <- file.path(parent_dir, "validation-figures")

  env$.acsCytofValidationSavePlots(comparison_tbl, target_dir)

  for (population in c("CD4_T_cells", "B_cells")) {
    expect_setequal(
      list.files(file.path(target_dir, "scatter-plots", population)),
      paste0(c("stimgate", "fbeta", "tailgate"), ".pdf")
    )
  }
  expect_length(
    list.files(file.path(target_dir, "scatter-plots"), recursive = TRUE), 6L
  )
  expect_length(list.files(file.path(target_dir, "heatmaps")), 3L * 2L * 2L)
  expect_false(dir.exists(file.path(target_dir, "no_tailgate")))
})

test_that("validation scatter uses stimulus shapes and independent cytokine axes", {
  rows <- .acs_validation_fixture() |>
    dplyr::filter(.data$pop == "CD4 T cells")
  p <- env$.acsCytofValidationPlotScatter(rows, "stimgate")
  built <- ggplot2::ggplot_build(p)
  points <- built$data[[4]]
  expect_setequal(points$shape, c(17, 15, 3, 16))
  expect_true(all(points$alpha == 0.5))
  expect_true(all(points$size == 0.9))
  expect_equal(length(unique(built$layout$layout$SCALE_X)), 2L)
  expect_equal(length(unique(built$layout$layout$SCALE_Y)), 2L)
  expect_equal(p$facet$params$ncol, 3)
  expect_match(p$labels$x, "%", fixed = TRUE)
  expect_match(p$labels$y, "%", fixed = TRUE)
  expect_identical(
    built$plot$scales$get_scales("colour")$get_labels(),
    built$plot$scales$get_scales("shape")$get_labels()
  )
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("heatmaps distinguish excluded strata from unavailable eligible correlations", {
  rows <- .acs_validation_fixture() |>
    dplyr::filter(.data$method == "stimgate", .data$pop == "CD4 T cells")
  # Low manual signal excludes mtb without changing the manuscript rules.
  rows$freq_stim_man[rows$stim == "mtb"] <- 0.01
  # Eligible p1 strata have no finite automated pairs.
  rows$freq_bs_auto[rows$stim == "p1"] <- NA_real_
  tbl <- env$.acsCytofValidationCorrelationTable(rows)
  expect_false("mtb" %in% tbl$stim)
  expect_true(all(is.na(tbl$pcc[tbl$stim == "p1"])))
  expect_setequal(attr(tbl, "excludedStrata")$stim, "mtb")
  p <- env$.acsCytofValidationPlotCorrelation(tbl, "stimgate")
  expect_true(all(p$data$cellLabel[p$data$excluded] == ""))
  expect_true(all(p$data$cellLabel[as.character(p$data$stim) == "Secreted Mtb proteins"] == "NA"))
  expect_true(all(p$data$textColour[p$data$excluded] == "black"))
  expect_true(all(p$data$textColour[as.character(p$data$stim) == "EBV and CMV"] == "white"))
  built <- ggplot2::ggplot_build(p)
  expect_true(all(built$data[[2]]$fill == "white"))
  expect_true(all(built$data[[3]]$size == 3))
  expect_null(p$labels$caption)
  expect_no_error(ggplot2::ggplotGrob(p))

  # Dark endpoint fills need white labels; the pale middle needs black.
  tbl$pcc <- rep(c(-1, 0), length.out = nrow(tbl))
  p <- env$.acsCytofValidationPlotCorrelation(tbl, "stimgate")
  expect_true(all(p$data$textColour[!p$data$excluded & p$data$pcc == -1] == "white"))
  expect_true(all(p$data$textColour[!p$data$excluded & p$data$pcc == 0] == "black"))
})

test_that("thesis heatmaps stack methods with shared scales and stimulus columns", {
  rows <- .acs_validation_fixture()
  expect_identical(env$.acsCytofValidationMethodLabels()[["fbeta"]], "F-beta (zeros removed)")
  expect_identical(env$.acsCytofValidationMethodLabels()[["fbeta_default"]], "F-beta")
  expect_error(
    env$.acsCytofThesisPlotCorrelation(env$.acsCytofValidationCorrelationTable(rows)),
    "F-beta with all cells"
  )
  rows <- dplyr::bind_rows(
    rows,
    dplyr::mutate(dplyr::filter(rows, .data$method == "fbeta"),
      method = "fbeta_default", freq_bs_auto = .data$freq_bs_auto * 2)
  )
  rows$freq_stim_man[rows$stim == "mtb"] <- 0.01
  tbl <- env$.acsCytofValidationCorrelationTable(rows)
  for (metric in c("pcc", "ccc")) {
    p <- env$.acsCytofThesisPlotCorrelation(tbl, metric)
    built <- ggplot2::ggplot_build(p)
    expect_identical(levels(p$data$method), c("StimGate", "F-beta", "Tailgate"))
    expect_equal(
      p$data[[metric]][as.character(p$data$method) == "F-beta" & !p$data$excluded],
      tbl[[metric]][tbl$method == "fbeta_default" & tbl$stim != "p4" &
        tbl$pop %in% env$.acsCytofValidationRealPopulations()]
    )
    expect_setequal(unique(built$layout$layout$ROW), 1:3)
    expect_setequal(unique(built$layout$layout$COL), 1:3)
    expect_equal(built$plot$scales$get_scales("fill")$limits, c(-1, 1))
    expect_true(all(built$data[[2]]$fill == "white"))
    expect_true(all(p$data$cellLabel[p$data$excluded] == ""))
    expect_null(p$labels$caption)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})


test_that("thesis palette variants preserve scores and separate method rows", {
  rows <- .acs_validation_fixture()
  rows <- dplyr::bind_rows(
    rows,
    dplyr::mutate(dplyr::filter(rows, .data$method == "fbeta"), method = "fbeta_default")
  )
  rows$freq_stim_man[rows$stim == "mtb"] <- 0.01
  tbl <- env$.acsCytofValidationCorrelationTable(rows)
  tbl$ccc <- rep(c(-0.5, 0, 0.5, 1), length.out = nrow(tbl))
  baseline <- env$.acsCytofThesisPlotCorrelation(tbl, "ccc")
  for (palette in c("standard", "midpoint_0_5", "f1", "viridis")) {
    p <- env$.acsCytofThesisPlotCorrelation(tbl, "ccc", palette = palette)
    built <- ggplot2::ggplot_build(p)
    fill <- built$plot$scales$get_scales("fill")
    expect_identical(p$data$ccc, baseline$data$ccc)
    expect_identical(p$data$cellLabel, baseline$data$cellLabel)
    expect_equal(fill$limits, c(-1, 1))
    separators <- p$layers[[4]]$data
    expect_identical(as.character(separators$method), c("StimGate", "F-beta"))
    expect_true(all(separators$yintercept == 0.5))
    expect_equal(nrow(built$data[[4]]), 6L)
    expect_true(all(built$layout$layout$ROW[
      match(built$data[[4]]$PANEL, built$layout$layout$PANEL)
    ] %in% 1:2))
    expect_true(all(is.na(p$data$tileColour[p$data$excluded])))
    expect_true(all(p$data$cellLabel[p$data$excluded] == ""))
    if (palette == "standard") {
      expect_true(all(is.na(p$data$tileColour)))
    } else {
      expect_true(all(p$data$tileColour[!p$data$excluded] == "grey80"))
      expect_identical(fill$map(c(-1, -0.5, 0)), rep(fill$map(0), 3))
      expect_true(all(p$data$textColour[!p$data$excluded & p$data$ccc <= 0] == "white"))
    }
    if (palette == "midpoint_0_5") {
      expect_identical(toupper(fill$map(0.5)), "#FFFFFF")
      endpoints <- rev(RColorBrewer::brewer.pal(11, "RdBu"))[c(1, 11)]
      expect_equal(grDevices::col2rgb(fill$map(c(0, 1))), grDevices::col2rgb(endpoints))
      expect_true(all(p$data$textColour[!p$data$excluded & p$data$ccc == 0.5] == "black"))
    }
    if (palette %in% c("f1", "viridis")) {
      f1tbl <- tbl
      f1tbl$f1 <- pmax(tbl$ccc, 0)
      f1plot <- env$.acsCytofThesisPlotCorrelation(f1tbl, "f1")
      f1fill <- ggplot2::ggplot_build(f1plot)$plot$scales$get_scales("fill")
      expect_equal(f1fill$limits, c(0, 1))
      probes <- seq(0, 1, length.out = 17)
      reference <- ggplot2::scale_fill_viridis_c(limits = c(0, 1))
      expect_equal(fill$map(probes), f1fill$map(probes))
      expect_equal(fill$map(probes), reference$map(probes))
      expect_equal(env$.acsCytofF1Palette(probes), reference$palette(probes))
    }
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})
