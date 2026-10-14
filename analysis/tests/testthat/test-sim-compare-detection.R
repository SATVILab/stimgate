.detection_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c(
    "analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-compare-detection.R", "sim-compare-null.R"
  )) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

.detection_fixture <- function() {
  tidyr::expand_grid(
    iter = 1:5, sample = 1:3,
    method = c("tailgate", "fbeta", "stimgate")
  ) |>
    dplyr::mutate(
      sim_seed = 20L, transformation = "gaussian", mean_pos_setting = "high",
      mean_pos = 5, n_cell = 20000, prob_background = 0.0001,
      prob_response = 0, propUnsTruth = 0.0001, propRespTruth = 0,
      sample_perturbation_sd = 0, condition_perturbation_sd = 0,
      cluster_perturbation_sd = 0, n_cell_uns_relative_to_stim = 1,
      nCellStim = 20000, nCellUns = 20000, nPosStim = 3 * .data$sample + .data$iter,
      nPosUns = 2, propRespEst = (.data$nPosStim - .data$nPosUns) / 20000,
      error = NA_character_
    )
}

test_that("log Fisher tails equal fisher.test and remain accurate below underflow", {
  env <- .detection_env()
  tables <- rbind(
    c(8, 1, 20, 20), c(0, 0, 100, 100), c(20, 20, 20, 20),
    c(1, 8, 20, 40), c(120, 2, 1000, 900), c(1000, 0, 1000, 1000)
  )
  d <- as.data.frame(tables)
  names(d) <- c("nPosStim", "nPosUns", "nCellStim", "nCellUns")
  d$propRespEst <- d$nPosStim / d$nCellStim - d$nPosUns / d$nCellUns
  out <- env$.simCompareFisherScores(d)
  for (i in 1:5) {
    tab <- rbind(
      c(d$nPosStim[i], d$nCellStim[i] - d$nPosStim[i]),
      c(d$nPosUns[i], d$nCellUns[i] - d$nPosUns[i])
    )
    expect_equal(out$fisher_p[i], stats::fisher.test(tab, alternative = "greater")$p.value,
      tolerance = 1e-12
    )
  }
  expect_equal(out$fisher_log_p[6], -lchoose(2000, 1000), tolerance = 1e-12)
  expect_equal(out$fisher_p[6], 0)
  expect_true(is.finite(out$fisher_log10[6]))
  expect_equal(out$fisher_log2, out$fisher_log10 * log2(10))
  expect_equal(out$sbb_b01_min[2], 1)
  expect_equal(out$sbb_log10[6], -(1 + out$fisher_log_p[6] + log(-out$fisher_log_p[6])) / log(10))
  expect_equal(out$sbb_b01_min[1], -exp(1) * out$fisher_p[1] * log(out$fisher_p[1]))
})

test_that("invalid gates and missing estimates are missing tests, never negatives", {
  env <- .detection_env()
  d <- .detection_fixture()[1:6, ]
  d$nPosStim <- c(NA, -1, 20001, 1.5, 4, 4)
  d$error[5] <- "estimator failed"
  d$propRespEst[6] <- NA_real_
  scores <- env$.simCompareFisherScores(d)
  expect_false(any(scores$detection_valid))
  expect_true(all(is.na(scores$fisher_log_p)))
  out <- env$.simCompareDetectionSummary(scores, "method", cutoffs = 0.05)
  expect_true(all(is.na(out$detection_rate)))
  expect_true(all(out$n_valid == 0))
  expect_equal(sum(out$n_missing), 6)
})

test_that("dataset bootstrap rates preserve paired methods and pooled denominators", {
  env <- .detection_env()
  d <- .detection_fixture()
  d$propRespEst[d$sample == 1] <- -0.01
  d$nPosUns[d$iter == 1 & d$method == "fbeta"] <- NA_real_
  scores <- env$.simCompareFisherScores(d)
  out <- env$.simCompareDetectionSummary(scores, "method", cutoffs = c(0.01, 0.05, 0.2), mcse = TRUE, retain_draws = TRUE)
  off <- env$.simCompareDetectionSummary(scores, "method", cutoffs = c(0.01, 0.05, 0.2))
  expect_equal(out$detection_rate, off$detection_rate)
  expect_equal(levels(out$method), env$.analysis_method_order)
  expect_equal(
    out$.boot_detection_rate[out$method == "stimgate"],
    out$.boot_detection_rate[out$method == "tailgate"]
  )
  expect_true(all(!out$interval_available_detection_rate[out$method == "fbeta"]))
  row <- out[out$method == "stimgate" & out$alpha == 0.05, ]
  x <- with(scores[scores$method == "stimgate", ], as.numeric(propRespEst > 0 & fisher_log_p <= log(0.05)))
  expected <- env$.analysis_mcse_pooled_cols(x, as.character(1:5)[rep(1:5, each = 3)],
    mean, "detection_rate", "sim_seed:20",
    mcse = TRUE
  )
  expect_equal(row$.boot_detection_rate, expected$.boot_detection_rate)
})

test_that("cell labels, method order, evidence caps and missing heatmap tiles build", {
  env <- .detection_env()
  expect_equal(
    env$.analysis_label_cells(c(500, 1000, 5000, 20000, 100000)),
    c("500", "1k", "5k", "20k", "100k")
  )
  scores <- env$.simCompareFisherScores(.detection_fixture())
  evidence <- env$.simCompareEvidenceSummary(
    scores,
    c("transformation", "mean_pos_setting", "n_cell", "method")
  )
  evidence$fisher_log10_p90 <- 20
  p <- env$.simCompareNullPlot(evidence, "fisher_log10_p90", "Evidence",
    display_cap = 12, reference = env$.simCompareEvidenceReference("fisher_log10")
  )
  build <- ggplot2::ggplot_build(p)
  expect_equal(
    as.character(build$plot$scales$get_scales("colour")$get_breaks()),
    c("stimgate", "fbeta", "tailgate")
  )
  expect_true(all(p$data$value == 12))
  expect_true(all(evidence$fisher_log10_p90 == 20))
  expect_true(all(p$data$capped))
  evidence$n_cell <- rep(c(5000, 100000), length.out = nrow(evidence))
  two <- env$.simCompareNullPlot(evidence, "fisher_log10_p90", "Evidence")
  expect_true(is.factor(two$data$n_cell))
  expect_equal(ggplot2::ggplot_build(two)$layout$panel_params[[1]]$x$get_labels(), c("5k", "100k"))
  rates <- env$.simCompareDetectionSummary(scores,
    c("transformation", "mean_pos_setting", "n_cell", "prob_background", "method"),
    cutoffs = 0.05
  )
  extra <- rates[1, ]
  extra$n_cell <- 100000
  extra$prob_background <- 0.01
  rates <- dplyr::bind_rows(rates, extra)
  heat <- env$.simCompareNullDetectionHeatmap(rates)
  expect_true(any(is.na(heat$data$detection_rate)))
  expect_false(any(as.character(heat$data$method) == "tailgate_default"))
  expect_equal(heat$scales$get_scales("fill")$limits, c(0, 1))
  expect_s3_class(ggplot2::ggplot_build(heat), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(env$.simCompareNullDetectionPlot(rates)), "ggplot_built")
})

test_that("ROC matching rejects different requested or realised backgrounds", {
  env <- .detection_env()
  null <- env$.simCompareFisherScores(.detection_fixture())
  response <- null
  response$prob_response <- 0.0005
  response$background_relative_to_response <- 0.2
  response$propRespTruth <- 0.0005
  roc <- env$.simCompareDetectionRoc(null, response, cutoffs = c(0.001, 0.05, 0.2), mcse = TRUE)
  expect_true(all(roc$audit$matched))
  expect_gt(nrow(roc$curves), 0)
  expect_s3_class(ggplot2::ggplot_build(env$.simCompareDetectionRocPlot(roc$curves)), "ggplot_built")
  response$prob_response <- 0.002
  expect_equal(nrow(env$.simCompareDetectionRoc(null, response)$curves), 0)
  response$prob_response <- 0.0005
  response$propUnsTruth <- 0
  expect_equal(nrow(env$.simCompareDetectionRoc(null, response)$curves), 0)
  response$propUnsTruth <- 0.0001
  # No thresholds span 5%: sensitivity at 5% is unavailable, not interpolated.
  absent <- env$.simCompareDetectionRoc(null, response, cutoffs = 1e-12)$at_5
  expect_false(any(absent$region_bracketed))
  expect_true(all(is.na(absent$sensitivity_at_supported_cutoff)))
})

test_that("ROC compatibility follows Analysis 7 source settings and grid", {
  env <- .detection_env()
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  source(file.path(root, "scripts", "r", "sim-misc.R"), local = env)
  qmd <- file.path(root, "analysis", "7-sim-compare-freq_bs.qmd")
  params <- env$.simCompareDetectionResponseParams(qmd, root, 12345L)
  expect_equal(params$stimgate_bias_uns_factor, 1)
  expect_true(all(params$analysis_grid_spec$n_cell %in% c(1000, 5000, 20000, 100000)))
  expect_true(0.01 %in% params$analysis_grid_spec$prob_response)
  expect_equal(params$n_iter_sim, 20L)
  copy <- tempfile(fileext = ".qmd")
  withr::defer(unlink(copy))
  writeLines(sub("stimgate_bias_uns_factor <- 1", "stimgate_bias_uns_factor <- 2",
    readLines(qmd),
    fixed = TRUE
  ), copy)
  changed <- env$.simCompareDetectionResponseParams(copy, root, 12345L)
  expect_equal(changed$stimgate_bias_uns_factor, 2)
  expect_false(isTRUE(all.equal(params$analysis_grid_spec, changed$analysis_grid_spec)))
})

test_that("quick-mode detection keeps estimates with one dataset and suppresses intervals", {
  env <- .detection_env()
  scores <- env$.simCompareFisherScores(dplyr::filter(.detection_fixture(), .data$iter == 1))
  out <- env$.simCompareDetectionSummary(scores, "method", cutoffs = 0.05, mcse = TRUE)
  expect_true(all(is.finite(out$detection_rate)))
  expect_false(any(out$interval_available_detection_rate))
  expect_true(all(is.na(out$detection_rate_lower)))
})

test_that("linear cell-count views keep actual ticks and two-level views stay discrete", {
  env <- .detection_env()
  cells <- c(1000, 5000, 20000, 100000)
  d <- data.frame(cells = cells, value = seq_along(cells))
  p <- ggplot2::ggplot(d, ggplot2::aes(cells, value)) +
    ggplot2::geom_point() + env$.analysis_scale_cells(cells, log = FALSE)
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$data[[1]]$x, cells)
  expect_equal(p$scales$get_scales("x")$breaks, cells)
  expect_equal(built$layout$panel_params[[1]]$x$get_labels(), c("1k", "5k", "20k", "100k"))
  d$cells <- env$.analysis_cell_x(rep(c(5000, 100000), 2))
  p <- ggplot2::ggplot(d, ggplot2::aes(cells, value)) +
    ggplot2::geom_point() + env$.analysis_scale_cells(d$cells, log = FALSE)
  built <- ggplot2::ggplot_build(p)
  expect_equal(as.numeric(built$data[[1]]$x), rep(1:2, 2))
  expect_equal(built$layout$panel_params[[1]]$x$get_labels(), c("5k", "100k"))
})
