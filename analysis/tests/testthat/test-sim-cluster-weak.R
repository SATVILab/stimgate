.weak_test_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c("analysis-runtime.R", "analysis-plot-style.R", "sim-compare-freq_bs.R", "sim-cluster-weak.R")) {
    source(file.path(root, "scripts", "r", fn), local = env)
  }
  env
}

test_that("weak-response scoring compares the same cells with strict gates", {
  env <- .weak_test_env()
  cells <- tibble::tibble(ind = "2", condition = "Stimulated",
    expression = c(0, 1, 2, 3, 4), label = c("gn", "gn", "gp", "gp", "gp"))
  gates <- tibble::tibble(ind = "2", sample = "weak-1", response = "Weak",
    stage = c("Before", "After"), threshold = c(3, 0))
  scored <- env$.simClusterWeakScore(cells, gates)
  expect_equal(scored$nTruePos, c(1, 3))
  expect_equal(scored$nFalsePos, c(0, 1))
  expect_equal(scored$nPosStim, c(1, 4))
  expect_equal(scored$sensitivity, c(1 / 3, 1))
  expect_equal(scored$fdp, c(0, 1 / 4))
  expect_equal(scored$false_positive_rate, c(0, 1 / 2))
  changes <- env$.simClusterWeakChanges(scored)
  expect_equal(changes$responders_recovered, 2)
  expect_equal(changes$additional_false_positives, 1)
  expect_equal(changes$threshold_change, -3)
  expect_equal(changes$sensitivity_change, 2 / 3)
  expect_error(env$.simClusterWeakChanges(scored[1, ]), "pair")

  empty <- env$.simClusterWeakScore(cells, transform(gates[1, ], threshold = 10))
  expect_equal(empty$sensitivity, 0)
  expect_true(is.na(empty$fdp))
  missing <- env$.simClusterWeakScore(cells, transform(gates[1, ], threshold = NA_real_))
  expect_true(is.na(missing$nPosStim))
  expect_true(is.na(missing$sensitivity))
  negative_cells <- transform(cells, label = "gn")
  expect_true(is.na(env$.simClusterWeakScore(negative_cells, gates)$sensitivity[[1]]))
})

test_that("weak-response simulation preserves labels, seed and applied counts", {
  testthat::skip_if_not_installed("simcyto")
  env <- .weak_test_env()
  settings <- list(seed = 559L, n_strong = 3L, n_weak = 1L,
    strong_prob = 0.01, weak_prob = 0.002, n_cell = 500L,
    mean_pos = 4.5, variance = 1.5, bw = 0.25, bias_uns = 0.25)
  withr::local_seed(31)
  rng <- .Random.seed
  data <- env$.analysis_with_seed(settings$seed, env$.simClusterWeakCells(settings))
  data_again <- env$.analysis_with_seed(settings$seed, env$.simClusterWeakCells(settings))
  expect_identical(data$cells, data_again$cells)
  expect_identical(.Random.seed, rng)
  expect_equal(nrow(data$cells), 2 * 4 * 500)
  expect_equal(data$settings_table$response_probability, c(0.01, 0.01, 0.01, 0.002))
  expect_true(all(data$cells$label[data$cells$condition == "Unstimulated"] == "gn"))
  expect_equal(sum(data$cells$label[data$cells$sample == "weak-1"] == "gp"), 1L)
  expect_error(env$.simClusterWeakCells(utils::modifyList(settings, list(weak_prob = 0.02))), "weak_prob")

  # Run the real public gate workflow; do not require it to lower a threshold.
  withr::local_envvar(c(STIMGATE_INTERMEDIATE = "2"))
  before_env <- Sys.getenv("STIMGATE_INTERMEDIATE")
  result <- NULL
  utils::capture.output(result <- suppressMessages(env$.simClusterWeakRun(settings)))
  expect_identical(Sys.getenv("STIMGATE_INTERMEDIATE"), before_env)
  expect_identical(.Random.seed, rng)
  expect_equal(nrow(result$allocations), 4L)
  expect_equal(nrow(result$scores), 8L)
  expect_equal(result$scores$nPosStim, result$scores$package_nPosStim)
  expect_equal(result$scores$nTruePos + result$scores$nFalsePos, result$scores$nPosStim)
  expect_identical(result$cells, data$cells)
  expect_silent(env$.simClusterWeakValidate(result, settings))
  altered <- result
  altered$scores$nTruePos[[1]] <- altered$scores$nTruePos[[1]] + 1L
  expect_error(env$.simClusterWeakValidate(altered, settings), "counts")
  expect_error(env$.simClusterWeakValidate(result, utils::modifyList(settings, list(seed = 560L))), "settings")

  plots <- list(env$.simClusterWeakDistributions(result),
    env$.simClusterWeakDistributions(result, tail = TRUE), env$.simClusterWeakOutcomes(result))
  for (plot in plots) {
    expect_s3_class(plot, "ggplot")
    expect_null(plot$labels$title)
    expect_null(plot$labels$subtitle)
    expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  }
  # Overlapping before/after lines remain separately grouped.
  lines <- ggplot2::ggplot_build(plots[[1]])$data[[2]]
  expect_equal(nrow(lines), 8L)
})

test_that("weak QMD guards all optional tables and figures", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  qmd <- readLines(file.path(root, "analysis", "12-sim-cluster-gates.qmd"))
  env <- new.env(parent = baseenv())
  env$run_plots <- FALSE
  for (label in c("weak-tables", "weak-distributions", "weak-outcomes")) {
    start <- match(paste0("#| label: ", label), qmd)
    expect_false(is.na(start))
    end <- start + which(qmd[(start + 1L):length(qmd)] == "```")[[1]]
    code <- qmd[(start + 1L):(end - 1L)]
    expect_length(utils::capture.output(eval(parse(text = code), env)), 0L)
  }
  expect_true(any(grepl("weak_prob = 0.002", qmd, fixed = TRUE)))
  expect_true(any(grepl(".write_rds_atomic", qmd, fixed = TRUE)))
})
