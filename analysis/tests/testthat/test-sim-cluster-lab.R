test_that("lab shifts move both tubes and retain paired control definitions", {
  testthat::skip_if_not_installed("simcyto")
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "sim-cluster-lab.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  withr::local_seed(19)
  old_seed <- .Random.seed
  base <- env$.analysis_with_seed(558L, env$.simClusterLabData(2L, 40L, 0))
  shifted <- env$.analysis_with_seed(558L, env$.simClusterLabData(2L, 40L, 4))
  expect_identical(.Random.seed, old_seed)
  expect_length(shifted$matrices, 8L)
  expect_equal(shifted$batch_list, base$batch_list)
  expect_equal(unname(shifted$batch_list), list(c(1L, 2L), c(3L, 4L), c(5L, 6L), c(7L, 8L)))
  for (i in seq_along(base$matrices)) {
    expect_equal(shifted$matrices[[i]] - base$matrices[[i]],
      matrix(if (i > 4L) 4 else 0, nrow = 40L, ncol = 1L,
        dimnames = list(NULL, "Marker")))
  }
  expect_equal(nrow(shifted$expression), 320L)
  expect_error(env$.simClusterLabData(1L, 40L, 4))
  expect_error(env$.simClusterLabData(2L, 40L, -4))
})

test_that("lab summaries distinguish separated, mixed and missing allocations", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  env <- new.env(parent = getNamespace("stimgate"))
  source(file.path(root, "scripts", "r", "sim-cluster-lab.R"), local = env)
  samples <- tibble::tibble(
    sample = letters[1:4], lab = c("Lab A", "Lab A", "Lab B", "Lab B"),
    ind = as.character(2L * seq_len(4L))
  )
  details <- tibble::tibble(
    ind = samples$ind, grp = c("7", "9", "2", "2"),
    detailObject = "locClusterQuantileTbl", detailPathStage = "init",
    cpOrigQuantMin = c(14, 15, 18, 19), cpJoinTgOrig = c(14, 15, 18, 18),
    locGeneratedDirect = TRUE, locClusterAction = "direct_retained",
    locClusterReason = "cluster_direct_threshold_quantile_transfer"
  )
  allocations <- env$.simClusterLabAllocations(samples, details)
  expect_equal(allocations$threshold_original, details$cpOrigQuantMin)
  expect_equal(allocations$threshold_adjusted, details$cpJoinTgOrig)
  summary <- env$.simClusterLabSummary(allocations)
  expect_true(summary$labs_separated)
  expect_equal(summary$n_clusters, 3L) # More than one cluster in Lab A is allowed.
  allocations$cluster[[4L]] <- "7"
  expect_false(env$.simClusterLabSummary(allocations)$labs_separated)
  expect_equal(env$.simClusterLabSummary(allocations)$n_mixed_lab_clusters, 1L)
  incomplete <- env$.simClusterLabAllocations(samples, details[-1L, ])
  expect_equal(nrow(incomplete), nrow(samples))
  expect_equal(env$.simClusterLabSummary(incomplete)$n_unassigned, 1L)
  expect_true(is.na(env$.simClusterLabSummary(incomplete)$labs_separated))
  absent <- env$.simClusterLabAllocations(samples, tibble::tibble())
  expect_equal(env$.simClusterLabSummary(absent)$n_unassigned, 4L)
  expect_error(env$.simClusterLabAllocations(samples, dplyr::bind_rows(details, details[1L, ])),
    "one initial cluster row")
})

test_that("lab runner reads current package details and restores temporary state", {
  testthat::skip_if_not_installed("simcyto")
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "sim-cluster-lab.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  withr::local_seed(12)
  withr::local_envvar(c(STIMGATE_INTERMEDIATE = "unchanged"))
  before_env <- Sys.getenv("STIMGATE_INTERMEDIATE")
  before_seed <- .Random.seed
  result <- NULL
  utils::capture.output(suppressMessages(
    result <- env$.simClusterLabRun(558L, 2L, 200L, 4)
  ))
  expect_identical(Sys.getenv("STIMGATE_INTERMEDIATE"), before_env)
  expect_identical(.Random.seed, before_seed)
  expect_equal(nrow(result$allocations), 4L)
  expect_setequal(result$allocations$ind, c("2", "4", "6", "8"))
  expect_equal(result$summary$n_assigned + result$summary$n_unassigned, 4L)
  # getStimGates() returns the original ("min") and clustered ("minClust") rows.
  expect_equal(dplyr::n_distinct(result$final_gates$ind), 4L)
  cluster_details <- result$details |>
    dplyr::filter(.data$detailObject == "locClusterQuantileTbl",
      .data$detailPathStage == "init")
  if (nrow(cluster_details) > 0L) {
    expect_equal(result$summary$n_assigned, sum(!is.na(cluster_details$grp)))
  }
  # Lab separation is the result of the demonstration, not a test assumption.
})

test_that("lab plot builders return plots and QMD honours disabled plotting", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "sim-cluster-lab.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  expression <- tibble::tibble(
    sample = rep(c("a", "b"), each = 40L),
    lab = rep(c("Lab A", "Lab B"), each = 40L),
    condition = rep(rep(c("Unstimulated", "Stimulated"), each = 20L), 2L),
    expression = rep(seq(8, 12, length.out = 40L), 2L) + rep(c(0, 4), each = 40L)
  )
  allocations <- tibble::tibble(sample = c("a", "b"), lab = c("Lab A", "Lab B"),
    cluster = c("1", NA_character_))
  plots <- list(env$.simClusterLabDistributionPlot(expression),
    env$.simClusterLabAllocationPlot(allocations))
  for (plot in plots) {
    expect_s3_class(plot, "ggplot")
    expect_null(plot$labels$title)
    expect_null(plot$labels$subtitle)
    expect_no_error(ggplot2::ggplot_build(plot))
  }
  expect_true("Unassigned" %in% plots[[2L]]$data$cluster_display)
  lines <- readLines(file.path(root, "analysis", "12-sim-cluster-gates.qmd"))
  starts <- which(startsWith(lines, "```{r"))
  env$run_plots <- FALSE
  env$run_simulations <- FALSE
  env$.analysis_print_save_fig <- function(...) stop("Unexpected figure output")
  env$.analysis_report_table <- function(...) stop("Unexpected table output")
  env$.simClusterLabRun <- function(...) stop("Unexpected simulation")
  for (start in starts) {
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1L]
    code <- parse(text = lines[(start + 1L):(end - 1L)])
    label <- lines[start + 1L]
    if (label %in% paste0("#| label: ", c("lab-simulate", "lab-distributions", "lab-allocations"))) {
      expect_no_error(eval(code, envir = env))
    }
  }
})
