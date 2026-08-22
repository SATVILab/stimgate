root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_bw_io <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")

test_that("adaptive bandwidth simulation helpers source cleanly without legacy functionsForBenchmarking-Cyt.R", {
  for (f in c(script_misc, script_bw, script_bw_io, script_bw_plot)) {
    if (!file.exists(f)) stop("Expected analysis helper not found: ", f)
  }

  env <- new.env(parent = getNamespace("stimgate"))
  expect_no_error(source(script_misc, local = env))
  expect_no_error(source(script_bw, local = env))
  expect_no_error(source(script_bw_io, local = env))
  expect_no_error(source(script_bw_plot, local = env))

  expect_false(exists("simCytExperiment", envir = env, inherits = FALSE))
})

test_that("analysis/5-sim-bw-est-adaptive.qmd does not source functionsForBenchmarking-Cyt.R", {
  qmd_path <- file.path(root_dir, "analysis", "5-sim-bw-est-adaptive.qmd")
  expect_true(file.exists(qmd_path))

  lines <- readLines(qmd_path, warn = FALSE)
  expect_false(
    any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
    info = "analysis/5-sim-bw-est-adaptive.qmd should not source functionsForBenchmarking-Cyt.R"
  )
})

test_that(".simBandwidthEstBwDirectAdaptive preserves simcyto simulation boundary and adaptive outputs", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)

  run_case <- function(
    seed,
    transformation,
    mean_pos,
    bias_uns,
    expected_bw_means
  ) {
    n_sample <- 2L
    n_condition <- 2L
    n_cell_stim <- 240L
    n_cell_uns <- 240L
    prob_response <- 0.05
    background_relative_to_response <- 0.2
    prob_response_uns <- prob_response * background_relative_to_response

    captured_args <- NULL
    captured_sim <- NULL
    orig_simcyto_experiment <- simcyto::simCytExperiment

    set.seed(seed)
    res <- testthat::with_mocked_bindings(
      simCytExperiment = function(...) {
        captured_args <<- list(...)
        out <- orig_simcyto_experiment(...)
        captured_sim <<- out
        out
      },
      .package = "simcyto",
      env$.simBandwidthEstBwDirectAdaptive(
        nSample = n_sample,
        nMarker = 1L,
        nCondition = n_condition,
        nCluster = 2L,
        nIter = 1L,
        biasUns = bias_uns,
        bwMtd = "hpi1Norm",
        bwFallback = 0.234,
        bwMin = -Inf,
        bwMax = Inf,
        bwNcellMax = 500L,
        nCellStim = n_cell_stim,
        probResponse = prob_response,
        probExact = TRUE,
        meanPos = mean_pos,
        transformation = transformation,
        backgroundRelativeToResponse = background_relative_to_response,
        ncellUnsRelativeToStim = 1,
        covEvMin = 1.5,
        covEvMax = 1.5,
        summarise = FALSE
      )
    )

    expect_s3_class(res, "tbl_df")
    expect_equal(nrow(res), n_sample)

    expect_type(captured_args, "list")
    expect_equal(captured_args$nCellByCondition, c(n_cell_uns, n_cell_stim))
    expect_equal(captured_args$meanExprMat, matrix(c(0, mean_pos), byrow = TRUE, ncol = 1))
    expect_equal(captured_args$clusterLabelVec, c("gn", "gp"))
    expect_equal(captured_args$probVecUns, c(1 - prob_response_uns, prob_response_uns))
    expect_equal(captured_args$probResponseVecByStimCondition, list(c(-prob_response, prob_response)))
    expect_equal(captured_args$samplePerturbationSd, 0)
    expect_equal(captured_args$conditionPerturbationSd, 0)
    expect_equal(captured_args$clusterPerturbationSd, 0)
    expect_equal(captured_args$covEvMin, 1.5)
    expect_equal(captured_args$covEvMax, 1.5)

    expected_trans <- env$.simBandwidthGetTrans(transformation)
    expect_equal(captured_args$transformationFunc(c(-1, 0, 1)), expected_trans(c(-1, 0, 1)))

    set.seed(seed)
    sim_direct <- do.call(simcyto::simCytExperiment, captured_args)

    truth_from_labels <- function(sim_obj) {
      purrr::map_df(seq_len(n_sample), function(sample_curr) {
        ind_uns <- (sample_curr - 1L) * n_condition + 1L
        ind_stim <- ind_uns + 1L
        labels_uns <- sim_obj$labelsList[[ind_uns]]
        labels_stim <- sim_obj$labelsList[[ind_stim]]
        prop_uns_truth <- sum(labels_uns %in% "gp") / length(labels_uns)
        prop_stim_truth <- sum(labels_stim %in% "gp") / length(labels_stim)

        tibble::tibble(
          sample = as.character(sample_curr),
          ind = as.character(ind_stim),
          propStimTruth = prop_stim_truth,
          propUnsTruth = prop_uns_truth,
          propRespTruth = prop_stim_truth - prop_uns_truth
        )
      }) |>
        dplyr::arrange(.data$sample, .data$ind)
    }

    expect_equal(truth_from_labels(captured_sim), truth_from_labels(sim_direct), tolerance = 1e-12)

    expr_means_helper <- vapply(
      captured_sim$flowFrameList,
      function(ff) mean(flowCore::exprs(ff)[, 1]),
      numeric(1)
    )
    expr_sds_helper <- vapply(
      captured_sim$flowFrameList,
      function(ff) stats::sd(flowCore::exprs(ff)[, 1]),
      numeric(1)
    )
    expr_means_direct <- vapply(
      sim_direct$flowFrameList,
      function(ff) mean(flowCore::exprs(ff)[, 1]),
      numeric(1)
    )
    expr_sds_direct <- vapply(
      sim_direct$flowFrameList,
      function(ff) stats::sd(flowCore::exprs(ff)[, 1]),
      numeric(1)
    )
    expect_equal(unname(expr_means_helper), unname(expr_means_direct), tolerance = 1e-12)
    expect_equal(unname(expr_sds_helper), unname(expr_sds_direct), tolerance = 1e-12)

    expect_equal(res$n_cell_stim, rep(n_cell_stim, n_sample))
    expect_equal(res$n_cell_uns, rep(n_cell_uns, n_sample))
    expect_true(all(res$n_uns_bw_core <= n_cell_uns))
    expect_true(all(res$n_stim_bw_core <= n_cell_stim))

    bw_cols <- c("bw_uns_core", "bw_stim_core", "bw_uns_extra", "bw_stim_extra")
    expect_true(all(vapply(res[bw_cols], function(x) all(is.finite(x)), logical(1))))
    expect_true(all(vapply(res[bw_cols], function(x) all(x > 0), logical(1))))

    bw_means <- colMeans(res[bw_cols], na.rm = TRUE)
    expect_equal(unname(bw_means), expected_bw_means, tolerance = 1e-7)
  }

  run_case(
    seed = 2026L,
    transformation = "gaussian",
    mean_pos = 8,
    bias_uns = 0.05,
    expected_bw_means = c(
      0.3302348212,
      0.3571534616,
      0.5248110771,
      0.4150524274
    )
  )

  run_case(
    seed = 2027L,
    transformation = "gamma",
    mean_pos = 4,
    bias_uns = 0.0025,
    expected_bw_means = c(
      0.0050264283,
      0.0070008414,
      0.0058989428,
      0.0238334830
    )
  )
})
