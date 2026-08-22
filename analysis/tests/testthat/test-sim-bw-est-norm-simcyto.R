root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_bw_io <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")

test_that("normalised bandwidth simulation helpers source cleanly without legacy functionsForBenchmarking-Cyt.R", {
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

test_that("analysis/4-sim-bw-est-norm.qmd does not source functionsForBenchmarking-Cyt.R", {
  qmd_path <- file.path(root_dir, "analysis", "4-sim-bw-est-norm.qmd")
  expect_true(file.exists(qmd_path))

  lines <- readLines(qmd_path, warn = FALSE)
  expect_false(
    any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
    info = "analysis/4-sim-bw-est-norm.qmd should not source functionsForBenchmarking-Cyt.R"
  )
})

test_that(".simBandwidthEstBwDirect normalised path has fixed-seed simcyto parity and stable ordinary-vs-normalised outputs", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)

  seed <- 2026L
  n_sample <- 2L
  n_condition <- 2L
  n_cell_stim <- 400L
  n_cell_uns <- 200L
  prob_response <- 0.05
  background_relative_to_response <- 0.2
  mean_pos <- 6

  common_args <- list(
    nSample = n_sample,
    nMarker = 1L,
    nCondition = n_condition,
    nCluster = 2L,
    nIter = 1L,
    biasUns = 0.05,
    bwFallback = 0.234,
    bwMin = -Inf,
    bwMax = Inf,
    bwNcellMax = 500L,
    bwCluster = 0.5,
    nCellStim = n_cell_stim,
    probResponse = prob_response,
    probExact = TRUE,
    meanPos = mean_pos,
    transformation = "skew",
    backgroundRelativeToResponse = background_relative_to_response,
    ncellUnsRelativeToStim = 0.5,
    covEvMin = 1.5,
    covEvMax = 1.5,
    tolClust = NULL,
    summarise = FALSE
  )

  captured_sim <- NULL
  captured_args <- NULL
  orig_simcyto_experiment <- simcyto::simCytExperiment

  set.seed(seed)
  res_norm <- testthat::with_mocked_bindings(
    simCytExperiment = function(...) {
      captured_args <<- list(...)
      out <- do.call(orig_simcyto_experiment, captured_args)
      captured_sim <<- out
      out
    },
    .package = "simcyto",
    do.call(env$.simBandwidthEstBwDirect, c(common_args, list(bwMtd = "hpi1Norm")))
  )

  expect_true(is.list(captured_args))
  expect_equal(captured_args$nCellByCondition, c(n_cell_uns, n_cell_stim))
  expect_equal(captured_args$meanExprMat, matrix(c(0, mean_pos), byrow = TRUE, ncol = 1))
  expect_equal(captured_args$clusterLabelVec, c("gn", "gp"))
  expect_equal(captured_args$probVecUns, c(0.99, 0.01), tolerance = 1e-12)
  expect_equal(captured_args$probResponseVecByStimCondition, list(c(-0.05, 0.05)))
  expect_true(is.function(captured_args$transformationFunc))
  expect_equal(attr(captured_args$transformationFunc, "sim_transformation"), "skew")
  expect_equal(captured_args$samplePerturbationSd, 0)
  expect_equal(captured_args$conditionPerturbationSd, 0)
  expect_equal(captured_args$clusterPerturbationSd, 0)

  set.seed(seed)
  sim_direct <- simcyto::simCytExperiment(
    nSample = n_sample,
    nMarker = 1L,
    nCondition = n_condition,
    nCluster = 2L,
    nCellByCondition = c(n_cell_uns, n_cell_stim),
    transformationFunc = env$.simMiscGetTrans("skew"),
    mixtureType = "gaussianOnly",
    meanExprMat = matrix(c(0, mean_pos), byrow = TRUE, ncol = 1),
    clusterLabelVec = c("gn", "gp"),
    probVecUns = c(0.99, 0.01),
    probExact = TRUE,
    probResponseVecByStimCondition = list(c(-0.05, 0.05)),
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    covEvMin = 1.5,
    covEvMax = 1.5
  )

  expect_type(captured_sim, "list")
  expect_true(all(c("flowFrameList", "labelsList") %in% names(captured_sim)))
  expect_length(captured_sim$flowFrameList, 4L)
  expect_length(captured_sim$labelsList, 4L)

  for (i in seq_along(captured_sim$flowFrameList)) {
    expr_helper <- flowCore::exprs(captured_sim$flowFrameList[[i]])[, "F1"]
    expr_direct <- flowCore::exprs(sim_direct$flowFrameList[[i]])[, "F1"]
    expect_equal(expr_helper, expr_direct, tolerance = 1e-12)
    expect_equal(captured_sim$labelsList[[i]], sim_direct$labelsList[[i]])
  }

  expr_means_helper <- vapply(
    captured_sim$flowFrameList,
    function(ff) mean(flowCore::exprs(ff)[, "F1"]),
    numeric(1)
  )
  expr_sds_helper <- vapply(
    captured_sim$flowFrameList,
    function(ff) stats::sd(flowCore::exprs(ff)[, "F1"]),
    numeric(1)
  )
  expr_means_direct <- vapply(
    sim_direct$flowFrameList,
    function(ff) mean(flowCore::exprs(ff)[, "F1"]),
    numeric(1)
  )
  expr_sds_direct <- vapply(
    sim_direct$flowFrameList,
    function(ff) stats::sd(flowCore::exprs(ff)[, "F1"]),
    numeric(1)
  )
  expect_equal(unname(expr_means_helper), unname(expr_means_direct), tolerance = 1e-12)
  expect_equal(unname(expr_sds_helper), unname(expr_sds_direct), tolerance = 1e-12)

  truth_from_helper <- purrr::map_dfr(seq_len(n_sample), function(sample_curr) {
    ind_uns <- (sample_curr - 1L) * n_condition + 1L
    ind_stim <- ind_uns + 1L
    labels_uns <- captured_sim$labelsList[[ind_uns]]
    labels_stim <- captured_sim$labelsList[[ind_stim]]

    tibble::tibble(
      sample = as.character(sample_curr),
      prop_uns_truth = sum(labels_uns == "gp") / length(labels_uns),
      prop_stim_truth = sum(labels_stim == "gp") / length(labels_stim),
      n_uns = length(labels_uns),
      n_stim = length(labels_stim)
    )
  })

  expect_true(all(truth_from_helper$prop_uns_truth >= 0 & truth_from_helper$prop_uns_truth <= 1))
  expect_true(all(truth_from_helper$prop_stim_truth >= 0 & truth_from_helper$prop_stim_truth <= 1))
  expect_equal(truth_from_helper$n_uns, rep(n_cell_uns, n_sample))
  expect_equal(truth_from_helper$n_stim, rep(n_cell_stim, n_sample))

  expect_s3_class(res_norm, "tbl_df")
  expect_equal(nrow(res_norm), n_sample)
  expect_equal(res_norm$n_cell_uns, rep(n_cell_uns, n_sample))
  expect_equal(res_norm$n_cell_stim, rep(n_cell_stim, n_sample))

  set.seed(seed)
  res_hpi1 <- do.call(env$.simBandwidthEstBwDirect, c(common_args, list(bwMtd = "hpi1")))

  expect_true(all(is.finite(res_hpi1$bw)))
  expect_true(all(is.finite(res_norm$bw)))
  expect_true(all(res_hpi1$bw > 0))
  expect_true(all(res_norm$bw > 0))

  expect_equal(res_hpi1$bw, c(0.529433696153644, 0.48829492916954), tolerance = 1e-12)
  expect_equal(res_norm$bw, c(0.570753196215079, 0.580931838697578), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(res_hpi1$bw, res_norm$bw, tolerance = 0)))
})
