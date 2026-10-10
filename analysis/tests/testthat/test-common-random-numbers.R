test_that("full simulation grids seed biological scenarios rather than tuning rows", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  for (id in c("2a", "2b", "3", "4", "5", "6", "7", "8")) {
    env <- new.env(parent = getNamespace("stimgate"))
    for (fn in c("analysis-runtime.R", "sim-misc.R", "sim-bandwidth.R",
                 "sim-bandwidth-analysis-run.R", "sim-compare-freq_bs.R",
                 "sim-compare-qmd7-presentation.R")) {
      source(file.path(root, "scripts", "r", fn), local = env)
    }
    env$root_dir <- root
    env$run_simulations <- FALSE
    env$analysis_dev <- env$analysis_quick <- FALSE
    env$sim_size <- "final"
    env$simulation_seed <- 12345L
    env$sim_grid_shuffle_seed <- 8L
    env$sim_grid_chunk_index <- env$sim_grid_n_chunks <- 1L
    # The adaptive-bandwidth analyses 5 and 6 are archived.
    dir <- file.path(root, "analysis")
    if (id %in% c("5", "6")) dir <- file.path(dir, "_archive")
    lines <- readLines(list.files(dir, paste0("^", id, "-.*qmd$"), full.names = TRUE))
    chunk <- function(label) {
      start <- which(lines == paste0("#| label: ", label))
      end <- start + which(lines[(start + 1L):length(lines)] == "```")[1L]
      lines[(start + 1L):(end - 1L)]
    }
    eval(parse(text = chunk(if (id == "2b") "bias-uns-settings" else "actual-settings")), env)
    label <- switch(id, `2a` = "bw-manual-grid", `2b` = "bias-uns-grid",
      `6` = "bw-manual-grid", `7` = "compare-grid", `8` = "mismatch-grids",
      "bw-estimate-grid")
    # Only grid construction is needed, before diagnostics/filtering/shuffling.
    code <- chunk(label)
    end <- which(grepl("sim_grid_full <- sim_grid", code, fixed = TRUE))[1L]
    if (!is.na(end)) code <- code[seq_len(end)]
    invisible(utils::capture.output(eval(parse(text = code), env)))
    grid <- env$sim_grid_full
    if (is.null(grid)) grid <- env$sim_grid
    bio <- intersect(c("transformation", "mean_pos", "prob_response", "n_cell",
      "sample_perturbation_sd", "condition_perturbation_sd", "cluster_perturbation_sd",
      "background_relative_to_response", "n_cell_uns_relative_to_stim"), names(grid))
    # Analysis 8 deliberately shares baseline draws across deterministic mismatches.
    groups <- grid |>
      dplyr::group_by(dplyr::across(dplyr::all_of(bio))) |>
      dplyr::summarise(n_seed = dplyr::n_distinct(.data$sim_seed),
        seed = dplyr::first(.data$sim_seed), .groups = "drop")
    expect_true(all(groups$n_seed == 1L), info = id)
    expect_equal(dplyr::n_distinct(groups$seed), nrow(groups), info = id)
    expect_equal(grid$sim_seed, 12345L + grid$base_scenario_id - 1L, info = id)
  }
})

test_that("real bandwidth wrappers pair every tube including later replicates", {
  withr::local_preserve_seed()
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c("analysis-runtime.R", "sim-misc.R", "sim-bandwidth.R",
               "sim-bandwidth-analysis-run.R")) {
    source(file.path(root, "scripts", "r", fn), local = env)
  }
  original <- simcyto::simCytExperiment
  tubes <- list()
  capture <- function(...) {
    out <- original(...)
    tubes[[length(tubes) + 1L]] <<- lapply(out$flowFrameList, flowCore::exprs)
    out
  }
  settings <- list(nSample = 2L, nMarker = 1L, nCondition = 2L,
    nCluster = 2L, nIter = 2L, probExact = TRUE)
  row <- tibble::tibble(sim_id = 1L, base_scenario_id = 1L, sim_seed = 17L,
    transformation = "gaussian", n_cell = 200L, prob_response = 0.2, mean_pos = 8,
    bw = 0.2, bias_uns = 0.05, sample_perturbation_sd = 0,
    condition_perturbation_sd = 0, cluster_perturbation_sd = 0,
    background_relative_to_response = 0.2, n_cell_uns_relative_to_stim = 1)
  testthat::with_mocked_bindings(simCytExperiment = capture, .package = "simcyto", {
    for (bw in c(0.2, 0.5)) {
      row$bw <- bw
      env$.simBandwidthRunRow(row, env$.simBandwidthFreqBsGlobalScenario, settings)
    }
    expect_length(tubes, 4L)
    expect_identical(tubes[[1]], tubes[[3]])
    expect_identical(tubes[[2]], tubes[[4]])
    expect_false(identical(tubes[[1]], tubes[[2]]))
    tubes <- list()
    # Deliberately consume different RNG amounts in the estimator to expose drift.
    real_bw <- env$.simBandwidthBwOne
    env$.simBandwidthBwOne <- function(..., bwMtd) {
      stats::runif(if (bwMtd == "nrd0") 5L else 31L)
      real_bw(..., bwMtd = bwMtd)
    }
    for (mtd in c("nrd0", "sj")) {
      env$.analysis_with_seed(17L, env$.simBandwidthEstBwDirect(
        nSample = 2L, nIter = 2L, nCellStim = 200L, probResponse = 0.2,
        meanPos = 8, transformation = "gaussian", bwMtd = mtd, summarise = FALSE
      ))
    }
    expect_identical(tubes[[1]], tubes[[3]])
    expect_identical(tubes[[2]], tubes[[4]])
  })
})
