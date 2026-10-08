root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

test_that("analysis 7 uses run-specific progress and validates full nested collation", {
  qmd_path <- file.path(root_dir, "analysis", "7-sim-compare-freq_bs.qmd")
  expect_true(file.exists(qmd_path))

  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl('stimgate_bw_scope <- "cytokine"', content, fixed = TRUE))
  expect_true(grepl("bw_scope = stimgate_bw_scope", content, fixed = TRUE))
  expect_true(grepl("stimgate_bw_scope = stimgate_bw_scope", content, fixed = TRUE))
  expect_true(grepl("stimgate_bw_ncell_min = stimgate_bw_ncell_min", content, fixed = TRUE))
  expect_true(grepl("simulation_seed:\\s*12345", content))
  expect_true(grepl(
    'comparison_semantics_version <- "corrected-comparison-v21"',
    content, fixed = TRUE
  ))
  expect_false(grepl("sim_grid_shuffle_seed", content, fixed = TRUE))
  helper_content <- paste(readLines(script_comp), collapse = "\n")
  expect_true(grepl(".simComparePromoteIfReady(", content, fixed = TRUE))
  expect_true(grepl(
    "path_progress_file <- run_ctx$progress_file",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("recursive\\s*=\\s*TRUE", helper_content))
  expect_true(grepl("pathList\\s*=\\s*scenario_paths", helper_content))
  expect_true(grepl("expected_sim_ids", helper_content, fixed = TRUE))
  expect_true(grepl(
    "Refusing to promote comparison",
    helper_content,
    fixed = TRUE
  ))
  expect_true(grepl(".analysis_current_file", content, fixed = TRUE))
  expect_true(grepl(
    "required_params = analysis_result_params",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    'normal.kind = "Inversion", sample.kind = "Rejection"',
    helper_content,
    fixed = TRUE
  ))
  expect_false(grepl(
    ".simCompareRunScenarioUnseeded",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(".analysis_results_context", content, fixed = TRUE))
  expect_true(grepl(
    ".simCompareGridOutputStatus",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("failed_ids", content, fixed = TRUE))
  expect_true(grepl("missing_ids", content, fixed = TRUE))
  expect_true(grepl("extra_ids", content, fixed = TRUE))
  expect_true(grepl("analysis_dev", content, fixed = TRUE))
  expect_false(grepl("bias_uns == 0.05", content, fixed = TRUE))
  expect_true(grepl("knitr::knit_exit()", content, fixed = TRUE))
  expect_true(grepl(
    '.analysis_require_packages(c("cytoUtils", "openCyto", "flowStats", "simcyto"))',
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "F-beta comparison script not found:",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "condition_perturbation_sd == 0",
    content,
    fixed = TRUE
  ))
  # Threshold densities are drawn by the QMD 7 plot helper.
  expect_true(grepl("after_stat(density)", helper_content, fixed = TRUE))
  expect_true(grepl(
    "Median absolute relative error",
    content,
    fixed = TRUE
  ))
  expect_false(grepl("#| error: true", content, fixed = TRUE))
})

test_that("analysis 7 nested chunk outputs can be collated without simulations", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_comp, local = env)

  tmp <- withr::local_tempdir()
  chunk_1 <- file.path(tmp, "chunks", "001-of-002", "output")
  chunk_2 <- file.path(tmp, "chunks", "002-of-002", "output")
  dir.create(chunk_1, recursive = TRUE)
  dir.create(chunk_2, recursive = TRUE)

  out_1 <- tibble::tibble(
    sim_id = 1L,
    method = "stimgate",
    propRespEst = 0.05,
    propRespTruth = 0.05,
    error = NA_character_
  )
  out_2 <- tibble::tibble(
    sim_id = 2L,
    method = "stimgate",
    propRespEst = 0.04,
    propRespTruth = 0.05,
    error = NA_character_
  )

  saveRDS(
    out_1,
    file.path(
      chunk_1,
      "compare_raw-chunk_001-of_002-sim_id_000001.rds"
    )
  )
  saveRDS(
    out_2,
    file.path(
      chunk_2,
      "compare_raw-chunk_002-of_002-sim_id_000002.rds"
    )
  )

  scenario_paths <- list.files(
    tmp,
    pattern = "^(compare_raw.*|sim_scenario.*|sim_raw.*)sim_id_[0-9]+[.]rds$",
    recursive = TRUE,
    full.names = TRUE
  )

  expect_length(scenario_paths, 2L)

  collated <- env$.simCompareCollateScenarioOutputs(
    pathList = scenario_paths,
    sim_grid = tibble::tibble(sim_id = c(1L, 2L))
  )

  expect_equal(sort(unique(collated$sim_id)), c(1L, 2L))
  expect_equal(nrow(collated), 2L)
})

test_that("shared comparison promotion validates all nested scenario outputs", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_comp, local = env)
  tmp <- withr::local_tempdir()
  ctx <- list(
    staging_run_dir = tmp,
    staging_collated_dir = file.path(tmp, "collated")
  )
  dir.create(ctx$staging_collated_dir)
  can_promote <- FALSE
  promoted <- FALSE
  env$.analysis_can_promote <- function(...) can_promote
  env$.analysis_promote_run <- function(...) { promoted <<- TRUE; TRUE }
  env$.analysis_mark_chunk <- function(...) invisible(NULL)
  env$.write_rds_atomic <- function(object, path) saveRDS(object, path)
  grid <- tibble::tibble(sim_id = c(1L, 2L))
  promote <- function() env$.simComparePromoteIfReady(
    ctx, grid, total_sims = 1L, completed_sims = 1L, failed_sims = 0L,
    nSample = 1L, nIter = 1L
  )
  expect_false(promote())
  expect_false(promoted)
  can_promote <- TRUE
  for (id in grid$sim_id) {
    chunk <- file.path(tmp, "chunks", as.character(id), "output")
    dir.create(chunk, recursive = TRUE)
    result <- tibble::tibble(
      sim_id = id, iter = 1L, sample = 1L,
      method = c("stimgate", "fbeta", "tailgate"),
      propRespEst = 0.05, propRespTruth = 0.05, error = NA_character_,
      nCellStim = 100, nPosStim = 5L, nTruePos = 5L, nFalsePos = 0L,
      nFalseNeg = 0L, nTrueNeg = 95L, unsExprSum = 1
    )
    saveRDS(result, file.path(chunk, sprintf("sim_raw-sim_id_%06d.rds", id)))
  }
  ctx$read_only <- TRUE
  expect_false(promote())
  expect_false(promoted)
  ctx$read_only <- FALSE
  expect_true(promote())
  expect_true(promoted)
  expect_equal(nrow(readRDS(file.path(ctx$staging_collated_dir,
    "compare_raw.rds"))), 6L)
  promoted <- FALSE
  result$error[1] <- "failure"
  saveRDS(result, file.path(chunk, sprintf("sim_raw-sim_id_%06d.rds", id)))
  expect_error(promote(), "Refusing to promote comparison")
  expect_false(promoted)
  result$error <- NA_character_
  result$propRespEst[1] <- NA_real_
  saveRDS(result, file.path(chunk, sprintf("sim_raw-sim_id_%06d.rds", id)))
  expect_error(promote(), "Refusing to promote comparison")
  expect_false(promoted)
  unlink(file.path(chunk, sprintf("sim_raw-sim_id_%06d.rds", id)))
  expect_error(promote(), "Refusing to promote comparison")
})

test_that("rerun chunks use the production scientific arguments", {
  find_calls <- function(expr, name) {
    if (missing(expr)) return(list())
    if (!is.call(expr) && !is.expression(expr) && !is.pairlist(expr)) {
      return(list())
    }
    if (is.call(expr) && identical(expr[[1]], as.name(name))) {
      return(list(expr))
    }
    unlist(lapply(as.list(expr), find_calls, name = name), recursive = FALSE)
  }
  for (qmd_name in c(
    "7-sim-compare-freq_bs.qmd", "8-sim-compare-freq_bs-batch.qmd"
  )) {
    lines <- readLines(file.path(root_dir, "analysis", qmd_name))
    starts <- which(lines == "```{r}")
    chunks <- lapply(starts, function(start) {
      end <- which(seq_along(lines) > start & lines == "```")[1]
      parse(text = lines[seq.int(start + 1L, end - 1L)])
    })
    grid_calls <- unlist(lapply(chunks, find_calls,
      name = ".simCompareFreqBsGrid"), recursive = FALSE)
    rerun_calls <- unlist(lapply(chunks, find_calls,
      name = ".simCompareRunScenario"), recursive = FALSE)
    expect_length(grid_calls, 1L)
    expect_length(rerun_calls, 1L)
    grid_args <- as.list(grid_calls[[1]])[-1]
    rerun_args <- as.list(rerun_calls[[1]])[-1]
    scientific <- setdiff(names(rerun_args), c("row", "resume"))
    expect_identical(grid_args[scientific], rerun_args[scientific])
  }
})
