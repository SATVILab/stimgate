root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

script_comp <- file.path(root_dir, "scripts", "r", "sim-compare-freq_bs.R")

test_that("analysis 7 uses run-specific progress and validates full nested collation", {
  qmd_path <- file.path(root_dir, "analysis", "7-sim-compare-freq_bs.qmd")
  expect_true(file.exists(qmd_path))

  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("simulation_seed:\\s*1", content))
  expect_true(grepl("comparison_semantics_version", content, fixed = TRUE))
  expect_true(grepl(
    "sim_grid_shuffle_seed\\s*=\\s*simulation_seed",
    content
  ))
  expect_true(grepl(
    "path_progress_file <- run_ctx$progress_file",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("recursive\\s*=\\s*TRUE", content))
  expect_true(grepl("pathList\\s*=\\s*scenario_paths", content))
  expect_true(grepl("expected_sim_ids", content, fixed = TRUE))
  expect_true(grepl(
    "Refusing to promote analysis 7",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(".analysis_current_file", content, fixed = TRUE))
  expect_true(grepl(
    "required_params = list(",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "set.seed(as.integer(row$sim_seed[[1]]))",
    content,
    fixed = TRUE
  ))
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
