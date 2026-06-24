test_that("stimgate_gate runs", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")
  invisible(stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))
  expect_true(file.exists(file.path(path_project, "gate_stats.rds")))
})

test_that("stimgate_gate runs with STIMGATE_DEBUG environment variable", {
  # Run stimgate_gate with STIMGATE_DEBUG environment variable to ensure it works without error
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate_debug")

  # Set the environment variable for debug mode
  Sys.setenv("STIMGATE_DEBUG" = "TRUE")
  on.exit(Sys.unsetenv("STIMGATE_DEBUG"))

  # Test that the function can be called with debug enabled without error
  # and that it returns the expected path
  result_path <- stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  )

  # Verify the function returns the expected path
  expect_equal(result_path, path_project)

  # Verify basic output files still exist (same as non-debug version)
  expect_true(file.exists(file.path(path_project, "gate_stats.rds")))

  # Clean up debug test directory
  if (dir.exists(path_project)) {
    unlink(path_project, recursive = TRUE)
  }
})

test_that("stimgate_gate runs with STIMGATE_INTERMEDIATE environment variable", {
  # Run stimgate_gate with STIMGATE_INTERMEDIATE environment variable to ensure
  # intermediate data is saved correctly
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(
    dirname(example_data$path_gs),
    "stimgate_intermediate"
  )

  # Set the environment variable for intermediate data saving
  Sys.setenv("STIMGATE_INTERMEDIATE" = "true")
  on.exit(Sys.unsetenv("STIMGATE_INTERMEDIATE"))

  # Test that the function can be called with intermediate saving enabled
  result_path <- stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  )

  # Verify the function returns the expected path
  expect_equal(result_path, path_project)

  # Verify basic output files exist
  expect_true(file.exists(file.path(path_project, "gate_stats.rds")))

  # Verify intermediate data directory exists
  intermediate_dir <- file.path(path_project, "intermediate_data")
  expect_true(dir.exists(intermediate_dir))

  # Verify that intermediate files were created
  intermediate_files <- list.files(intermediate_dir, recursive = TRUE)
  expect_true(length(intermediate_files) > 0)

  # Verify that files exist for the "init" stage
  init_files <- list.files(
    file.path(intermediate_dir, "init"),
    recursive = TRUE
  )
  expect_true(length(init_files) > 0)

  # Clean up intermediate test directory
  if (dir.exists(path_project)) {
    unlink(path_project, recursive = TRUE)
  }
})
