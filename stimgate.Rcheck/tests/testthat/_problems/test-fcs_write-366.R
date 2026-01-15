# Extracted from test-fcs_write.R:366

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "stimgate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(testthat)
example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)
path_project <- file.path(dirname(example_data$path_gs), "stimgate")
invisible(stimgate::stimgate_gate(
  .data = gs,
  path_project = path_project,
  pop_gate = "root",
  batch_list = example_data$batch_list,
  marker = example_data$marker
))

# test -------------------------------------------------------------------------
path_dir_save <- file.path(tempdir(), "fcs_output_exclusions")
if (length(example_data$marker) > 1) {
    combn_exc <- list(example_data$marker[[1]])

    result <- stimgate::stimgate_fcs_write(
      path_project = path_project,
      .data = gs,
      ind_batch_list = example_data$batch_list,
      path_dir_save = path_dir_save,
      chnl = example_data$marker,
      combn_exc = combn_exc
    )

    expect_equal(result, path_dir_save)
    expect_true(dir.exists(path_dir_save))
  } else {
    # Test with NULL exclusions
    result <- stimgate::stimgate_fcs_write(
      path_project = path_project,
      .data = gs,
      ind_batch_list = example_data$batch_list,
      path_dir_save = path_dir_save,
      chnl = example_data$marker[[1]],
      combn_exc = NULL
    )

    expect_equal(result, path_dir_save)
    expect_true(dir.exists(path_dir_save))
  }
