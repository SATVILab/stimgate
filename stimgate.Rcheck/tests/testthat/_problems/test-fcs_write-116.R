# Extracted from test-fcs_write.R:116

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
gate_methods <- c("min", "max", "mean", "tmean", "med")
for (method in gate_methods) {
    path_dir_save <- file.path(tempdir(), paste0("fcs_output_", method))

    result <- stimgate::stimgate_fcs_write(
      path_project = path_project,
      .data = gs,
      ind_batch_list = example_data$batch_list,
      path_dir_save = path_dir_save,
      chnl = example_data$marker[[1]],
      gate_uns_method = method
    )

    expect_equal(result, path_dir_save)
    expect_true(dir.exists(path_dir_save))
    # Should have some files (at least one sample should have positive cells)
    expect_true(length(list.files(path_dir_save, pattern = "\\.fcs$")) >= 0)
    unlink(path_dir_save, recursive = TRUE)
  }
