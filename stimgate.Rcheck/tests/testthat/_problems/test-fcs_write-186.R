# Extracted from test-fcs_write.R:186

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
path_dir_save <- file.path(tempdir(), "fcs_output_validation")
stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]]
  )
