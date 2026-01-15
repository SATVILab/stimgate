# Extracted from test-run.R:14

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "stimgate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)
path_project <- file.path(dirname(example_data$path_gs), "stimgate")
debugonce(.chnl_lab)
invisible(stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))
