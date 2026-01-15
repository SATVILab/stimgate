# Extracted from test-gate-get.R:29

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "stimgate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
skip_if_not_installed("stimgate")
library(stimgate)
example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)
path_project <- file.path(dirname(example_data$path_gs), "stimgate_gate_get_test")
invisible(stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))
gate_tbl <- stimgate_gate_get(path_project)
expect_true(is.data.frame(gate_tbl))
expect_true(nrow(gate_tbl) > 0)
