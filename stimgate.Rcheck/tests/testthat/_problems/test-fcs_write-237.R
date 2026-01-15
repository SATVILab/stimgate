# Extracted from test-fcs_write.R:237

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
gate_tbl <- data.frame(
    chnl = rep(example_data$marker[[1]], length(unlist(example_data$batch_list))),
    marker = rep("BC1", length(unlist(example_data$batch_list))),
    batch = paste0("batch_", rep(seq_along(example_data$batch_list),
                                times = sapply(example_data$batch_list, length))),
    ind = as.character(unlist(example_data$batch_list)),
    gate = rep(0.5, length(unlist(example_data$batch_list))),
    gate_cyt = rep(0.5, length(unlist(example_data$batch_list))),
    gate_single = rep(0.5, length(unlist(example_data$batch_list))),
    stringsAsFactors = FALSE
  )
path_dir_save <- file.path(tempdir(), "fcs_output_custom_gate")
result <- stimgate::stimgate_fcs_write(
    path_project = tempdir(), # Not used when gate_tbl provided
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]],
    gate_tbl = gate_tbl
  )
