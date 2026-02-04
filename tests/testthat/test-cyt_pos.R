library(testthat)

test_that("cyt_pos gates actually happen", {
  test_data <- get_test_data(
    scenario = "cyt_pos",
    dir_cache = testthat::test_path("cache", "test_data", "default"),
    clear = TRUE,
    n_ind = 1
  )
  gs <- flowWorkspace::load_gs(test_data$path_gs)
  path_project <- file.path(dirname(test_data$path_gs), "stimgate")
  invisible(stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = test_data$batch_list,
    marker = test_data$marker
  ))
  expect_true(file.exists(file.path(path_project, "gate_stats.rds")))
})
