library(testthat)

test_that("stimgate_plot function exists", {
  # Just test that the function exists and is callable
  expect_true(exists("stimgate_plot", envir = asNamespace("stimgate")))
  expect_true(is.function(stimgate_plot))
})

test_that("stimgate_gate runs", {
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
  p <- stimgate_plot(
    ind = example_data$batch_list[[1]], # indices in `gs` to plot
    .data = gs, # GatingSet
    path_project = path_project,
    marker = example_data$marker,
    grid = TRUE
  )
  expect_true(inherits(p, "ggplot"))
})

# Test edge cases for plot_gate functionality

test_that("stimgate_plot returns NULL when p_list is empty", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")

  # Test with empty ind list to generate empty p_list
  result <- stimgate_plot(
    ind = list(),
    .data = gs,
    path_project = path_project,
    marker = example_data$marker,
    grid = TRUE
  )
  expect_null(result)
})

test_that(".plot_gate_bv returns NULL for single marker", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")

  # Test single marker scenario
  single_marker <- example_data$marker[1]
  result <- stimgate:::.plot_gate_bv(
    marker = single_marker,
    ind = example_data$batch_list[[1]],
    ind_lab = NULL,
    .data = gs,
    marker_lab = NULL,
    path_project = path_project,
    exc_min = TRUE,
    limits_expand = NULL,
    limits_equal = FALSE,
    show_gate = TRUE,
    min_cell = 10
  )
  expect_null(result)
})

test_that("hexbin installation check works", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")

  # First run gating to create necessary gate data
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))

  # Test with two markers to trigger hexbin requirement
  # The function should handle hexbin installation automatically
  expect_no_error({
    result <- stimgate:::.plot_gate_bv(
      marker = example_data$marker,
      ind = example_data$batch_list[[1]],
      ind_lab = NULL,
      .data = gs,
      marker_lab = NULL,
      path_project = path_project,
      exc_min = TRUE,
      limits_expand = NULL,
      limits_equal = FALSE,
      show_gate = TRUE,
      min_cell = 10
    )
  })
})

test_that("plot functions handle min_cell threshold correctly", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")

  # First run gating to create necessary gate data
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))

  # Test with very high min_cell to trigger early return
  result <- stimgate:::.plot_gate_bv(
    marker = example_data$marker,
    ind = example_data$batch_list[[1]],
    ind_lab = NULL,
    .data = gs,
    marker_lab = NULL,
    path_project = path_project,
    exc_min = TRUE,
    limits_expand = NULL,
    limits_equal = FALSE,
    show_gate = TRUE,
    min_cell = 999999  # Very high threshold
  )
  # Should return a list with NULLs filtered out, or NULL
  expect_true(is.null(result) || (is.list(result) && length(result) == 0))
})

test_that(".plot_get_lab handles various val_lab configurations", {
  # Test with NULL val_lab
  result1 <- stimgate:::.plot_get_lab(
    val = c("A", "B"),
    val_lab = NULL,
    i = NULL
  )
  expect_equal(result1, c("A", "B"))

  # Test with named val_lab
  result2 <- stimgate:::.plot_get_lab(
    val = c("A", "B"),
    val_lab = c("A" = "Label A", "B" = "Label B"),
    i = NULL
  )
  expect_equal(result2, c("Label A", "Label B"))
  expect_null(names(result2))

  # Test with unnamed val_lab and i NULL
  result3 <- stimgate:::.plot_get_lab(
    val = c("A", "B"),
    val_lab = c("Label A", "Label B"),
    i = NULL
  )
  expect_equal(result3, c("Label A", "Label B"))
  expect_null(names(result3))

  # Test with unnamed val_lab and i non-null
  result4 <- stimgate:::.plot_get_lab(
    val = c("A", "B"),
    val_lab = c("Label A", "Label B"),
    i = 1
  )
  expect_equal(result4, "Label A")
  expect_null(names(result4))
})

test_that(".plot_gate_uv returns NULL when all markers return NULL", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")

  # First run gating to create necessary gate data
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))

  # Test with empty ind to generate NULL results
  result <- stimgate:::.plot_gate_uv(
    ind = list(),
    ind_lab = NULL,
    .data = gs,
    marker = example_data$marker,
    exc_min = TRUE,
    marker_lab = NULL,
    show_gate = TRUE,
    path_project = path_project,
    min_cell = 10
  )
  expect_null(result)
})

test_that(".plot_gate_uv_marker returns NULL when plot_tbl is NULL", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")

  # First run gating to create necessary gate data
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))

  # Test with empty ind list to generate NULL plot_tbl
  result <- stimgate:::.plot_gate_uv_marker(
    marker = example_data$marker[1],
    ind = list(),
    .data = gs,
    exc_min = TRUE,
    ind_lab = NULL,
    marker_lab = NULL,
    show_gate = TRUE,
    path_project = path_project,
    min_cell = 10
  )
  expect_null(result)
})

test_that(".plot_gate_uv_marker_get_plot_tbl returns NULL for insufficient cells", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)

  # Test with very high min_cell
  result <- stimgate:::.plot_gate_uv_marker_get_plot_tbl(
    ind = example_data$batch_list[[1]],
    .data = gs,
    marker = example_data$marker[1],
    exc_min = TRUE,
    ind_lab = NULL,
    min_cell = 999999  # Very high threshold
  )
  expect_null(result)
})

test_that(".plot_gate_uv_marker_plot_init handles different condition branches", {
  # Create mock plot_tbl for testing
  plot_tbl <- data.frame(
    x = 1:10,
    y = 1:10,
    type = rep(c("raw", "adj"), 5),
    ind_lab = rep(c("Sample1", "Sample2"), 5)
  )

  # Test exc_min = TRUE, multiple ind
  p1 <- stimgate:::.plot_gate_uv_marker_plot_init(
    plot_tbl = plot_tbl,
    exc_min = TRUE,
    ind = c(1, 2),
    ind_lab = c("Sample1", "Sample2")
  )
  expect_s3_class(p1, "ggplot")

  # Test exc_min = TRUE, single ind
  p2 <- stimgate:::.plot_gate_uv_marker_plot_init(
    plot_tbl = plot_tbl,
    exc_min = TRUE,
    ind = c(1),
    ind_lab = c("Sample1")
  )
  expect_s3_class(p2, "ggplot")

  # Test exc_min = FALSE, multiple ind
  p3 <- stimgate:::.plot_gate_uv_marker_plot_init(
    plot_tbl = plot_tbl,
    exc_min = FALSE,
    ind = c(1, 2),
    ind_lab = c("Sample1", "Sample2")
  )
  expect_s3_class(p3, "ggplot")

  # Test exc_min = FALSE, single ind
  p4 <- stimgate:::.plot_gate_uv_marker_plot_init(
    plot_tbl = plot_tbl,
    exc_min = FALSE,
    ind = c(1),
    ind_lab = c("Sample1")
  )
  expect_s3_class(p4, "ggplot")
})

test_that(".plot_grid returns p_list when plot = FALSE", {
  # Create mock plot list
  p_list <- list(
    plot1 = ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 1, y = 1)),
    plot2 = ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 2, y = 2))
  )

  # Test with plot = FALSE
  result <- stimgate:::.plot_grid(
    plot = FALSE,
    p_list = p_list,
    n_col = 2
  )
  expect_identical(result, p_list)

  # Test with plot = TRUE (should return combined plot)
  result2 <- stimgate:::.plot_grid(
    plot = TRUE,
    p_list = p_list,
    n_col = 2
  )
  expect_s3_class(result2, "ggplot")
})

# Additional comprehensive edge case tests
test_that("comprehensive edge case coverage for plot_gate functions", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")

  # Run gating once for multiple tests
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))

  # Test .plot_gate_uv_marker_get_plot_tbl_ind with insufficient cells
  expect_null(stimgate:::.plot_gate_uv_marker_get_plot_tbl_ind(
    ind = c(1),  # Single index
    .data = gs,
    marker = example_data$marker[1],
    exc_min = TRUE,
    min_cell = 999999  # Impossible threshold
  ))

  # Test .plot_get_ex_tbl with single index
  ex_tbl <- stimgate:::.plot_get_ex_tbl(
    ind = c(1),
    .data = gs,
    marker = example_data$marker[1],
    exc_min = TRUE
  )
  expect_true(is.data.frame(ex_tbl))
  expect_true(nrow(ex_tbl) > 0)

  # Test .plot_add_axis_title with single and multiple markers
  p_base <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 1, y = 1))
  p_single <- stimgate:::.plot_add_axis_title(p_base, example_data$marker[1], NULL)
  expect_s3_class(p_single, "ggplot")

  p_double <- stimgate:::.plot_add_axis_title(p_base, example_data$marker, NULL)
  expect_s3_class(p_double, "ggplot")

  # Test .plot_add_title
  p_titled <- stimgate:::.plot_add_title(p_base, c(1, 2), 1, NULL)
  expect_s3_class(p_titled, "ggplot")

  # Test .plot_add_gate when show_gate = FALSE
  p_no_gate <- stimgate:::.plot_add_gate(
    p_base, gs, c(1), example_data$marker[1], path_project, show_gate = FALSE
  )
  expect_identical(p_no_gate, p_base)
})

test_that("test plot_cyto import and dependencies", {
  # Test that plot_cyto function is accessible
  expect_true(exists("plot_cyto", envir = asNamespace("stimgate")))

  # Test hexbin namespace checking
  hexbin_available <- requireNamespace("hexbin", quietly = TRUE)
  expect_true(is.logical(hexbin_available))
})
