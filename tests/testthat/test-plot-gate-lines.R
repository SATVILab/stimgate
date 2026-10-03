test_that("vectorised gate lines preserve axes, overlaps, styles and limits", {
  gates <- tibble::tibble(
    gateName = c("base", "base", "extra", "base", "base", "base"),
    chnl = c("A", "A", "A", "B", "B", "A"),
    marker = c("a", "a", "a", "b", "b", "a"),
    ind = c("1", "2", "1", "1", "2", "99"),
    batch = "b", gate = c(2, 2, -4, 3, 5, 99)
  )
  testthat::local_mocked_bindings(
    .gateGetChnl = function(...) c("A", "B"),
    getStimGates = function(...) gates
  )
  base <- ggplot2::ggplot(
    data.frame(x = 0:1, y = 0:1), ggplot2::aes(x = x, y = y)
  ) + ggplot2::geom_point()
  plot <- .plotAddGate(
    base, ind = c("1", "2"), marker = NULL, chnl = c("A", "B"),
    pop = "root", pathProject = "unused", showGate = TRUE
  )
  expect_length(plot$layers, 5L)
  built <- ggplot2::ggplot_build(plot)
  vertical <- built$data[[2L]]
  horizontal <- built$data[[4L]]
  expect_equal(vertical$xintercept, c(2, 2, -4))
  expect_equal(horizontal$yintercept, c(3, 5))
  expect_identical(vertical$colour, rep("red", 3L))
  expect_identical(vertical$alpha, rep(0.5, 3L))
  expect_identical(horizontal$colour, rep("red", 2L))
  expect_identical(horizontal$alpha, rep(0.5, 2L))

  # Verify that the drawing stage retains both lines at threshold 2.
  grob <- plot$layers[[2L]]$geom$draw_panel(
    vertical, built$layout$panel_params[[1L]], plot$coordinates
  )
  expect_length(grob$x0, 3L)

  reference <- ggplot2::ggplot_build(base + ggplot2::expand_limits(
    x = c(2.2, -4.4), y = c(3.3, 5.5)
  ))
  expect_equal(
    built$layout$panel_params[[1L]]$x.range,
    reference$layout$panel_params[[1L]]$x.range
  )
  expect_equal(
    built$layout$panel_params[[1L]]$y.range,
    reference$layout$panel_params[[1L]]$y.range
  )
  expect_identical(.plotAddGate(
    base, ind = "missing", marker = NULL, chnl = c("A", "B"),
    pop = "root", pathProject = "unused", showGate = TRUE
  ), base)
})
