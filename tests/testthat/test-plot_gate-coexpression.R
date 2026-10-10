local({
  # Mock saved tables and expression so these tests need no gating run or getter.
  coex <- tibble::tibble(
    pop = c(rep("root", 6), "other"), batch = "batch1",
    ind = c(rep("1", 5), "99", "1"),
    chnlCond = c("A", "A", "B", "C", "A", "A", "A"),
    markerCond = c("a", "a", "b", "c", "a", "a", "a"),
    chnl = c("B", "B", "A", "B", "B", "B", "B"),
    marker = c("b", "b", "a", "b", "b", "b", "b"),
    gate = 10, cut = c(2, 2, 1, 3, 9, 8, 8),
    condCut = c(5, 5, 6, 7, 5, 5, 5), floor = 0, floorCond = 0,
    z = 5, purityDp = 0.9, lowered = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
  )
  gates <- tibble::tibble(
    gateName = "base", chnl = c("A", "B"), marker = c("a", "b"),
    ind = "1", batch = "batch1", gate = 10
  )
  mockPlots <- function(table = coex, error = FALSE) {
    # Give the reader a mocked getter even before that binding exists in stimgate.
    reader <- .plotCoexGates
    environment(reader) <- list2env(list(
      getStimGatesCoexpression = function(pathProject, pop = NULL) {
        if (error) stop("No coexpression table")
        table
      }
    ), parent = environment(reader))
    testthat::local_mocked_bindings(
      .plotCoexGates = reader,
      .gateGetChnl = function(...) c("A", "B"),
      stimgateMetaReadMarkerLab = function(...) c(a = "A", b = "B"),
      getStimGates = function(...) gates,
      getStimExpr = function(ind, marker = NULL, chnl = NULL, ...) {
        vars <- chnl %||% marker
        tibble::as_tibble(stats::setNames(
          lapply(vars, function(v) seq(0, 12, length.out = 30)), vars
        ))
      },
      .env = parent.frame()
    )
  }
  plots <- function(vars, showGateCyt = TRUE, ind = 1, marker = FALSE) {
    args <- list(
      ind = ind, .data = NULL, pathProject = "unused", pop = "root",
      excMin = FALSE, grid = FALSE, showGateCyt = showGateCyt
    )
    args[[if (marker) "marker" else "chnl"]] <- vars
    do.call(plotStim, args)
  }
  layerIndex <- function(p, geom) {
    which(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1)))
  }

  test_that("coexpression segments start at raised conditioning cuts on either axis", {
    skip_if_not_installed("hexbin")
    mockPlots()
    p <- plots(c("a", "b"), marker = TRUE)[[1]]
    i <- layerIndex(p, "GeomSegment")
    expect_length(i, 1L)
    seg <- ggplot2::layer_data(p, i)
    expect_equal(seg$x, c(5, 5, 1))
    expect_equal(seg$y, c(2, 2, 6))
    expect_equal(seg$xend, c(Inf, Inf, 1))
    expect_equal(seg$yend, c(2, 2, Inf))
    expect_identical(seg$colour, rep("blue", 3))
    expect_equal(seg$alpha, rep(0.6, 3))
    expect_length(unique(seg$group), 3L)

    # The drawing stage must retain the two coincident horizontal segments.
    built <- ggplot2::ggplot_build(p)
    grob <- p$layers[[i]]$geom$draw_panel(
      seg, built$layout$panel_params[[1]], p$coordinates
    )
    expect_length(grob$x0, 3L)

    reversed <- plots(c("B", "A"))[[1]]
    seg <- ggplot2::layer_data(reversed, layerIndex(reversed, "GeomSegment"))
    expect_equal(seg$x, c(2, 2, 6))
    expect_equal(seg$y, c(5, 5, 1))
    expect_equal(seg$xend, c(2, 2, Inf))
    expect_equal(seg$yend, c(Inf, Inf, 1))

    ordinary <- plots(c("A", "B"), showGateCyt = FALSE)[[1]]
    expect_length(layerIndex(ordinary, "GeomSegment"), 0L)
    expect_equal(
      ggplot2::ggplot_build(ordinary)$data,
      built$data[-i]
    )
    unstim <- plots(c("A", "B"), ind = 2)[[1]]
    expect_length(layerIndex(unstim, "GeomSegment"), 0L)
  })

  test_that("univariate lowered gates include every conditioning marker and keep overlaps", {
    mockPlots()
    p <- plots("b", marker = TRUE)[[1]]
    i <- layerIndex(p, "GeomVline")
    expect_length(i, 2L)
    lowered <- ggplot2::layer_data(p, i[[2]])
    expect_equal(lowered$xintercept, c(2, 2, 3))
    expect_identical(lowered$colour, rep("blue", 3))
    expect_equal(lowered$alpha, rep(0.6, 3))
    expect_identical(lowered$linetype, rep("dashed", 3))
    expect_length(unique(lowered$group), 3L)
    built <- ggplot2::ggplot_build(p)
    grob <- p$layers[[i[[2]]]]$geom$draw_panel(
      lowered, built$layout$panel_params[[1]], p$coordinates
    )
    expect_length(grob$x0, 3L)

    ordinary <- plots("B", showGateCyt = FALSE)[[1]]
    expect_length(layerIndex(ordinary, "GeomVline"), 1L)
    expect_equal(ggplot2::ggplot_build(ordinary)$data, built$data[-i[[2]]])
    expect_length(layerIndex(plots("B", ind = 2)[[1]], "GeomVline"), 0L)
  })

  test_that("absent, empty or unavailable coexpression tables leave plots unchanged", {
    for (table in list(NULL, coex[0, ])) {
      mockPlots(table)
      expect_equal(
        ggplot2::ggplot_build(plots("B")[[1]])$data,
        ggplot2::ggplot_build(plots("B", showGateCyt = FALSE)[[1]])$data
      )
    }
    mockPlots(error = TRUE)
    expect_null(.plotCoexGates("unused", "root"))
    expect_equal(
      ggplot2::ggplot_build(plots("B")[[1]])$data,
      ggplot2::ggplot_build(plots("B", showGateCyt = FALSE)[[1]])$data
    )
  })

  test_that("disabled cytokine gates skip the optional reader", {
    mockPlots()
    testthat::local_mocked_bindings(
      .plotCoexGates = function(...) stop("Should not read coexpression gates")
    )
    expect_s3_class(plots("B", showGateCyt = FALSE)[[1]], "ggplot")
    hidden <- plotStim(
      ind = 1, .data = NULL, pathProject = "unused", pop = "root",
      chnl = "B", excMin = FALSE, grid = FALSE, showGate = FALSE
    )[[1]]
    expect_length(layerIndex(hidden, "GeomVline"), 0L)
    expect_error(plots("B", showGateCyt = NA), "single logical")
  })
})
