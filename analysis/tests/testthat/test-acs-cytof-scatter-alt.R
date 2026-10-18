scatter_alt_root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
scatter_alt_env <- new.env(parent = getNamespace("stimgate"))
for (script in c("analysis-runtime.R", "analysis-plot-style.R", "acs_cytof-plot_cyt.R", "acs_cytof-scatter-alt.R")) {
  source(file.path(scatter_alt_root, "scripts", "r", script), local = scatter_alt_env)
}

# Automated = 2 * manual + 0.1 for p1 and = manual for p4, so per-stimulation
# and pooled fits are known. Sample "c" has a failed threshold.
.acs_scatter_alt_fixture <- function() {
  tidyr::expand_grid(
    method = c("tailgate", "fbeta", "fbeta_default", "stimgate"),
    pop = c("NK cells", "CD4 T cells"),
    cyt = c("IL6", "IFNg"),
    stim = c("p4", "p1"), SampleID = c("a", "b", "c")
  ) |>
    dplyr::mutate(
      freq_bs_man = match(.data$SampleID, c("a", "b", "c")) / 10,
      freq_bs_auto = ifelse(.data$stim == "p1", 2 * .data$freq_bs_man + 0.1, .data$freq_bs_man),
      freq_stim_man = 0.5, freq_uns_man = 0.1,
      thresholdFailed = .data$SampleID == "c"
    )
}

test_that("alternative scatter data keeps thesis methods and drops failed thresholds", {
  tbl <- .acs_scatter_alt_fixture()
  tbl$freq_bs_auto[1] <- NA_real_
  out <- scatter_alt_env$.acsCytofScatterAltData(tbl, c("stimgate", "fbeta_default", "tailgate"))
  expect_equal(levels(out$method), c("stimgate", "fbeta_default", "tailgate"))
  expect_false(any(out$SampleID == "c"))
  expect_equal(nrow(out), 3 * 2 * 2 * 2 * 2 - 1)
  expect_equal(levels(out$pop), c("CD4 T cells", "NK cells"))
  expect_equal(levels(out$cyt)[1:2], c("IFNg", "IL2"))
  expect_equal(levels(out$stim), c("mtb", "p1", "p4", "ebv"))
})

test_that("lines of best fit are least-squares lines over each group's manual range", {
  tbl <- .acs_scatter_alt_fixture()
  out <- scatter_alt_env$.acsCytofScatterAltData(tbl, "stimgate")
  by_stim <- scatter_alt_env$.acsCytofScatterAltFits(out, c("pop", "cyt", "stim"))
  expect_equal(nrow(by_stim), 2 * 2 * 2 * 2)
  p1 <- by_stim[by_stim$stim == "p1", ]
  expect_equal(p1$slope, rep(2, nrow(p1)))
  expect_equal(p1$intercept, rep(0.1, nrow(p1)), tolerance = 1e-12)
  expect_equal(range(by_stim$freq_bs_man), c(0.1, 0.2))
  pooled <- scatter_alt_env$.acsCytofScatterAltFits(out, c("pop", "cyt"))
  expected <- stats::coef(stats::lm(freq_bs_auto ~ freq_bs_man,
    data = out[out$pop == "NK cells" & out$cyt == "IL6", ]))
  nk <- pooled[pooled$pop == "NK cells" & pooled$cyt == "IL6", ]
  expect_equal(unique(nk$slope), unname(expected[[2]]))
  # No spread in manual frequency: no line.
  flat <- out
  flat$freq_bs_man <- 0
  expect_equal(nrow(scatter_alt_env$.acsCytofScatterAltFits(flat, c("pop", "cyt"))), 0)
})

test_that("all-method figures colour by method and draw the requested lines", {
  tbl <- .acs_scatter_alt_fixture()
  layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
  for (fit in c("stim", "panel", "none")) {
    res <- scatter_alt_env$.acsCytofScatterAlt(tbl, fit = fit)
    built <- ggplot2::ggplot_build(res$plot)
    points <- built$data[[which(layer_geoms(res$plot) == "GeomPoint")]]
    expect_setequal(unique(points$colour), unname(scatter_alt_env$.analysis_method_colours[c("stimgate", "fbeta", "tailgate")]))
    expect_true(all(points$shape %in% c(22, 23)))
    expect_true(all(points$size >= 1.5))
    expect_equal(nrow(points), 3 * 2 * 2 * 2 * 2)
    expect_true("GeomAbline" %in% layer_geoms(res$plot))
    n_lines <- sum(layer_geoms(res$plot) == "GeomLine")
    expect_equal(n_lines, as.integer(fit != "none"))
    if (fit != "none") {
      lines <- built$data[[which(layer_geoms(res$plot) == "GeomLine")]]
      expect_equal(length(unique(lines$group)), 3 * 2 * 2 * if (fit == "stim") 2 else 1)
      expect_equal(length(unique(lines$linetype)), if (fit == "stim") 2 else 1)
    }
    expect_equal(nrow(ggplot2::ggplot_build(res$plot)$layout$layout), 4)
    expect_equal(res$width, scatter_alt_env$.analysis_thesis_width)
  }
})

test_that("one-method figures colour by stimulation and validate the method", {
  tbl <- .acs_scatter_alt_fixture()
  res <- scatter_alt_env$.acsCytofScatterAlt(tbl, method = "fbeta_default", fit = "panel")
  expect_true(all(as.character(res$data$method) == "fbeta_default"))
  built <- ggplot2::ggplot_build(res$plot)
  geoms <- vapply(res$plot$layers, function(l) class(l$geom)[1], character(1))
  points <- built$data[[which(geoms == "GeomPoint")]]
  expect_setequal(unique(points$colour), unname(scatter_alt_env$.acsCytofValidationStimColours()[c("p1", "p4")]))
  lines <- built$data[[which(geoms == "GeomLine")]]
  expect_equal(unique(lines$colour), "black")
  expect_equal(length(unique(lines$group)), 4)
  expect_error(scatter_alt_env$.acsCytofScatterAlt(tbl, method = "fbeta"), "method must")
  expect_error(
    scatter_alt_env$.acsCytofScatterAlt(dplyr::mutate(tbl, thresholdFailed = TRUE), method = "stimgate"),
    "No matched"
  )
})
