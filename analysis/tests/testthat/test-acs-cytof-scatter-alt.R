scatter_alt_root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
scatter_alt_env <- new.env(parent = getNamespace("stimgate"))
for (script in c("analysis-runtime.R", "analysis-plot-style.R", "acs_cytof-plot_cyt.R", "acs_cytof-scatter-alt.R")) {
  source(file.path(scatter_alt_root, "scripts", "r", script), local = scatter_alt_env)
}

# Automated = 2 * manual + 0.1 for p1 and = manual for p4, so per-stimulation
# and pooled fits are known. Sample "f" has a failed threshold.
.acs_scatter_alt_fixture <- function() {
  tidyr::expand_grid(
    method = c("tailgate", "fbeta", "fbeta_default", "stimgate"),
    pop = c("NK cells", "CD4 T cells"),
    cyt = c("IL6", "IFNg"),
    stim = c("p4", "p1"), SampleID = letters[1:6]
  ) |>
    dplyr::mutate(
      freq_bs_man = match(.data$SampleID, letters) / 10,
      freq_bs_auto = ifelse(.data$stim == "p1", 2 * .data$freq_bs_man + 0.1, .data$freq_bs_man),
      freq_stim_man = 0.5, freq_uns_man = 0.1,
      thresholdFailed = .data$SampleID == "f"
    )
}

.acs_scatter_alt_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))

test_that("alternative scatter data keeps thesis methods, drops failed thresholds and draws StimGate last", {
  tbl <- .acs_scatter_alt_fixture()
  tbl$freq_bs_auto[1] <- NA_real_
  out <- scatter_alt_env$.acsCytofScatterAltData(tbl, c("stimgate", "fbeta_default", "tailgate"))
  expect_equal(levels(out$method), c("stimgate", "fbeta_default", "tailgate"))
  expect_false(any(out$SampleID == "f"))
  expect_equal(nrow(out), 3 * 2 * 2 * 2 * 5 - 1)
  expect_equal(levels(out$pop), c("CD4 T cells", "NK cells"))
  expect_equal(levels(out$cyt)[1:2], c("IFNg", "IL2"))
  expect_equal(levels(out$stim), c("mtb", "p1", "p4", "ebv"))
  expect_equal(unique(as.character(out$method)), c("tailgate", "fbeta_default", "stimgate"))
})

test_that("signed square-root axes keep signs and use round breaks in original units", {
  trans <- scatter_alt_env$.acsCytofSignedSqrtTrans()
  x <- c(-4, -0.01, 0, 0.25, 9)
  expect_equal(trans$transform(x), c(-2, -0.1, 0, 0.5, 3))
  expect_equal(trans$inverse(trans$transform(x)), x)
  brk <- scatter_alt_env$.acsCytofSignedSqrtBreaks
  expect_equal(brk(c(0, 15)), c(0, 2, 10))
  expect_equal(brk(c(0, 0.6)), c(0, 0.1, 0.6))
  expect_equal(brk(c(-2, 8)), c(-2, 0, 2, 8))
  expect_equal(brk(c(-0.001, 8)), c(0, 2, 8))
})

test_that("each panel has one range for both axes, including zero", {
  tbl <- .acs_scatter_alt_fixture()
  out <- scatter_alt_env$.acsCytofScatterAltData(tbl, "stimgate")
  lim <- scatter_alt_env$.acsCytofScatterAltLimits(out)
  expect_equal(nrow(lim), 4)
  expect_equal(lim$lo, rep(0, 4))
  expect_equal(lim$hi, rep(1.1, 4))
})

test_that("lines of best fit are least-squares lines on the original scale", {
  tbl <- .acs_scatter_alt_fixture()
  out <- scatter_alt_env$.acsCytofScatterAltData(tbl, "stimgate")
  by_stim <- scatter_alt_env$.acsCytofScatterAltFits(out, c("pop", "cyt", "stim"), n = 50)
  expect_equal(nrow(by_stim), 2 * 2 * 2 * 50)
  p1 <- by_stim[by_stim$stim == "p1", ]
  expect_equal(p1$slope, rep(2, nrow(p1)))
  expect_equal(p1$intercept, rep(0.1, nrow(p1)), tolerance = 1e-12)
  expect_equal(p1$freq_bs_auto, 2 * p1$freq_bs_man + 0.1, tolerance = 1e-12)
  expect_equal(range(by_stim$freq_bs_man), c(0.1, 0.5))
  pooled <- scatter_alt_env$.acsCytofScatterAltFits(out, c("pop", "cyt"))
  expected <- stats::coef(stats::lm(freq_bs_auto ~ freq_bs_man,
    data = out[out$pop == "NK cells" & out$cyt == "IL6", ]))
  nk <- pooled[pooled$pop == "NK cells" & pooled$cyt == "IL6", ]
  expect_equal(unique(nk$slope), unname(expected[[2]]))
  # Lines stop at the panel's range.
  lim <- scatter_alt_env$.acsCytofScatterAltLimits(out)
  lim$hi <- 0.6
  cut <- scatter_alt_env$.acsCytofScatterAltFits(out, c("pop", "cyt", "stim"), limits = lim)
  expect_true(all(cut$freq_bs_auto <= 0.6))
  expect_lt(nrow(cut), nrow(scatter_alt_env$.acsCytofScatterAltFits(out, c("pop", "cyt", "stim"))))
})

test_that("groups with too few points or too little spread get no line", {
  tbl <- .acs_scatter_alt_fixture()
  out <- scatter_alt_env$.acsCytofScatterAltData(tbl, "stimgate")
  few <- out[out$SampleID %in% c("a", "b", "c"), ]
  expect_equal(nrow(scatter_alt_env$.acsCytofScatterAltFits(few, c("pop", "cyt", "stim"))), 0)
  narrow <- out
  narrow$freq_bs_man[narrow$SampleID %in% c("c", "d", "e")] <- 0.2
  expect_equal(nrow(scatter_alt_env$.acsCytofScatterAltFits(narrow, c("pop", "cyt", "stim"))), 0)
})

test_that("all-method figures have one square row per population, coloured by method", {
  tbl <- .acs_scatter_alt_fixture()
  for (fit in c("stim", "panel", "none")) {
    res <- scatter_alt_env$.acsCytofScatterAlt(tbl, fit = fit)
    expect_s3_class(res$plot, "ggplot")
    expect_equal(names(res$rows), c("CD4 T cells", "NK cells"))
    row <- res$rows[["CD4 T cells"]]
    geoms <- .acs_scatter_alt_geoms(row)
    built <- ggplot2::ggplot_build(row)
    points <- built$data[[which(geoms == "GeomPoint")]]
    expect_setequal(unique(points$colour), unname(scatter_alt_env$.analysis_method_colours[c("stimgate", "fbeta", "tailgate")]))
    expect_true(all(points$shape %in% c(22, 23)))
    expect_equal(nrow(points), 3 * 2 * 2 * 5)
    expect_true("GeomAbline" %in% geoms)
    expect_equal(row$theme$aspect.ratio, 1)
    # Every cytokine keeps its column, so rows line up.
    expect_equal(nrow(built$layout$layout), 6)
    for (panel in built$layout$panel_params) {
      expect_equal(panel$x$continuous_range, panel$y$continuous_range)
    }
    expect_equal(row$scales$get_scales("x")$get_transformation()$name, "signed_sqrt")
    expect_equal(sum(geoms == "GeomLine"), as.integer(fit != "none"))
    if (fit != "none") {
      lines <- built$data[[which(geoms == "GeomLine")]]
      expect_equal(length(unique(lines$group)), 3 * 2 * if (fit == "stim") 2 else 1)
      expect_equal(length(unique(lines$linetype)), if (fit == "stim") 2 else 1)
    }
    expect_equal(res$width, scatter_alt_env$.analysis_thesis_width)
  }
})

test_that("one-method figures colour by stimulation and validate the method", {
  tbl <- .acs_scatter_alt_fixture()
  res <- scatter_alt_env$.acsCytofScatterAlt(tbl, method = "fbeta_default", fit = "panel")
  expect_true(all(as.character(res$data$method) == "fbeta_default"))
  row <- res$rows[["NK cells"]]
  geoms <- .acs_scatter_alt_geoms(row)
  built <- ggplot2::ggplot_build(row)
  points <- built$data[[which(geoms == "GeomPoint")]]
  expect_setequal(unique(points$colour), unname(scatter_alt_env$.acsCytofValidationStimColours()[c("p1", "p4")]))
  lines <- built$data[[which(geoms == "GeomLine")]]
  expect_equal(unique(lines$colour), "black")
  expect_equal(length(unique(lines$group)), 2)
  expect_error(scatter_alt_env$.acsCytofScatterAlt(tbl, method = "fbeta"), "method must")
  expect_error(
    scatter_alt_env$.acsCytofScatterAlt(dplyr::mutate(tbl, thresholdFailed = TRUE), method = "stimgate"),
    "No matched"
  )
})
