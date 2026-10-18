scatter_root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
scatter_env <- new.env(parent = getNamespace("stimgate"))
for (script in c("analysis-runtime.R", "analysis-plot-style.R", "acs_cytof-plot_cyt.R", "acs_cytof-scatter-stacked.R")) {
  source(file.path(scatter_root, "scripts", "r", script), local = scatter_env)
}

.acs_stacked_scatter_fixture <- function() {
  tidyr::expand_grid(
    method = c("tailgate", "fbeta", "fbeta_default", "stimgate"),
    pop = c("NK cells", "CD4 T cells"),
    cyt = c("IL6", "IFNg", "IL2", "TNF", "IL17", "IL22"),
    stim = c("p1", "p4"), SampleID = c("a", "b")
  ) |>
    dplyr::mutate(
      freq_bs_man = -0.005, freq_bs_auto = -0.004,
      freq_stim_man = 0.001, freq_uns_man = 0.006,
      thresholdFailed = .data$SampleID == "b"
    )
}

test_that("stacked scatter retains raw scatter cohorts and masks failed thresholds", {
  tbl <- .acs_stacked_scatter_fixture()
  plot <- scatter_env$.acsCytofStackedScatterBlock(tbl, "fbeta_default", "CD4 T cells", 6)
  expected <- tbl |> dplyr::filter(.data$method == "fbeta_default", .data$pop == "CD4 T cells")
  expect_equal(nrow(plot$data), nrow(expected))
  expect_equal(as.character(plot$data$stim), expected$stim)
  expect_equal(plot$data$freq_bs_man, expected$freq_bs_man)
  expect_equal(plot$data$freq_bs_auto[!expected$thresholdFailed], expected$freq_bs_auto[!expected$thresholdFailed])
  expect_true(all(is.na(plot$data$freq_bs_auto[expected$thresholdFailed])))
  expect_equal(levels(plot$data$cyt), c("IFNg", "IL2", "TNF", "IL17", "IL22", "IL6"))
  expect_equal(plot$facet$params$ncol, 6)
  expect_false(plot$facet$params$drop)
})

test_that("stacked layouts order thesis methods and populations with readable dimensions", {
  tbl <- .acs_stacked_scatter_fixture()
  for (columns in c(3L, 6L)) {
    result <- scatter_env$.acsCytofStackedScatter(tbl, population = "CD4 T cells", ncol = columns)
    expect_equal(result$groups, c("stimgate", "fbeta_default", "tailgate"))
    expect_length(result$blocks, 3)
    expect_true(all(vapply(result$blocks, function(p) p$facet$params$ncol == columns, logical(1))))
    expect_gte(result$width, columns * 60)
  }
  result <- scatter_env$.acsCytofStackedScatter(tbl, method = "fbeta_default")
  expect_equal(result$groups, c("CD4 T cells", "NK cells"))
  expect_length(result$blocks, 2)
  expect_error(scatter_env$.acsCytofStackedScatter(tbl), "exactly one")
  expect_error(scatter_env$.acsCytofStackedScatter(tbl, method = "fbeta"), "method must")
})

test_that("a missing method remains an explicitly labelled block", {
  tbl <- .acs_stacked_scatter_fixture() |> dplyr::filter(.data$method != "tailgate")
  result <- scatter_env$.acsCytofStackedScatter(tbl, population = "CD4 T cells")
  expect_length(result$blocks, 3)
  expect_equal(result$blocks[[3]]$labels$subtitle, "No matched rows available")
})
