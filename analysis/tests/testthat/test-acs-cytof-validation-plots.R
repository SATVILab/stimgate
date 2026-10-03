root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)
script_plot <- file.path(root_dir, "scripts", "r", "acs_cytof-plot_cyt.R")

env <- new.env(parent = getNamespace("stimgate"))
source(script_plot, local = env)

.acs_validation_fixture <- function() {
  tidyr::expand_grid(
    method = c("stimgate", "fbeta", "tailgate"),
    pop = c("CD4 T cells", "B cells"),
    cyt = c("IFNg", "IL2"),
    stim = c("p1", "mtb", "ebv", "p4"),
    SampleID = paste0("sample", 1:4)
  ) |>
    dplyr::group_by(method, pop, cyt, stim) |>
    dplyr::mutate(
      freq_bs_man = seq(0.1, 0.4, length.out = dplyr::n()),
      freq_bs_auto = .data$freq_bs_man * 1.1,
      freq_stim_man = .data$freq_bs_man + 0.02,
      freq_uns_man = 0.02
    ) |>
    dplyr::ungroup()
}

test_that("CCC is one for identical vectors", {
  expect_equal(env$.acsCytofValidationCcc(1:5, 1:5), 1)
  expect_true(is.na(env$.acsCytofValidationCcc(rep(1, 5), 1:5)))
})

test_that("correlation table contains PCC and CCC by method and stratum", {
  comparison_tbl <- .acs_validation_fixture()
  correlation_tbl <- env$.acsCytofValidationCorrelationTable(comparison_tbl)

  expect_true(all(c("method", "pop", "cyt", "stim", "pcc", "ccc") %in%
    names(correlation_tbl)))
  expect_equal(nrow(correlation_tbl), 3 * 2 * 2 * 4)
  expect_true(all(abs(correlation_tbl$pcc - 1) < 1e-10))
})

test_that("correlation table keeps only manuscript-eligible signal strata", {
  comparison_tbl <- tibble::tibble(
    method = "stimgate",
    pop = "CD4 T cells",
    cyt = "IFNg",
    stim = rep(c("p1", "mtb", "ebv"), each = 4L),
    SampleID = rep(paste0("sample", 1:4), 3L),
    freq_bs_man = c(
      0.10, 0.20, 0.30, 0.40,
      0.10, 0.20, 0.30, 0.40,
      0.001, 0.002, 0.003, 0.004
    ),
    freq_bs_auto = c(
      0.11, 0.22, 0.33, 0.44,
      0.11, 0.22, 0.33, 0.44,
      0.0011, 0.0022, 0.0033, 0.0044
    ),
    freq_stim_man = c(
      0.12, 0.22, 0.32, 0.42,
      0.03, 0.04, 0.05, 0.06,
      0.12, 0.22, 0.32, 0.42
    ),
    freq_uns_man = 0.02
  )

  correlation_tbl <- env$.acsCytofValidationCorrelationTable(comparison_tbl)

  expect_equal(correlation_tbl$stim, "p1")
  expect_equal(correlation_tbl$n, 4L)
  expect_equal(correlation_tbl$pcc, 1, tolerance = 1e-12)
})


test_that("validation plotting helpers return ggplot objects", {
  comparison_tbl <- .acs_validation_fixture()
  correlation_tbl <- env$.acsCytofValidationCorrelationTable(comparison_tbl)

  expect_s3_class(
    env$.acsCytofValidationPlotScatter(comparison_tbl, "stimgate"),
    "ggplot"
  )
  expect_s3_class(
    env$.acsCytofValidationPlotCorrelation(
      correlation_tbl,
      method = "stimgate",
      metric = "ccc",
      realPopulationsOnly = TRUE
    ),
    "ggplot"
  )
})
