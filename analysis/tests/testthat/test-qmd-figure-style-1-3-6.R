.style_qmd_files <- c(
  "1-sim-trans.qmd", "3-sim-bw-est-base.qmd", "4-sim-bw-est-norm.qmd",
  "5-sim-bw-est-adaptive.qmd", "6-sim-bw-freq_bs-adaptive.qmd"
)

.style_qmd_lines <- function(file) {
  readLines(file.path(
    testthat::test_path(), "..", "..", "..", "analysis", file
  ), warn = FALSE)
}

test_that("figure QMDs 1 and 3-6 source the shared style and use it for plots", {
  for (file in .style_qmd_files) {
    lines <- .style_qmd_lines(file)
    code <- lines[!grepl("^\\s*#", lines)]
    expect_true(any(grepl("analysis-plot-style.R", lines, fixed = TRUE)), info = file)
    expect_true(any(grepl(".analysis_theme(", code, fixed = TRUE)), info = file)
    expect_true(any(grepl("\\.analysis_(save|print_save)_fig\\(", code)), info = file)
    expect_false(
      any(grepl("theme_cowplot|ggsave\\(|labs\\(title|title = paste0|ggtitle", code)),
      info = file
    )
    expect_false(any(grepl("c(gamma = \"Gamma\"", code, fixed = TRUE)), info = file)
  }
})

test_that("figure loops in QMDs 3-6 write headings and avoid fig- labels", {
  for (file in .style_qmd_files[-1]) {
    lines <- .style_qmd_lines(file)
    expect_false(any(grepl("^#\\| label: fig-", lines)), info = file)
    expect_true(any(grepl("^#\\| results: asis", lines)), info = file)
    expect_true(any(grepl(".analysis_heading(", lines, fixed = TRUE)), info = file)
    expect_true(any(grepl(".analysis_print_save_fig(", lines, fixed = TRUE)), info = file)
    expect_false(any(grepl("^\\s*print\\(p\\)", lines)), info = file)
  }
})

test_that("transformation labels in sim-misc follow the standard order", {
  env <- new.env(parent = globalenv())
  source(
    file.path(testthat::test_path(), "../../../scripts/r/analysis-plot-style.R"),
    local = env
  )
  source(
    file.path(testthat::test_path(), "../../../scripts/r/sim-misc.R"),
    local = env
  )
  expect_identical(env$.simMiscGetTransPretty(), env$.analysis_trans_labels)
})
