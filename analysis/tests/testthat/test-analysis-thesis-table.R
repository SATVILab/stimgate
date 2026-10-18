root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

.thesisTableEnv <- function() {
  env <- new.env(parent = baseenv())
  source(file.path(root_dir, "scripts", "r", "analysis-thesis-table.R"), local = env)
  env
}

test_that("thesis tables are captioned pipe tables inside scoped Typst rules", {
  env <- .thesisTableEnv()
  cells <- cbind(c("A", "", "B", "", ""), letters[1:5], sprintf("%.3f", 1:5 / 10))
  bold <- matrix(FALSE, 5, 3)
  bold[2, 2] <- TRUE
  md <- env$.analysis_md_table(
    cells, header = c("Group", "Setting", "Value"),
    align = c("left", "left", "right"), bold = bold, groups = c(2, 3),
    caption = c("First sentence.", "Second sentence."), label = "tbl-test"
  )
  lines <- strsplit(md, "\n", fixed = TRUE)[[1]]
  expect_identical(lines[[1]], "```{=typst}")
  expect_identical(utils::tail(lines, 3), c("```{=typst}", "]", "```"))
  expect_true("    bottom: if y == 0 or y == 5 {" %in% lines)
  expect_true("    } else if y == 2 {" %in% lines)
  rows <- grep("^\\| ", lines, value = TRUE)
  expect_length(rows, 7L)
  expect_true(all(lengths(regmatches(rows, gregexpr("|", rows, fixed = TRUE))) == 4L))
  expect_identical(rows[[2]], "| :--- | :--- | ---: |")
  expect_identical(rows[[4]], "|  | **b** | 0.200 |")
  expect_true(": First sentence." %in% lines)
  expect_true("Second sentence. {#tbl-test}" %in% lines)
  expect_error(env$.analysis_md_table(cells, c("a", "b"), "x", "tbl-x"), "do not match")
  expect_error(env$.analysis_md_table(cells, c("a", "b", "c"), "x", "tbl-x", groups = 4), "do not match")
})

test_that("one group has no separator, and numbers keep their precision", {
  env <- .thesisTableEnv()
  md <- env$.analysis_md_table(matrix("x", 3, 1), "H", "Cap.", "tbl-one")
  expect_false(grepl("D6D6D6", md, fixed = TRUE))
  expect_match(md, "bottom: if y == 0 or y == 3 {", fixed = TRUE)
  expect_identical(env$.analysis_md_num(c(0.7351381, 0.05, NA)), c("0.735", "0.050", "NA"))
  expect_identical(env$.analysis_md_num(c(1000, 0, NA), digits = 0, big.mark = ","), c("1,000", "0", "NA"))
  expect_identical(env$.analysis_md_pow10(c(0.01, 10^-2.5, 1e-3)), c("10^{-2}", "10^{-2.5}", "10^{-3}"))
})
