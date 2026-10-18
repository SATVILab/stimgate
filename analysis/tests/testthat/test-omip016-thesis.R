.omip016ThesisTestEnv <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  for (file in c("omip016-methods.R", "omip016-thesis.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

test_that("OMIP-016 thesis rows preserve true-mask scores and missing outcomes", {
  env <- .omip016ThesisTestEnv()
  scored <- data.frame(
    method = c("stimgate", "fbeta", "tailgate", "manual"),
    cyt = "IL2", stim = "gag", ind = 2L,
    freq_bs_man = 2, freq_bs = c(1, 0, NA, 2),
    tp = c(2, 0, NA, 3), fp = c(1, 0, NA, 0),
    fn = c(1, 3, NA, 0), tn = c(6, 7, NA, 7),
    nCellStim = 10, countStim = c(3, 0, NA, 3), nManualStim = 3
  )
  scored$f1 <- env$.omipClassificationMetrics(scored)$f1
  rows <- env$.omip016ThesisRows(scored)
  expect_identical(rows$method, c("stimgate", "fbeta", "tailgate"))
  expect_equal(rows$f1, c(2 / 3, 0, NA))
  expect_identical(rows$population, rep("CD4 T cells", 3))
  expect_identical(rows$marker, rep("IL2", 3))
  expect_identical(rows$condition, rep("gag", 3))
  expect_equal(rows$estimate, c(1, 0, NA))
  expect_true(all(rows$scatterEligible))
  scored$f1[1] <- 1
  expect_error(env$.omip016ThesisRows(scored), "F1")
  scored$f1 <- env$.omipClassificationMetrics(scored)$f1
  scored$countStim[1] <- 4
  expect_error(env$.omip016ThesisRows(scored), "cell counts")
})

test_that("OMIP-016 thesis rows retain both StimGate settings and undefined F1", {
  env <- .omip016ThesisTestEnv()
  scored <- data.frame(
    method = c("StimGate", "Tailgate"),
    cyt = "TNFa", stim = "seb", freq_bs_man = 0, freq_bs = 0,
    tp = 0, fp = 0, fn = 0, tn = 10,
    nCellStim = 10, countStim = 0, nManualStim = 0, f1 = NA_real_
  )
  rows <- env$.omip016ThesisRows(scored)
  expect_identical(rows$method, scored$method)
  expect_true(all(is.na(rows$f1)))
  expect_error(env$.omip016ThesisRows(scored[, -1]), "complete existing score table")
})
