root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
env_ccc <- new.env(parent = getNamespace("stimgate"))
for (script in c("analysis-thesis-table.R", "acs_cytof-plot_cyt.R", "acs_cytof-ccc-summary.R")) {
  source(file.path(root_dir, "scripts/r", script), local = env_ccc)
}
.ccc_summary_fixture <- function(values) {
  data.frame(
    method = rep(c("stimgate", "fbeta_default", "tailgate"), nrow(values)),
    pop = "CD4 T cells", cyt = "IFNg", stim = rep(c("p1", "mtb", "ebv")[seq_len(nrow(values))], each = 3),
    ccc = as.vector(t(values))
  )
}

test_that("CCC counts follow literal pairwise rules and include exact best ties", {
  tbl <- .ccc_summary_fixture(rbind(c(0.9, 0.6, 0.2), c(0.8, 0.8, 0.55), c(-0.4, -0.2, -0.5)))
  result <- env_ccc$.acsCytofCccSummary(tbl)
  expect_equal(result$counts$best_ccc, c(2L, 2L, 0L))
  # F-beta loses the first stratum, so beating Tailgate does not count as a win.
  expect_equal(result$counts$win_one, c(2L, 2L, 0L))
  expect_equal(result$counts$win_both, c(1L, 1L, 0L))
  expect_equal(result$counts$lose_0_1, c(1L, 1L, 3L))
  expect_equal(result$counts$lose_0_25, c(0L, 1L, 3L))
  expect_equal(result$counts$fail, c(0L, 0L, 1L))
  expect_equal(result$coverage$exact_top_ties, 1L)
  expect_equal(result$coverage$complete_three_finite, 3L)
  second_place <- env_ccc$.acsCytofCccSummary(tbl[tbl$stim == "p1", ])
  expect_equal(second_place$counts$win_one, c(1L, 0L, 0L))
})

test_that("missing CCC and undisplayed strata never become zero counts", {
  tbl <- .ccc_summary_fixture(rbind(c(0.9, NA, 0.2), c(0.8, 0.7, 0.6)))
  extra <- tbl[tbl$stim == "mtb", ]
  extra$stim <- "p4"
  extra$ccc <- 1
  result <- env_ccc$.acsCytofCccSummary(rbind(tbl, extra))
  expect_equal(result$coverage$eligible, 2L)
  expect_equal(result$coverage$complete_three_finite, 1L)
  expect_equal(result$coverage$incomplete, 1L)
  expect_equal(result$counts$best_ccc, c(1L, 0L, 0L))
  expect_equal(result$method_coverage$missing_or_nonfinite_ccc, c(0L, 1L, 0L))
  missing <- env_ccc$.acsCytofCccSummary(tbl[tbl$stim == "p1", ])
  expect_true(all(is.na(as.matrix(missing$counts[-1L]))))
  expect_equal(missing$coverage$complete_three_finite, 0L)
  empty <- env_ccc$.acsCytofCccSummary(tbl[FALSE, ])
  expect_true(all(is.na(as.matrix(empty$counts[-1L]))))
  expect_equal(empty$coverage$eligible, 0L)
})

test_that("CCC thesis source uses shared formatting and no row separators", {
  summary <- env_ccc$.acsCytofCccSummary(.ccc_summary_fixture(matrix(c(0.9, 0.6, 0.2), nrow = 1)))
  md <- env_ccc$.acsCytofCccSummaryMarkdown(summary)
  expect_match(md, "$F_{\\beta}$ (all cells)", fixed = TRUE)
  expect_match(md, "Win $\\geq 0.1$", fixed = TRUE)
  expect_false(grepl("\\\\", md, fixed = TRUE))
  expect_match(md, "{#tbl-stimgate-results-acs-ccc-wins}", fixed = TRUE)
  expect_match(md, "including exact ties", fixed = TRUE)
  expect_false(grepl("#D6D6D6", md, fixed = TRUE))
  summary$counts$best_ccc[1] <- 1000L
  expect_match(env_ccc$.acsCytofCccSummaryMarkdown(summary), "1,000", fixed = TRUE)
})

test_that("F1 winners reuse CCC comparison rules on a small finite cohort", {
  tbl <- .ccc_summary_fixture(rbind(c(0.9, 0.6, 0.2), c(0.8, 0.8, 0.55), c(0.8, NA, 0.6)))
  names(tbl)[names(tbl) == "ccc"] <- "f1"
  excluded <- tbl[tbl$stim == "p1", ]
  excluded$cyt <- "IL2"
  excluded$f1 <- NA_real_
  attr(tbl, "excludedStrata") <- excluded
  result <- env_ccc$.acsCytofCccSummary(tbl, metric = "f1")
  expect_equal(result$coverage$displayed_observed, 4L)
  expect_equal(result$coverage$eligible, 3L)
  expect_equal(result$coverage$complete_three_finite, 2L)
  expect_equal(result$coverage$incomplete, 1L)
  expect_equal(result$coverage$exact_top_ties, 1L)
  expect_equal(result$counts$best_f1, c(2L, 1L, 0L))
  expect_equal(result$counts$win_one, c(2L, 1L, 0L))
  expect_equal(result$counts$win_both, c(1L, 0L, 0L))
  expect_equal(result$counts$lose_0_1, c(0L, 1L, 2L))
  expect_equal(result$counts$lose_0_25, c(0L, 1L, 2L))
  expect_equal(result$counts$fail, c(0L, 0L, 1L))
  expect_equal(result$method_coverage$missing_or_nonfinite_f1, c(0L, 1L, 0L))
  second_place <- env_ccc$.acsCytofCccSummary(tbl[tbl$stim == "p1", ], metric = "f1")
  expect_equal(second_place$counts$win_one, c(1L, 0L, 0L))
  missing <- env_ccc$.acsCytofCccSummary(tbl[tbl$stim == "ebv", ], metric = "f1")
  expect_true(all(is.na(as.matrix(missing$counts[-1L]))))
  expect_equal(missing$coverage$complete_three_finite, 0L)
  md <- env_ccc$.acsCytofCccSummaryMarkdown(result)
  expect_match(md, "{#tbl-stimgate-results-acs-f1-wins}", fixed = TRUE)
  expect_match(md, "mean finite per-tube F1", fixed = TRUE)
  expect_match(md, "strict upper-tail reference closest to the manual", fixed = TRUE)
  expect_match(md, "does not use manual event labels", fixed = TRUE)
  expect_false(grepl("signed CCC", md, fixed = TRUE))
})
