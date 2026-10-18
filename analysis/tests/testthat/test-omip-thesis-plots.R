.omipThesisTestEnv <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-plot-style.R", "acs_cytof-plot_cyt.R", "omip016-methods.R",
                 "omip-thesis-plots.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

.omipThesisFixture <- function() {
  data.frame(method = rep(c("StimGate", "F-beta", "Tailgate"), each = 3),
    population = "CD4", marker = "IFNg", condition = "C57",
    reference = rep(c(-1, 0, 2), 3), estimate = rep(c(-2, 0, 4), 3),
    f1 = rep(c(0.2, 0.8, NA_real_), 3), scatterEligible = rep(c(TRUE, FALSE, TRUE), 3))
}

test_that("OMIP heatmap summaries use tube means and explicit missing coverage", {
  env <- .omipThesisTestEnv()
  rows <- .omipThesisFixture()
  summary <- env$.omipThesisSummary(rows)
  expect_equal(summary$strata$f1, rep(0.5, 3))
  expect_equal(summary$strata$n_f1_total, rep(3L, 3))
  expect_equal(summary$strata$n_f1_finite, rep(2L, 3))
  expect_equal(summary$strata$n_f1_missing, rep(1L, 3))
  expect_equal(summary$strata$n_scatter_eligible, rep(2L, 3))
  expect_equal(summary$strata$n, rep(3L, 3))
  expect_equal(summary$strata$pearson, rep(1, 3))
  expect_true(all(summary$strata$ccc < 1))
  rows$f1 <- NA_real_
  expect_true(all(is.na(env$.omipThesisSummary(rows)$strata$f1)))
  rows$f1[[1]] <- 2
  expect_error(env$.omipThesisSummary(rows), "Invalid")
})

test_that("OMIP unmeasured heatmap cells stay blank and undefined observed scores stay NA", {
  env <- .omipThesisTestEnv()
  rows <- .omipThesisFixture()
  extra <- rows
  extra$marker <- "IL2"
  extra$population <- "CD8"
  extra$f1 <- NA_real_
  summary <- env$.omipThesisSummary(rbind(rows, extra))$strata
  labels <- stats::setNames(unique(rows$method), unique(rows$method))
  plot <- env$.omipThesisHeatmap(summary, "f1", "viridis", labels)
  expect_equal(nrow(plot$data), 6L)
  expect_equal(sum(plot$data$label == "NA"), 3L)
  expect_false(any(plot$data$marker == "IL2" & plot$data$population == "CD4"))
  expect_equal(plot$scales$get_scales("fill")$limits, c(0, 1))
  expect_equal(plot$theme$panel.background$fill, "white")
  expect_equal(plot$theme$legend.position, "bottom")
  expect_equal(plot$labels$x, "Cytokine")
  signed <- env$.omipThesisHeatmap(summary, "ccc", "viridis", labels)
  fill <- signed$scales$get_scales("fill")
  expect_equal(fill$limits, c(-1, 1))
  expect_equal(fill$palette(c(0, 0.5)), rep(env$.acsCytofF1Palette(0), 2))
  expect_equal(nrow(signed$layers[[3]]$data), 2L)
})

test_that("OMIP111 adaptor joins by mouse keys without losing rows or gate settings", {
  env <- .omipThesisTestEnv()
  results <- data.frame(strain = "C57", mouse = c("M1", "M2"), population = "CD4",
    marker = "IFNg", method = c("StimGate", "Tailgate"),
    propBs = c(-0.01, 0.02), manualBs = c(0.01, 0.02), referenceAccepted = c(TRUE, FALSE))
  classification <- results[2:1, c("strain", "mouse", "population", "marker", "method")]
  classification$f1 <- c(0.9, 0.2)
  classification$tp <- c(9, 2)
  classification$fp <- c(1, 8)
  classification$fn <- c(1, 8)
  classification$tn <- c(4, 4)
  classification$n_f1 <- c(20, 20)
  normalized <- env$.omip111ThesisRows(results, classification)
  expect_equal(normalized$f1, c(0.2, 0.9))
  expect_equal(normalized$tp, c(2, 9))
  expect_equal(normalized$n_f1, c(20, 20))
  expect_false(any(grepl("[.][xy]$", names(normalized))))
  expect_equal(normalized$estimate, c(-1, 2))
  expect_equal(normalized$reference, c(1, 2))
  expect_identical(normalized$scatterEligible, c(TRUE, FALSE))
  expect_equal(normalized$method, results$method)
  expect_error(env$.omip111ThesisRows(results, classification[1, ]), "do not match")
  expect_error(env$.omip111ThesisRows(rbind(results, results[1, ]), classification), "duplicate")
  results$n_f1 <- c(1, 1)
  expect_error(env$.omip111ThesisRows(results, classification), "audit fields disagree")
})

test_that("OMIP scatter preserves signed frequencies and its existing eligibility cohort", {
  env <- .omipThesisTestEnv()
  rows <- .omipThesisFixture()
  plot <- env$.omipThesisScatterPanel(rows, "IFNg", list(IFNg = c(-2.5, 4.5)),
    c(C57 = "#0072B2"), c(C57 = 16), "Projected reference")
  expect_equal(nrow(plot$data), 6L)
  expect_true(any(plot$data$reference < 0))
  expect_equal(plot$scales$get_scales("x")$limits, c(-2.5, 4.5))
  expect_equal(plot$scales$get_scales("y")$limits, c(-2.5, 4.5))
  expect_equal(plot$theme$panel.background$fill, "white")
})
