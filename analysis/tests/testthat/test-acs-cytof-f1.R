.acs_f1_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "acs_cytof-plot_cyt.R", "acs_cytof-f1.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

test_that("manual net frequency gives strict attainable gates and cell F1", {
  env <- .acs_f1_env()
  out <- env$.acsCytofF1One(1:10, 30, c(7, 5, 9, 10, NA, 7), c(rep(FALSE, 5), TRUE))
  expect_equal(out$referenceGate, rep(7, 6))
  expect_equal(out$nReference, rep(3L, 6))
  expect_equal(out$f1, c(1, 0.75, 0.5, 0, NA, NA))
  expect_equal(out$tp[1:4], c(3L, 3L, 1L, 0L))
  expect_equal(out$fp[1:4], c(0L, 2L, 0L, 0L))
  expect_equal(out$fn[1:4], c(0L, 0L, 2L, 3L))
  tied <- env$.acsCytofF1One(c(0, 0, 1, 1), 25, 0, FALSE)
  expect_equal(tied$referenceGate, 0)
  expect_equal(tied$nReference, 2L)
  expect_equal(tied$referenceFrequencyError, 25)
  zero <- env$.acsCytofF1One(1:10, 0, c(10, 0), c(FALSE, FALSE))
  expect_equal(zero$f1, c(NA, 0))
  expect_equal(zero$referenceGate, c(10, 10))
  full <- env$.acsCytofF1One(1:10, 100, 0, FALSE)
  expect_equal(full$referenceGate, -Inf)
  expect_equal(full$f1, 1)
  expect_true(is.na(env$.acsCytofF1One(numeric(), 10, 0, FALSE)$f1))
  expect_error(env$.acsCytofF1One(c(1, NA), 10, 0, FALSE), "Non-finite")
  for (target in c(-1, NA_real_, 101)) {
    invalid <- env$.acsCytofF1One(1:10, target, 7, FALSE)
    expect_true(is.na(invalid$f1))
    expect_true(is.na(invalid$nReference))
    expect_equal(invalid$nPredicted, 3L)
    expect_false(invalid$referenceReason == "valid")
    expect_equal(invalid$freqManualNet, target)
  }
  expect_equal(env$.acsCytofF1One(1:10, -1, 7, FALSE)$referenceReason,
               "negative_manual_net_frequency")
})

test_that("F1 summaries retain failures, eligibility and a common 0 to 1 scale", {
  env <- .acs_f1_env()
  d <- tidyr::expand_grid(
    method = c("stimgate", "fbeta", "fbeta_default", "tailgate"),
    pop = c("CD4 T cells", "B cells"), cyt = "IFNg", stim = c("p1", "p4"),
    SampleID = letters[1:4]
  )
  d$freq_stim_man <- 1
  d$freq_uns_man <- 0.01
  d$freq_bs_man <- 0.99
  d$f1 <- rep(c(0.5, 1, NA, NA), nrow(d) / 4)
  d$savedCountAvailable <- TRUE
  d$thresholdFailed <- rep(c(FALSE, FALSE, TRUE, FALSE), nrow(d) / 4)
  d$referenceReason <- rep(c("valid", "valid", "valid", "negative_manual_net_frequency"), nrow(d) / 4)
  d$freq_bs_man[d$pop == "B cells"] <- 0
  summary <- env$.acsCytofF1Summary(d)
  expect_equal(summary$f1, rep(0.75, nrow(summary)))
  expect_true(all(summary$n_total == 4 & summary$n_failed == 1 &
    summary$n_scored == 2 & summary$n_undefined == 1 & summary$n_invalid_reference == 1))
  expect_true(all(attr(summary, "excludedStrata")$pop == "B cells"))
  expect_true(all(is.na(as.matrix(attr(summary, "excludedStrata")[c("f1", "f1_median", "f1_q10")]))))
  plot <- env$.acsCytofThesisPlotCorrelation(summary, "f1", FALSE)
  expect_no_error(ggplot2::ggplotGrob(plot))
  expect_equal(plot$scales$get_scales("fill")$limits, c(0, 1))
  expect_equal(length(unique(plot$data$method)), 3L)
  expect_true(all(as.character(plot$data$stim) != "Non-secreted Mtb proteins"))
  expect_true(all(plot$data$cellLabel[plot$data$excluded] == ""))
})

test_that("F1 mean, median and tenth percentile summarise finite tube scores", {
  env <- .acs_f1_env()
  d <- tibble::tibble(
    method = "stimgate", pop = "CD4 T cells", cyt = "IFNg", stim = "p1",
    SampleID = letters[1:4], freq_stim_man = 1, freq_uns_man = 0.01,
    freq_bs_man = 0.99, f1 = c(0.1, 0.2, 0.9, NA_real_),
    thresholdFailed = c(FALSE, FALSE, FALSE, TRUE), savedCountAvailable = TRUE,
    referenceReason = "valid"
  )
  summary <- env$.acsCytofF1Summary(d)
  expect_equal(summary$f1, 0.4)
  expect_equal(summary$f1_median, 0.2)
  expect_equal(summary$f1_q10, 0.12)
  expect_equal(summary$n_total, 4L)
  expect_equal(summary$n_scored, 3L)
  expect_equal(summary$n_failed, 1L)
  expect_equal(summary$n_undefined, 0L)
  d$f1 <- NA_real_
  d$referenceReason[1] <- "negative_manual_net_frequency"
  empty <- env$.acsCytofF1Summary(d)
  expect_true(all(is.na(as.matrix(empty[c("f1", "f1_median", "f1_q10")]))))
  expect_equal(empty$n_total, 4L)
  expect_equal(empty$n_scored, 0L)
  expect_equal(empty$n_failed, 1L)
  expect_equal(empty$n_invalid_reference, 1L)
  expect_equal(empty$n_undefined, 3L)
})


test_that("reconstructed F1 gates reproduce saved counts and percent frequencies", {
  env <- .acs_f1_env()
  scores <- env$.acsCytofF1One(1:10, 30, c(7, NA), c(FALSE, TRUE))
  d <- tibble::tibble(nCellStim = c(10, 10), countStim = c(3, NA),
                      freq_stim_auto = c(30, NA))
  expect_equal(env$.acsCytofF1ValidateCounts(d, scores, 10), c(0, NA))
  wrong <- d
  wrong$nCellStim[1] <- NA
  expect_error(env$.acsCytofF1ValidateCounts(wrong, scores, 10), "cell count")
  wrong <- d
  wrong$countStim[1] <- NA
  expect_equal(env$.acsCytofF1ValidateCounts(wrong, scores, 10), c(NA_real_, NA_real_))
  wrong <- d
  wrong$countStim[1] <- 4
  expect_error(env$.acsCytofF1ValidateCounts(wrong, scores, 10), "positive counts")
  wrong <- d
  wrong$freq_stim_auto[1] <- 0.3
  expect_error(env$.acsCytofF1ValidateCounts(wrong, scores, 10), "stimulated frequencies")
  wrong$freq_stim_auto[1] <- NA
  expect_error(env$.acsCytofF1ValidateCounts(wrong, scores, 10), "stimulated frequencies")
})


test_that("F1 table requires saved denominator and stimulated frequency fields", {
  env <- .acs_f1_env()
  d <- tibble::tibble(
    method = "stimgate", pop = "CD4 T cells", cyt = "IFNg", stim = "p1",
    SampleID = "a", freq_bs_auto = 30, freq_bs_man = 30,
    freq_stim_man = 31, freq_uns_man = 1, popCode = "cd4", ind = "1",
    threshold = 7, thresholdFailed = FALSE, nCellStim = 10, countStim = 3,
    freq_stim_auto = 30
  )
  for (field in c("nCellStim", "countStim", "freq_stim_auto")) {
    missing <- d[, setdiff(names(d), field)]
    expect_error(env$.acsCytofF1Table(missing, "unused"), "requires saved gates")
  }
})
