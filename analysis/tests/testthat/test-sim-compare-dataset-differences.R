.dataset_difference_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  source(file.path(root, "scripts", "r", "sim-compare-freq_bs.R"), local = env)
  env
}

.dataset_error_fixture <- function() {
  tidyr::expand_grid(
    scenario = c("a", "b"), iter = 1:5, sample = c("1", "2"),
    method = c("stimgate", "fbeta", "tailgate")
  ) |>
    dplyr::mutate(
      approach = .data$method,
      propRespTruth = 0.5,
      propRespEst = 0.5 + dplyr::case_when(
        .data$method == "stimgate" ~ .data$iter *
          ifelse(.data$sample == "1", 0.1, 0.3),
        .data$method == "fbeta" ~ 0.1,
        TRUE ~ 0.2
      )
    )
}

test_that("frequency comparisons pair dataset means and use between-dataset spread", {
  env <- .dataset_difference_env()
  tbl <- .dataset_error_fixture()
  out <- env$.simCompareDatasetDifferences(tbl, c("scenario", "method", "approach"))
  a <- out[out$scenario == "a" & out$method == "fbeta" & out$outcome == "abs_error", ]
  # Each dataset's mean StimGate error is .2 * iter; F-beta's is .1.
  delta <- 0.2 * (1:5) - 0.1
  half <- 1.96 * stats::sd(delta) / sqrt(5)
  expect_equal(a$mean_difference, mean(delta))
  expect_equal(a$lower, mean(delta) - half)
  expect_equal(a$upper, mean(delta) + half)
  expect_equal(a$n_dataset, 5L)
  expect_equal(a$n_pair, 5L)
  expect_equal(a$pair_coverage, 1)
  expect_equal(a$tube_coverage_stimgate, 1)
  rel <- out[out$scenario == "a" & out$method == "fbeta" & out$outcome == "abs_rel_error", ]
  expect_equal(rel$mean_difference, 2 * a$mean_difference)
  expect_equal(rel$lower, 2 * a$lower)
  expect_equal(nrow(out), 8L)

  # Row order and dependent tube counts do not change the dataset weights.
  extra <- tbl[tbl$iter == 1L & tbl$sample == "1", ]
  extra$sample <- "3"
  extra$propRespEst[extra$method == "stimgate"] <- 0.7
  more <- dplyr::bind_rows(tbl, extra)
  expect_equal(env$.simCompareDatasetDifferences(more, "scenario")$mean_difference,
    out$mean_difference)
  expect_equal(env$.simCompareDatasetDifferences(tbl[nrow(tbl):1, ], "scenario"), out)
  expect_error(env$.simCompareDatasetDifferences(dplyr::bind_rows(tbl, tbl[1, ]), "scenario"),
    "one primary row")
})

test_that("incomplete pairs retain coverage and suppress small-D intervals", {
  env <- .dataset_difference_env()
  tbl <- .dataset_error_fixture()
  tbl$propRespEst[tbl$method == "fbeta" & tbl$iter == 5L] <- NA_real_
  out <- env$.simCompareDatasetDifferences(tbl, "scenario")
  fb <- out[out$method == "fbeta", ]
  expect_true(all(fb$n_dataset == 5L))
  expect_true(all(fb$n_pair == 4L))
  expect_true(all(fb$pair_coverage == 0.8))
  expect_true(all(fb$tube_coverage_competitor == 0.8))
  expect_true(all(is.na(fb$lower) & is.na(fb$upper)))
  tbl$propRespTruth <- 0
  zero <- env$.simCompareDatasetDifferences(tbl, "scenario", outcomes = "abs_rel_error")
  expect_true(all(zero$n_pair == 0L))
  expect_true(all(is.na(zero$mean_difference)))
  expect_true(all(zero$pair_coverage == 0))
  expect_error(env$.simCompareDatasetDifferences(tbl, "scenario", outcomes = "unknown"),
    "Unknown")
  expect_error(env$.simCompareDatasetDifferences(tbl[setdiff(names(tbl), "iter")], "scenario"),
    "require")
})

test_that("classification comparisons use dataset medians and preserve undefined FDP", {
  env <- .dataset_difference_env()
  tbl <- .dataset_error_fixture() |>
    dplyr::mutate(
      nFalsePos = ifelse(.data$method == "stimgate",
        .data$iter + ifelse(.data$sample == "1", -1, 1), 2),
      nTruePos = 10 - .data$nFalsePos,
      nFalseNeg = .data$nFalsePos,
      nTrueNeg = 10 - .data$nFalsePos
    )
  out <- env$.simCompareDatasetDifferences(tbl, "scenario", c("fdp", "sensitivity"))
  fdp <- out[out$scenario == "a" & out$method == "fbeta" & out$outcome == "fdp", ]
  delta <- (1:5) / 10 - 0.2
  expect_equal(fdp$mean_difference, mean(delta))
  expect_equal(fdp$upper, mean(delta) + 1.96 * stats::sd(delta) / sqrt(5))
  sens <- out[out$scenario == "a" & out$method == "fbeta" & out$outcome == "sensitivity", ]
  expect_equal(sens$mean_difference, -fdp$mean_difference)
  expect_equal(sens$lower, -fdp$upper)
  empty <- tbl$iter == 5 & tbl$method == "stimgate"
  tbl$nTruePos[empty] <- 0
  tbl$nFalsePos[empty] <- 0
  tbl$nFalseNeg[empty] <- 10
  out <- env$.simCompareDatasetDifferences(tbl, "scenario", c("fdp", "sensitivity"))
  expect_true(all(out$n_pair[out$outcome == "fdp"] == 4L))
  expect_true(all(out$n_pair[out$outcome == "sensitivity"] == 5L))
})

test_that("both reports call dataset comparisons and label pooled tubes", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  for (document in c("7-sim-compare-freq_bs.qmd", "8-sim-compare-freq_bs-batch.qmd")) {
    text <- paste(readLines(file.path(root, "analysis", document)), collapse = "\n")
    expect_match(text, "tube-level distributions", fixed = TRUE)
    expect_match(text, ".simCompareDatasetDifferences(", fixed = TRUE)
    expect_match(text, "Dataset-level paired differences: StimGate minus competitor.", fixed = TRUE)
  }
})
