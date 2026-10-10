.omipReportTestEnv <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("acs_cytof-plot_cyt.R", "omip016-methods.R", "omip111.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

.omip111ReportFixture <- function() {
  x <- rbind(c(1, 3, 0), c(1.5, 3, 1.5), c(1.5, 1.5, 0), c(3, 0, 1.5))
  colnames(x) <- c("IFNg", "IL2", "TNF")
  manual <- rbind(c(FALSE, TRUE, FALSE), c(TRUE, FALSE, FALSE),
                  c(FALSE, TRUE, FALSE), c(TRUE, FALSE, TRUE))
  colnames(manual) <- colnames(x)
  prepared <- list(input = list(stim = list(CD8 = list(
    expression = x, manual = manual, author = !manual
  ))))
  results <- data.frame(
    strain = "C57", mouse = "M1", population = "CD8", marker = colnames(x),
    method = "StimGate", sampleStim = "stim", threshold = 2, gate = 2,
    gateCyt = 1, countStim = 2, nCellStim = 4
  )
  list(prepared = prepared, results = results)
}

test_that("classification metrics retain denominators and undefined proportions", {
  env <- .omipReportTestEnv()
  counts <- data.frame(tp = c(2, 0, 0, 4), fp = c(1, 0, 0, 0),
                       fn = c(3, 2, 0, 0), tn = c(4, 5, 5, 0))
  scored <- env$.omipClassificationMetrics(counts)
  expect_equal(scored$precision, c(2 / 3, NA, NA, 1))
  expect_equal(scored$sensitivity, c(2 / 5, 0, NA, 1))
  expect_equal(scored$specificity, c(4 / 5, 1, 1, NA))
  expect_equal(scored$f1, c(4 / 8, 0, NA, 1))
  expect_equal(scored$fdp, c(1 / 3, NA, NA, 0))
  expect_equal(scored$n_selected, counts$tp + counts$fp)
  expect_equal(scored$n_manual_positive, counts$tp + counts$fn)
  expect_equal(scored$n_manual_negative, counts$tn + counts$fp)
  expect_equal(scored$n_f1, 2 * counts$tp + counts$fp + counts$fn)
})

test_that("OMIP-016 score adds TN and render-time classification metrics", {
  env <- .omipReportTestEnv()
  x <- matrix(c(0, 1, 2, 3), ncol = 1, dimnames = list(NULL, "A"))
  env$.omip016Expr <- function(gs, ind, channels) if (ind == 1) x * 0 else x
  prep <- list(gs = NULL, pre = list(responseChannels = c(A = "marker"),
    sampleMap = data.frame(file = c("uns", "stim"), stim = c("uns", "stim")),
    batchList = list(c(1, 2))),
    labels = list(uns = data.frame(marker = rep(FALSE, 4)),
                  stim = data.frame(marker = c(FALSE, TRUE, FALSE, TRUE))))
  gates <- data.frame(method = "test", ind = 2, chnl = "A", gate = 1)
  out <- env$.omip016Score(prep, gates)
  expect_equal(unname(unlist(out[c("tp", "fp", "fn", "tn")])), rep(1, 4))
  expect_equal(out$precision, 0.5)
  expect_equal(out$specificity, 0.5)
  expect_equal(out$n, 4)
  expect_equal(out$n_defined, 4)
  expect_equal(out$freq_bs, 50)
})

test_that("OMIP-111 reconstructs strict conditional calls against primary masks", {
  env <- .omipReportTestEnv()
  fixture <- .omip111ReportFixture()
  out <- env$.omip111Classification(fixture$prepared, fixture$results)
  expect_equal(out$tp, c(2, 1, 1))
  expect_equal(out$fp, c(0, 1, 1))
  expect_equal(out$fn, c(0, 1, 0))
  expect_equal(out$tn, c(2, 1, 2))
  expect_equal(out$tp + out$fp, fixture$results$countStim)
  expect_equal(out$n, rep(4, 3))
  package_calls <- stimgate:::.getPosIndByChnl(
    as.data.frame(fixture$prepared$input$stim$CD8$expression),
    data.frame(chnl = fixture$results$marker, gate = 2, gateCyt = 1),
    gateTypeCytPos = "cyt"
  )
  expect_identical(do.call(cbind, package_calls), env$.omip016Classify(
    fixture$prepared$input$stim$CD8$expression, rep(2, 3), rep(1, 3)
  ))
  # A conditional-positive marker cannot supply ordinary-positive context.
  expect_identical(env$.omip016Classify(
    fixture$prepared$input$stim$CD8$expression, rep(2, 3), rep(1, 3)
  )[3, ], c(IFNg = FALSE, IL2 = FALSE, TNF = FALSE))
  comparator <- fixture$results
  comparator$method <- "Tailgate"
  comparator$threshold <- 1.5
  out_comp <- env$.omip111Classification(fixture$prepared, comparator)
  expect_equal(out_comp$tp + out_comp$fp, c(1, 2, 0))
  comparator$threshold <- NA_real_
  failed <- env$.omip111Classification(fixture$prepared, comparator)
  expect_true(all(is.na(failed$precision)))
  expect_true(all(is.na(failed$specificity)))
})

test_that("OMIP-111 refuses missing fields, incomplete gates and count mismatches", {
  env <- .omipReportTestEnv()
  fixture <- .omip111ReportFixture()
  bad <- fixture$results
  bad$countStim[1] <- 1
  expect_error(env$.omip111Classification(fixture$prepared, bad), "marginal counts")
  bad <- fixture$results
  bad$gateCyt <- NULL
  expect_error(env$.omip111Classification(fixture$prepared, bad), "essential")
  expect_error(env$.omip111Classification(fixture$prepared, fixture$results[-1, ]),
               "Incomplete marker")
  fixture$prepared$input$stim$CD8$manual <- NULL
  expect_error(env$.omip111Classification(fixture$prepared, fixture$results), "manual masks")
})

test_that("OMIP-111 pools counts and reports defined mouse metric coverage", {
  env <- .omipReportTestEnv()
  fixture <- .omip111ReportFixture()
  out <- env$.omip111Classification(fixture$prepared, fixture$results)
  second <- out
  second$mouse <- "M2"
  second[c("tp", "fp", "fn", "tn")] <- list(rep(0, 3), rep(0, 3), rep(0, 3), rep(4, 3))
  second[names(env$.omipClassificationMetrics(second))] <- env$.omipClassificationMetrics(second)
  summary <- env$.omip111ClassificationSummary(rbind(out, second))
  row <- summary[summary$marker == "IL2", ]
  expect_equal(row$n_mice, 2)
  expect_equal(row$n_defined, 2)
  expect_equal(row$n, 8)
  expect_equal(row$n_cells_defined, 8)
  expect_equal(row$tn, 5)
  expect_equal(row$precision, 0.5)
  expect_equal(row$precision_median, 0.5)
  expect_equal(row$precision_n_defined, 1)
  expect_equal(row$specificity_median, 0.75)
  expect_equal(row$specificity, 5 / 6)
  failed <- out
  failed[c("tp", "fp", "fn", "tn")] <- NA_real_
  failed[names(env$.omipClassificationMetrics(failed))] <- env$.omipClassificationMetrics(failed)
  summary <- env$.omip111ClassificationSummary(failed)
  expect_true(all(is.na(summary$precision)))
  expect_true(all(is.na(summary$precision_median)))
  expect_equal(summary$n_defined, rep(0, 3))
})

test_that("OMIP correlations reuse population CCC and retain finite-pair coverage", {
  env <- .omipReportTestEnv()
  x <- c(1, 2, 4, 8)
  y <- 2 * x + 3
  expected <- 2 * mean((x - mean(x)) * (y - mean(y))) /
    (mean((x - mean(x))^2) + mean((y - mean(y))^2) + (mean(x) - mean(y))^2)
  sample_ccc <- 2 * stats::cov(x, y) /
    (stats::var(x) + stats::var(y) + (mean(x) - mean(y))^2)
  expect_gt(abs(expected - sample_ccc), 0.01)
  data <- data.frame(method = "StimGate", stim = "s", cyt = "c",
                     est = c(x, NA, Inf), ref = c(y, 1, 2))
  table <- env$.omipFrequencyCorrelationTable(data, c("stim", "cyt"), "est", "ref")
  expect_equal(table$n, c(4, 4))
  expect_equal(table$n_total, c(6, 6))
  expect_equal(table$pearson, c(1, 1))
  expect_equal(table$ccc, c(expected, expected))
  expect_equal(env$.acsCytofValidationCcc(x, y), expected)
  for (values in list(rep(1, 6), c(1, 2, NA, NA, NA, NA))) {
    data$est <- values
    table <- env$.omipFrequencyCorrelationTable(data, c("stim", "cyt"), "est", "ref")
    expect_true(all(is.na(table$pearson)))
    expect_true(all(is.na(table$ccc)))
  }
  data <- data.frame(method = rep(c("StimGate default", "StimGate shifted peak"), each = 4),
    population = "CD8", strain = rep(c("C57", "BALB"), 4), marker = "TNF",
    est = rep(x, 2), ref = rep(y, 2))
  table <- env$.omipFrequencyCorrelationTable(data,
    c("strain", "population", "marker"), "est", "ref", "population")
  expect_equal(table$n[table$scope == "pooled"], c(4, 4))
  expect_true(all(is.na(table$ccc[table$scope == "stratum"])))
  expect_setequal(table$method, c("StimGate default", "StimGate shifted peak"))
})

test_that("new OMIP report chunks guard plotting and missing results", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  files <- c("14-real-compare-omip111.qmd", "14b-real-compare-omip111-shifted-peak.qmd",
             "15-real-compare-omip016.qmd", "15b-real-compare-omip016-shifted-peak.qmd")
  for (file in files) {
    lines <- readLines(file.path(root, "analysis", file), warn = FALSE)
    labels <- c(if (startsWith(file, "14")) "classification-tables" else "classification",
                "frequency-correlations")
    for (label in labels) {
      start <- which(lines == paste0("#| label: ", label))
      expect_length(start, 1)
      end <- start + which(lines[(start + 1L):length(lines)] == "```")[[1]]
      code <- parse(text = lines[(start + 1L):(end - 1L)])
      for (plots in c(FALSE, TRUE)) {
        env <- new.env(parent = baseenv())
        env$run_plots <- plots
        expect_no_error(eval(code, envir = env))
        expect_false(exists("cls", envir = env, inherits = FALSE))
        expect_false(exists("correlations", envir = env, inherits = FALSE))
      }
    }
  }
})

test_that("agreement overview pools counts before proportions and keeps failed gates out", {
  env <- .omipReportTestEnv()
  counts <- data.frame(
    population = "CD4", method = c("A", "A", "B"),
    tp = c(1, 3, NA), fp = c(1, 0, NA), fn = c(0, 1, NA), tn = c(8, 6, NA)
  )
  correlations <- data.frame(
    population = "CD4", method = c("A", "B", "A"), scope = c("pooled", "pooled", "stratum"),
    n = c(4, 0, 2), pearson = c(0.5, NA, 0.9), ccc = c(0.4, NA, 0.8)
  )
  out <- env$.omipAgreementOverview(counts, correlations, "population")
  a <- out[out$method == "A", ]
  expect_equal(c(a$n_outcomes, a$n_defined, a$tp, a$fp), c(2, 2, 4, 1))
  expect_equal(a$sensitivity, 4 / 5)
  expect_equal(a$precision, 4 / 5)
  expect_equal(a$specificity, 14 / 15)
  expect_equal(a$f1, 8 / 10)
  expect_equal(c(a$n_pairs, a$pearson), c(4, 0.5))
  b <- out[out$method == "B", ]
  expect_equal(b$n_defined, 0)
  expect_true(is.na(b$sensitivity) && is.na(b$precision) && is.na(b$pearson))
  printed <- capture.output(env$.omipAgreementKable(out))
  expect_true(any(grepl("80.0%", printed, fixed = TRUE)))
  expect_false(any(grepl("^\\|.*\\btp\\b", printed)))
})
