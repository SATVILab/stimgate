local({
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "scripts", "r", "omip111.R"), env)
  sys.source(file.path(root, "scripts", "r", "omip111-debug.R"), env)
  env$.analysis_method_colours <- c(fbeta = "blue", tailgate = "green")

  testthat::test_that("diagnostics overlay actual tube-specific author cutoffs and conditional gates", {
    rows <- data.frame(
      method = c("StimGate", "F-beta", "Tailgate"),
      threshold = c(2, 3, 4), gateCyt = c(1.5, NA, NA),
      sampleStim = "stim", sampleUns = "uns"
    )
    refs <- data.frame(sample = c("uns", "stim"), cytokineLowerRaw = c(150, 300))
    lines <- env$.omip111DebugLines(rows, refs, 150)
    testthat::expect_equal(lines$x[lines$line == "author stim"], asinh(2))
    testthat::expect_equal(lines$x[lines$line == "author uns"], asinh(1))
    testthat::expect_equal(lines$x[lines$line == "conditional gate"], 1.5)
    rows$gateCyt[1] <- 2
    testthat::expect_false("conditional gate" %in% env$.omip111DebugLines(rows, refs, 150)$line)
  })

  testthat::test_that("diagnostic pages print every page with its note", {
    dir <- withr::local_tempdir()
    writeBin(as.raw(1:8), file.path(dir, "page.png"))
    index <- expand.grid(
      mouse = c("M1", "M2"), marker = "TNF", population = c("CD4", "CD8"),
      strain = c("C57", "BALB"), stringsAsFactors = FALSE
    )
    index$page <- 1L
    index$png <- "page.png"
    index$note <- paste0(": note ", seq_len(nrow(index)))
    out <- paste(utils::capture.output(env$.omip111DiagnosticPages(index, dir)), collapse = "\n")
    testthat::expect_equal(lengths(regmatches(out, gregexpr("<details>", out))), nrow(index))
    testthat::expect_equal(lengths(regmatches(out, gregexpr("data:image/png", out))), nrow(index))
    testthat::expect_true(all(vapply(index$note, grepl, logical(1), x = out, fixed = TRUE)))
    testthat::expect_match(out, "### BALB CD44-high CD8", fixed = TRUE)
    # Both fit populations are grouped under their own headings.
    both <- rbind(transform(index, fit = "cd44"), transform(index, fit = "nonnaive"))
    out <- paste(utils::capture.output(env$.omip111DiagnosticPages(both, dir)), collapse = "\n")
    testthat::expect_equal(lengths(regmatches(out, gregexpr("<details>", out))), nrow(both))
    testthat::expect_match(out, "### Fitted on all non-naïve cells", fixed = TRUE)
    testthat::expect_match(out, "#### BALB non-naïve CD8", fixed = TRUE)
    testthat::expect_match(out, "#### C57 CD44-high CD4", fixed = TRUE)
  })

  testthat::test_that("OMIP-111 pairs by mouse and condition and rejects incomplete pairs", {
    files <- unlist(lapply(c("C57", "BALB"), function(strain) {
      unlist(lapply(1:5, function(mouse) {
        paste0(
          c("E02", "F02"), " ", strain, "_M", mouse,
          c("_BFA", "_Stim"), "_TS WLSM.fcs"
        )
      }))
    }))
    map <- env$.omip111SampleMap(rev(files))
    testthat::expect_equal(length(unique(map$mouse)), 10L)
    testthat::expect_equal(map$condition, rep(c("uns", "stim"), 10))
    testthat::expect_error(env$.omip111SampleMap(files[-1]), "20")
    files[2] <- files[1]
    testthat::expect_error(env$.omip111SampleMap(files), "exactly one")
  })

  testthat::test_that("OMIP-111 retains negative net frequencies and zero references", {
    counts <- env$.omip111Counts(1, 100, 2, 100, 0, 0)
    testthat::expect_equal(counts$propBs, -0.01)
    testthat::expect_equal(counts$manualBs, 0)
    failed <- env$.omip111Counts(NA_real_, 100, NA_real_, 100, 0, 0)
    testthat::expect_true(is.na(failed$propBs))
    testthat::expect_error(env$.omip111Counts(1, 0, 2, 100, 0, 0), "positive")
  })

  testthat::test_that("author figure labels follow the printed precision", {
    testthat::expect_equal(
      env$.omip111AuthorFigureLabel(c(5.843, 10.312, 0.0371, 0.0603, 0.117)),
      c("5.84", "10.3", "0.037", "0.060", "0.12")
    )
  })

  testthat::test_that("fit and evaluation populations rescore the same gates", {
    m <- function(x) matrix(x, dimnames = list(NULL, "IFNg"))
    tube <- function(x, parent, nonNaive, manual, cut) list(CD8 = list(
      expression = m(x), parent = m(parent), manual = m(manual), cytokineCut = m(cut),
      nonNaive = nonNaive
    ))
    # The CD44+ parent lies within the non-naive gate, as the authors gated.
    prepared <- list(input = list(
      stim = tube(c(3, 3, 3, 0, 0), c(TRUE, FALSE, FALSE, TRUE, FALSE),
        c(TRUE, TRUE, TRUE, TRUE, FALSE), c(TRUE, FALSE, FALSE, FALSE, FALSE),
        c(TRUE, TRUE, TRUE, FALSE, FALSE)),
      uns = tube(c(3, 0, 0, 0), c(FALSE, FALSE, TRUE, TRUE), rep(TRUE, 4), rep(FALSE, 4),
        c(TRUE, FALSE, FALSE, FALSE))
    ))
    base <- data.frame(strain = "C57", mouse = "C57_M1", population = "CD8", marker = "IFNg",
      method = "StimGate", sampleStim = "stim", sampleUns = "uns", threshold = 2)
    # Saved counts on each fit's own population.
    results <- cbind(base, countStim = 1, countUns = 0, referenceAccepted = FALSE)
    nonNaive <- cbind(base, countStim = 3, countUns = 1)
    score <- function(fit, evaluate, reference = "bivariate") {
      env$.omip111FitEval(prepared, results, nonNaive, fit, evaluate, reference)
    }
    # A and C: 1 of 2 CD44+ stimulated cells, none of the control's.
    testthat::expect_equal(score("cd44", "cd44")$propBs, 1 / 2)
    testthat::expect_equal(score("nonnaive", "cd44")$propBs, 1 / 2)
    testthat::expect_equal(score("cd44", "cd44")$manualBs, 1 / 2)
    # B and D: the same CD44+ positives over every non-naive cell; CD44-low
    # cells above the gate are never positive (unrestricted D would be 1/2).
    testthat::expect_equal(score("cd44", "nonnaive")$propBs, 1 / 4)
    testthat::expect_equal(score("nonnaive", "nonnaive")$propBs, 1 / 4)
    testthat::expect_equal(score("nonnaive", "nonnaive")$manualBs, 1 / 4)
    # The cytokine cutoff alone is also counted within the CD44+ parent.
    testthat::expect_equal(score("cd44", "cd44", "univariate")$manualBs, 1 / 2)
    testthat::expect_equal(score("nonnaive", "nonnaive", "univariate")$manualBs, 1 / 4)
    # Only the 2-D reconstruction carries the import tolerance.
    testthat::expect_false(score("cd44", "cd44")$referenceAccepted)
    testthat::expect_true(score("cd44", "cd44", "univariate")$referenceAccepted)
    results$countStim <- 3
    testthat::expect_error(score("cd44", "cd44"), "saved counts")
  })

  testthat::test_that("float32 inputs share exactly the package precision", {
    x <- c(-0.123456789, 0, 1 / 3, 1.23456789)
    y <- env$.omip111Float32(x)
    testthat::expect_identical(env$.omip111Float32(y), y)
    testthat::expect_true(y[1] < 0)
    testthat::expect_equal(y[2], 0)
  })

  testthat::test_that("OMIP-111 exposes failures and fallback coverage", {
    x <- data.frame(
      strain = "C57", population = "CD4", marker = "IFNg",
      method = "Tailgate", propBs = c(0.1, NA), error = c(NA, "error"),
      thresholdFallbackUsed = c(TRUE, FALSE), errorPp = c(2, NA), absoluteErrorPp = c(2, NA)
    )
    summary <- env$.omip111Summary(x)
    testthat::expect_equal(summary$nMice, 2L)
    testthat::expect_equal(summary$nFinite, 1L)
    testthat::expect_equal(summary$nAgreementFinite, 1L)
    testthat::expect_equal(summary$nErrors, 1L)
    testthat::expect_equal(summary$nFallback, 1L)
    testthat::expect_equal(summary$meanAbsoluteErrorPp, 2)
  })
})
