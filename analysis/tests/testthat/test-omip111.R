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
