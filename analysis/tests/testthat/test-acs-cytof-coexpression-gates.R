local({
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = globalenv())
  for (f in c("analysis-plot-style.R", "coexpression-gates.R", "acs_cytof-coexpression-gates.R")) {
    sys.source(file.path(root, "scripts", "r", f), env)
  }
  # Two donors' tube pairs: IFNg and TNF co-expressed by responders (d1) or
  # no response (d2); IL2 has no gate.
  tube <- function(donor, respond) {
    withr::with_seed(if (respond) 1 else 2, {
      neg <- function(n) data.frame(IFNg = abs(rnorm(n, 0.4, 0.3)), TNF = abs(rnorm(n, 0.6, 0.4)),
        IL2 = abs(rnorm(n, 0.5, 0.3)))
      resp <- data.frame(IFNg = rnorm(200, 3.5, 0.4), TNF = runif(200, 2, 5), IL2 = abs(rnorm(200, 0.5, 0.3)))
      list(stim = if (respond) rbind(neg(20000), resp) else neg(20000), uns = neg(20000),
        gate = c(IFNg = 2, TNF = 4.5, IL2 = NA))
    })
  }
  out <- lapply(list(d1 = TRUE, d2 = FALSE), function(r) NULL)
  out$d1 <- env$.acsCoexGatesTube(tube("d1", TRUE), "cd4", "p1", "d1")
  out$d2 <- env$.acsCoexGatesTube(tube("d2", FALSE), "cd4", "p1", "d2")
  low <- do.call(rbind, lapply(out, `[[`, "low"))
  counts <- do.call(rbind, lapply(out, `[[`, "counts"))

  testthat::test_that("tube results cover every pair, rule and combination", {
    testthat::expect_equal(nrow(low), 2L * 6L)
    testthat::expect_true(low$lowered[low$donor == "d1" & low$a == "IFNg" & low$b == "TNF"])
    testthat::expect_false(any(low$lowered[low$donor == "d2"]))
    # Cells added are new positives for b, so only where a gate was lowered.
    testthat::expect_true(all(low$addedStim[!low$lowered] == 0L))
    testthat::expect_gt(low$addedStim[low$donor == "d1" & low$a == "IFNg" & low$b == "TNF"], 0L)
    tot <- stats::aggregate(n ~ donor + rule + tube, counts, sum)
    testthat::expect_true(all(tot$n[tot$tube == "uns"] == 20000L))
    testthat::expect_equal(length(unique(counts$combn)), 8L)
  })

  testthat::test_that("single- and multi-positive frequencies move between rules", {
    paired <- env$.acsCoexGatesPaired(env$.acsCoexGatesSingleMulti(env$.acsCoexGatesFreq(counts)))
    d1 <- paired[paired$donor == "d1" & paired$cytokine == "any", ]
    testthat::expect_gt(d1$diffBs[d1$type == "multi"], 0)
    testthat::expect_lt(d1$diffBs[d1$type == "single"], 0)
    d2 <- paired[paired$donor == "d2", ]
    testthat::expect_equal(d2$diffBs, rep(0, nrow(d2)))
    flags <- env$.acsCoexGatesFlags(low, paired)
    testthat::expect_equal(sort(flags$donor), c("d1", "d2"))
    testthat::expect_false(any(flags$flagControl))
  })

  testthat::test_that("COMPASS input puts the all-negative category last", {
    inp <- env$.acsCoexGatesCompassInput(counts, "ordinary")
    testthat::expect_equal(dim(inp$n_s), c(2L, 8L))
    testthat::expect_equal(colnames(inp$n_s)[[8]], "!IFNg&!TNF&!IL2")
    testthat::expect_identical(colnames(inp$n_u), colnames(inp$n_s))
    testthat::expect_equal(sum(inp$n_u["d1", ]), 20000)
    keep <- env$.acsCoexGatesCompassKeep(counts, minCell = 5L, minDonor = 1L)
    testthat::expect_true("!IFNg&!TNF&!IL2" %in% keep)
    testthat::expect_true("IFNg&TNF&!IL2" %in% keep)
    kept <- env$.acsCoexGatesCompassInput(counts, "coexpr", keep)
    testthat::expect_equal(colnames(kept$n_s)[[ncol(kept$n_s)]], "!IFNg&!TNF&!IL2")
  })

  testthat::test_that("Analysis 18 plots build without titles", {
    paired <- env$.acsCoexGatesPaired(env$.acsCoexGatesSingleMulti(env$.acsCoexGatesFreq(counts)))
    scores <- data.frame(pop = "cd4", stim = "p1", donor = rep(c("d1", "d2"), 2),
      rule = rep(c("ordinary", "coexpr"), each = 2), fs = c(0.1, 0, 0.2, 0), pfs = c(0.1, 0, 0.15, 0))
    dat <- tube("d1", TRUE)
    lw <- low[low$donor == "d1" & low$a == "IFNg" & low$b == "TNF", ]
    plots <- list(
      env$.acsCoexGatesScatterPlot(paired), env$.acsCoexGatesDiffPlot(paired),
      env$.acsCoexGatesCompassPlot(scores, "fs"), env$.acsCoexGatesCompassPlot(scores, "pfs"),
      env$.acsCoexGatesTubeHexPlot(list(list(dat = dat, a = "IFNg", b = "TNF", cut = lw$cut,
        condCut = lw$condCut, label = "d1 p1")))
    )
    for (p in plots) {
      testthat::expect_s3_class(p, "ggplot")
      testthat::expect_no_error(ggplot2::ggplot_build(p))
      testthat::expect_null(p$labels$title)
    }
  })
})
