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
  out <- list(
    d1 = env$.acsCoexGatesTube(tube("d1", TRUE), "cd4", "p1", "d1"),
    d2 = env$.acsCoexGatesTube(tube("d2", FALSE), "cd4", "p1", "d2")
  )
  low <- do.call(rbind, lapply(out, `[[`, "low"))
  counts <- do.call(rbind, lapply(out, `[[`, "counts"))

  # Hand-made counts for two markers, A and B, in one donor and two rules.
  small_counts <- function() {
    combn <- c("A-B-", "A+B-", "A-B+", "A+B+")
    rbind(
      data.frame(pop = "p", stim = "s", donor = "d1", rule = "ordinary", tube = "stim",
        combn = combn, n = c(900L, 50L, 30L, 20L)),
      data.frame(pop = "p", stim = "s", donor = "d1", rule = "ordinary", tube = "uns",
        combn = combn, n = c(1980L, 10L, 10L, 0L)),
      data.frame(pop = "p", stim = "s", donor = "d1", rule = "unscreened", tube = "stim",
        combn = combn, n = c(900L, 40L, 30L, 30L)),
      data.frame(pop = "p", stim = "s", donor = "d1", rule = "unscreened", tube = "uns",
        combn = combn, n = c(1980L, 8L, 10L, 2L))
    )
  }

  testthat::test_that("tube results cover every pair, variant, rule and combination", {
    testthat::expect_equal(nrow(low), 2L * 2L * 6L)
    testthat::expect_setequal(unique(low$rule), c("screened", "unscreened"))
    sel <- low$donor == "d1" & low$a == "IFNg" & low$b == "TNF"
    testthat::expect_true(all(low$lowered[sel]))
    testthat::expect_false(any(low$lowered[low$donor == "d2" & low$rule == "screened"]))
    # A cytokine without a gate is never lowered nor used to lower.
    testthat::expect_false(any(low$lowered[low$a == "IL2" | low$b == "IL2"]))
    # Without the screen every pair the screen passes gives the same gates.
    scr <- low[low$rule == "screened", ]
    uns <- low[low$rule == "unscreened", ]
    m <- merge(scr, uns, by = c("donor", "a", "b"), suffixes = c(".s", ".u"))
    pass <- m$lowered.s
    testthat::expect_true(all(m$lowered.u[pass]))
    testthat::expect_equal(m$cut.s[pass], m$cut.u[pass])
    testthat::expect_equal(m$condCut.s[pass], m$condCut.u[pass])
    # Cells added are new positives for b, so only where a gate was lowered.
    testthat::expect_true(all(low$addedStim[!low$lowered] == 0L))
    testthat::expect_gt(low$addedStim[sel & low$rule == "screened"], 0L)
    testthat::expect_setequal(unique(counts$rule), c("ordinary", "screened", "unscreened"))
    tot <- stats::aggregate(n ~ donor + rule + tube, counts, sum)
    testthat::expect_true(all(tot$n[tot$tube == "uns"] == 20000L))
    # IL2 has no gate but stays in the combinations, so every tube shares them.
    testthat::expect_equal(length(unique(counts$combn)), 8L)
    testthat::expect_equal(sum(counts$n[grepl("IL2+", counts$combn, fixed = TRUE)]), 0)
  })

  testthat::test_that("single- and multi-positive frequencies move between rules", {
    paired <- env$.acsCoexGatesPaired(env$.acsCoexGatesSingleMulti(env$.acsCoexGatesFreq(counts)))
    testthat::expect_setequal(unique(paired$rule), c("screened", "unscreened"))
    d1 <- paired[paired$donor == "d1" & paired$cytokine == "any" & paired$rule == "screened", ]
    testthat::expect_gt(d1$diffBs[d1$type == "multi"], 0)
    testthat::expect_lt(d1$diffBs[d1$type == "single"], 0)
    d2 <- paired[paired$donor == "d2" & paired$rule == "screened", ]
    testthat::expect_equal(d2$diffBs, rep(0, nrow(d2)))
    flags <- env$.acsCoexGatesFlags(low, paired)
    testthat::expect_equal(nrow(flags), 4L)
    testthat::expect_setequal(paste(flags$donor, flags$rule),
      c("d1 screened", "d1 unscreened", "d2 screened", "d2 unscreened"))
    testthat::expect_false(any(flags$flagControl[flags$rule == "screened"]))
  })

  testthat::test_that("combination labels parse into a positivity matrix", {
    pos <- env$.acsCoexGatesCombnMatrix(c("IL2-IL22+", "IL2+IL22-", "IL2+IL22+"))
    testthat::expect_equal(colnames(pos), c("IL2", "IL22"))
    testthat::expect_equal(unname(pos[, "IL2"]), c(FALSE, TRUE, TRUE))
    testthat::expect_equal(unname(pos[, "IL22"]), c(TRUE, FALSE, TRUE))
    testthat::expect_error(env$.acsCoexGatesCombnMatrix(c("A+B-", "B+A-")), "same markers")
    testthat::expect_equal(env$.acsCoexGatesCombnLabel(c("A+B-C+", "A+"), c("exclusive", "marginal")),
      c("A+C+", "A+ (any)"))
  })

  testthat::test_that("test counts give exclusive combinations and marginals", {
    tc <- env$.acsCoexGatesTestCounts(small_counts())
    o <- tc[tc$rule == "ordinary", ]
    testthat::expect_false("A-B-" %in% o$combn)
    testthat::expect_setequal(o$combn[o$type == "exclusive"], c("A+B-", "A-B+", "A+B+"))
    testthat::expect_equal(o$nStim[o$combn == "A+"], 70)
    testthat::expect_equal(o$nUns[o$combn == "B+"], 10)
    testthat::expect_true(all(o$NStim == 1000) && all(o$NUns == 2000))
    testthat::expect_equal(o$propBs[o$combn == "A+B+"], 0.02)
  })

  testthat::test_that("Fisher, Wilcoxon and t-test helpers match base R", {
    cases <- list(c(20, 1000, 0, 2000), c(0, 1000, 0, 2000), c(0, 1000, 5, 2000), c(3, 50, 3, 50))
    for (x in cases) {
      tab <- matrix(c(x[[1]], x[[3]], x[[2]] - x[[1]], x[[4]] - x[[3]]), 2)
      testthat::expect_equal(env$.acsCoexGatesFisherGreater(x[[1]], x[[2]], x[[3]], x[[4]]),
        stats::fisher.test(tab, alternative = "greater")$p.value)
    }
    # Zero and tied differences: dropped and normal approximation, as wilcox.test().
    s <- c(0.010, 0.020, 0.015, 0.010, 0.030)
    u <- c(0.010, 0.010, 0.005, 0.000, 0.020)
    ref <- suppressWarnings(stats::wilcox.test(s, u, paired = TRUE, alternative = "greater")$p.value)
    testthat::expect_equal(env$.acsCoexGatesWilcoxGreater(s - u), ref)
    testthat::expect_true(is.na(env$.acsCoexGatesWilcoxGreater(c(0, 0, 0))))
    testthat::expect_equal(env$.acsCoexGatesTGreater(s - u),
      stats::t.test(s, u, paired = TRUE, alternative = "greater")$p.value)
    testthat::expect_true(is.na(env$.acsCoexGatesTGreater(c(0.01, 0.01, 0.01))))
    testthat::expect_true(is.na(env$.acsCoexGatesTGreater(c(0, 0))))
    testthat::expect_true(is.na(env$.acsCoexGatesTGreater(0.01)))
  })

  testthat::test_that("response tests adjust within each population, stimulation and rule", {
    withr::with_seed(3, {
      tc <- do.call(rbind, lapply(c("ordinary", "unscreened"), function(rule) {
        do.call(rbind, lapply(paste0("c", 1:4), function(combn) {
          shift <- c(c1 = 30, c2 = 5, c3 = 0, c4 = 2)[[combn]] + (rule == "unscreened") * 5
          nU <- rpois(8, 10)
          data.frame(pop = "p", stim = "s", donor = paste0("d", 1:8), rule = rule,
            type = "exclusive", combn = combn, nStim = nU + rpois(8, shift), NStim = 10000,
            nUns = nU, NUns = 10000)
        }))
      }))
    })
    # One combination with no difference at all: untestable.
    tc$nStim[tc$combn == "c3"] <- tc$nUns[tc$combn == "c3"]
    tc$propStim <- tc$nStim / tc$NStim
    tc$propUns <- tc$nUns / tc$NUns
    tc$propBs <- tc$propStim - tc$propUns
    tests <- env$.acsCoexGatesResponseTests(tc)
    testthat::expect_equal(nrow(tests), 8L)
    testthat::expect_true(all(is.na(tests$pWilcox[tests$combn == "c3"])))
    testthat::expect_false(any(tests$responding[tests$combn == "c3"]))
    for (rule in c("ordinary", "unscreened")) {
      d <- tests[tests$rule == rule, ]
      testthat::expect_equal(d$pWilcoxAdj, stats::p.adjust(d$pWilcox, "BH"))
      testthat::expect_equal(d$pTAdj, stats::p.adjust(d$pT, "BH"))
    }
    # Another rule's (or population's) tests do not change a family's adjustment.
    alone <- env$.acsCoexGatesResponseTests(tc[tc$rule == "ordinary", ])
    testthat::expect_equal(alone$pWilcoxAdj, tests$pWilcoxAdj[tests$rule == "ordinary"])
    other_pop <- transform(tc[tc$rule == "ordinary", ], pop = "q", nStim = nUns + 1L)
    other_pop$propBs <- other_pop$nStim / other_pop$NStim - other_pop$propUns
    both <- env$.acsCoexGatesResponseTests(rbind(tc[tc$rule == "ordinary", ], other_pop))
    testthat::expect_equal(both$pWilcoxAdj[both$pop == "p"], alone$pWilcoxAdj)
    p1 <- tc[tc$rule == "ordinary" & tc$combn == "c1", ]
    testthat::expect_equal(tests$nDonorFisher[tests$rule == "ordinary" & tests$combn == "c1"],
      sum(env$.acsCoexGatesFisherGreater(p1$nStim, p1$NStim, p1$nUns, p1$NUns) < 0.05))
    testthat::expect_true(tests$responding[tests$rule == "ordinary" & tests$combn == "c1"])

    # Compare and summarise: a combination only the unscreened rule calls.
    tests3 <- rbind(tests, transform(tests[tests$rule == "ordinary", ], rule = "screened"))
    tests3$responding[tests3$combn == "c4"] <- tests3$rule[tests3$combn == "c4"] == "unscreened"
    cmp <- env$.acsCoexGatesResponseCompare(tests3)
    testthat::expect_true(cmp$unscreenedOnly[cmp$combn == "c4"])
    testthat::expect_true(cmp$changed[cmp$combn == "c4"])
    testthat::expect_false(cmp$unscreenedOnly[cmp$combn == "c3"])
    sm <- env$.acsCoexGatesResponseSummary(tests3, cmp)
    testthat::expect_equal(nrow(sm), 3L)
    testthat::expect_true(all(sm$nCombn == 4L))
    testthat::expect_equal(sm$nTested, rep(3L, 3L))
    testthat::expect_gte(sm$nOnlyThisRule[sm$rule == "unscreened"], 1L)
    testthat::expect_s3_class(env$.acsCoexGatesResponseHeatmap(tests3), "ggplot")
  })

  testthat::test_that("control-tube additions are counted against the ordinary gates", {
    ct <- env$.acsCoexGatesContamination(env$.acsCoexGatesTestCounts(small_counts()))
    testthat::expect_equal(nrow(ct), 1L)
    testthat::expect_equal(ct$rule, "unscreened")
    # A+: 70 -> 70 stimulated, 10 -> 10 control; B+: 50 -> 60 and 10 -> 12.
    testthat::expect_equal(ct$addedStim, 10)
    testthat::expect_equal(ct$addedUns, 2)
    testthat::expect_equal(ct$addedUnsPerMillion, 1000)
    testthat::expect_equal(ct$controlShare, 0.1)
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
    for (rule in c("screened", "unscreened")) {
      kept <- env$.acsCoexGatesCompassInput(counts, rule, keep)
      testthat::expect_equal(colnames(kept$n_s)[[ncol(kept$n_s)]], "!IFNg&!TNF&!IL2")
    }
  })

  testthat::test_that("Analysis 18 plots build without titles", {
    paired <- env$.acsCoexGatesPaired(env$.acsCoexGatesSingleMulti(env$.acsCoexGatesFreq(counts)))
    scores <- data.frame(pop = "cd4", stim = "p1", donor = rep(c("d1", "d2"), 3),
      rule = rep(c("ordinary", "screened", "unscreened"), each = 2),
      fs = c(0.1, 0, 0.2, 0, 0.25, 0.01), pfs = c(0.1, 0, 0.15, 0, 0.2, 0.01))
    dat <- tube("d1", TRUE)
    lw <- low[low$donor == "d1" & low$rule == "screened" & low$a == "IFNg" & low$b == "TNF", ]
    tests <- env$.acsCoexGatesResponseTests(env$.acsCoexGatesTestCounts(counts))
    plots <- list(
      env$.acsCoexGatesScatterPlot(paired), env$.acsCoexGatesDiffPlot(paired),
      env$.acsCoexGatesCompassPlot(scores, "fs"), env$.acsCoexGatesCompassPlot(scores, "pfs"),
      env$.acsCoexGatesResponseHeatmap(tests),
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
