local({
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "scripts", "r", "analysis-plot-style.R"), env)
  sys.source(file.path(root, "scripts", "r", "acs_cytof-explore-score.R"), env)

  testthat::test_that("tail p-values are conformal and zeros score nothing", {
    ref <- c(0, 0, 1, 2, 3)
    # New cells: (#ref >= x + 1) / (n + 1).
    testthat::expect_equal(env$.acsScoreSurv(ref, c(2, 3.5)), c(3 / 6, 1 / 6))
    # Reference cells count themselves: #ref >= x / n.
    testthat::expect_equal(env$.acsScoreSurv(ref, c(1, 3), self = TRUE), c(3 / 5, 1 / 5))
    testthat::expect_equal(env$.acsScoreSurv(ref, c(0, -1)), c(1, 1))
  })

  testthat::test_that("scores add over markers and tails and nets are fractions", {
    dat <- list(
      stim = data.frame(A = c(0, 3, 3), B = c(0, 3, 0)),
      uns = data.frame(A = c(0, 1, 2), B = c(0, 1, 2)),
      gate = c(A = 2.5, B = 2.5)
    )
    tbl <- env$.acsScoreTbl(dat, c("A", "B"), "d1")
    stim <- tbl$score[tbl$tube == "stimulated"]
    testthat::expect_equal(stim, c(0, 2 * -log10(1 / 4), -log10(1 / 4)))
    tail <- env$.acsScoreTail(tbl, grid = c(0, 1))
    at1 <- tail$net[tail$net$threshold == 1, "net"]
    testthat::expect_equal(at1, 1 / 3 - 0)
    testthat::expect_equal(env$.acsScoreRectNet(dat, c("A", "B"), "d1")$net, 1 / 3)
    testthat::expect_true(env$.acsScoreKl(tbl)$klBits >= 0)
  })
})

local({
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "scripts", "r", "analysis-plot-style.R"), env)
  sys.source(file.path(root, "scripts", "r", "acs_cytof-explore-score.R"), env)

  testthat::test_that("residuals are zero under exact independence", {
    bins <- data.frame(donor = "d", tube = "stimulated",
      x = rep(c(1, 2), 2), y = rep(c(1, 2), each = 2), count = c(4, 2, 2, 1))
    res <- env$.acsCoexResiduals(bins)
    testthat::expect_equal(res$residual, rep(0, 4))
  })

  testthat::test_that("co-expression beyond independence is found only with a joint response", {
    withr::local_seed(1)
    neg <- function(n) data.frame(IFNg = abs(rnorm(n, 0.3, 0.3)), TNF = abs(rnorm(n, 0.3, 0.3)))
    resp <- data.frame(IFNg = rnorm(150, 3.5, 0.4), TNF = rnorm(150, 3.5, 0.4))
    gate <- c(IFNg = 1.5, TNF = 2)
    dat <- list(stim = rbind(neg(3000), resp), uns = neg(3000), gate = gate)
    ind <- env$.acsCoexIndependence(dat, "IFNg", "TNF", "d1")
    s <- env$.acsCoexIndependenceSummary(ind, dat, "IFNg", "TNF")
    testthat::expect_gt(s$stimCalled, 100)
    testthat::expect_lt(s$controlCalled, 5)
    testthat::expect_gt(s$netCalledPct, 3)
    none <- list(stim = neg(3000), uns = neg(3000), gate = gate)
    s0 <- env$.acsCoexIndependenceSummary(
      env$.acsCoexIndependence(none, "IFNg", "TNF", "d2"), none, "IFNg", "TNF")
    testthat::expect_lte(s0$stimCalled, 5)
  })
})

local({
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "scripts", "r", "analysis-plot-style.R"), env)
  sys.source(file.path(root, "scripts", "r", "acs_cytof-explore-score.R"), env)

  testthat::test_that("the floor sits about 1.5 SD above a normal negative peak", {
    withr::local_seed(2)
    x <- rnorm(20000, 1, 0.4)
    testthat::expect_equal(env$.acsCoexNegFloor(x), 1 + 1.5 * 0.4, tolerance = 0.08)
  })

  testthat::test_that("lowered gates follow a stimulation-specific diagonal", {
    withr::local_seed(3)
    neg <- function(n) data.frame(IFNg = abs(rnorm(n, 0.4, 0.3)), TNF = abs(rnorm(n, 0.6, 0.4)))
    # Responders: IFNg high, TNF spread from 2 to 5 (mostly below TNF's gate of 4.5).
    resp <- data.frame(IFNg = rnorm(200, 3.5, 0.4), TNF = runif(200, 2, 5))
    gate <- c(IFNg = 2, TNF = 4.5)
    dat <- list(stim = rbind(neg(20000), resp), uns = neg(20000), gate = gate)
    low <- env$.acsCoexLowerGate(dat, a = "IFNg", b = "TNF")
    testthat::expect_gt(low$z, 2)
    testthat::expect_lt(low$cut, 2.5)
    testthat::expect_gte(low$cut, env$.acsCoexNegFloor(dat$uns$TNF))
    testthat::expect_gte(low$condCut, gate[["IFNg"]])
    s <- env$.acsCoexLowerSummary(dat, "IFNg", "TNF", "d1")$summary
    testthat::expect_gt(s$netPct, s$netRectanglePct)
    testthat::expect_lte(s$controlCalled, 2L)
    # The same co-expression in the control tube: no double-positive response,
    # so nothing moves.
    both <- list(stim = rbind(neg(20000), resp), uns = rbind(neg(20000), resp), gate = gate)
    testthat::expect_equal(env$.acsCoexLowerGate(both, a = "IFNg", b = "TNF")$cut, 4.5)
  })

  testthat::test_that("purity is one with no control cells and undefined with none at all", {
    testthat::expect_equal(env$.acsCoexPurity(c(4, 4, 0, 0), c(0, 2, 1, 0), 100, 100),
      c(1, 0.5, -Inf, NA))
  })

  testthat::test_that("a background band just above the conditioning gate is trimmed", {
    withr::local_seed(4)
    neg <- function(n) data.frame(IFNg = abs(rnorm(n, 0.4, 0.3)), TNF = abs(rnorm(n, 0.6, 0.4)))
    resp <- data.frame(IFNg = rnorm(200, 4, 0.3), TNF = runif(200, 2, 5))
    # Background in both tubes: IFNg just above its gate with mid TNF.
    bg <- function() data.frame(IFNg = runif(40, 2, 2.6), TNF = runif(40, 2, 4))
    gate <- c(IFNg = 2, TNF = 4.5)
    dat <- list(stim = rbind(neg(20000), resp, bg()), uns = rbind(neg(20000), bg()), gate = gate)
    low <- env$.acsCoexLowerGate(dat, a = "IFNg", b = "TNF")
    testthat::expect_lt(low$cut, 4.5)
    testthat::expect_gt(low$condCut, 2.2)
  })
})
