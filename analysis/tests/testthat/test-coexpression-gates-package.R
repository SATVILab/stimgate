# The package's co-expression gates (stimControl(cytPosMethod =
# "coexpression")) must reproduce the analysis reference implementation in
# scripts/r/coexpression-gates.R, on data where gates are lowered.
local({
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "scripts", "r", "coexpression-gates.R"), env)

  testthat::test_that("package co-expression gates and counts match the reference", {
    withr::local_seed(11)
    neg <- function(n) cbind(A = abs(rnorm(n, 0.6, 0.4)), B = abs(rnorm(n, 0.8, 0.5)),
      C = abs(rnorm(n, 0.7, 0.4)))
    resp <- function(n) cbind(A = rnorm(n, 4, 0.4), B = runif(n, 1.8, 5.5), C = abs(rnorm(n, 0.7, 0.4)))
    tubes <- list(neg(6000), rbind(neg(6000), resp(250)), neg(6000), rbind(neg(6000), resp(200)))
    names(tubes) <- c("u1", "s1", "u2", "s2")
    path <- withr::local_tempdir("coexpkg")
    stimgate::gateStim(path, tubes, batchList = list(c(1L, 2L), c(3L, 4L)),
      chnl = c("A", "B", "C"),
      control = stimgate::stimControl(cytPosMethod = "coexpression", clusterGates = FALSE))
    low <- stimgate::getStimGatesCoexpression(path)
    gates <- stimgate::getStimGates(path)
    stats <- stimgate::getStimStats(path)
    testthat::expect_true(any(low$lowered))
    for (b in list(c(1L, 2L), c(3L, 4L))) {
      ind <- b[[2]]
      g <- gates[gates$ind == as.character(ind), ]
      gate <- stats::setNames(g$gate, g$chnl)[c("A", "B", "C")]
      ex <- function(i) as.data.frame(stimgate::getStimExpr(path, ind = i, chnl = c("A", "B", "C")))
      dat <- list(stim = ex(ind), uns = ex(b[[1]]), gate = gate)
      ref <- env$.coexLowerGates(dat)
      pkg <- low[low$ind == as.character(ind), ]
      pkg <- pkg[match(paste(ref$a, ref$b), paste(pkg$chnlCond, pkg$chnl)), ]
      testthat::expect_equal(pkg$cut, ref$cut)
      testthat::expect_equal(pkg$condCut, ref$condCut)
      testthat::expect_equal(pkg$lowered, ref$lowered)
      testthat::expect_equal(pkg$z, ref$z)
      rows <- stats[stats$ind == as.character(ind), ]
      for (tube in c("stim", "uns")) {
        n <- env$.coexCombnCounts(env$.coexPositive(dat[[tube]], gate, ref))
        lab <- gsub("([+-])", "~\\1~", names(n))
        got <- if (tube == "stim") rows$countStim else rows$countUns
        testthat::expect_equal(unname(n[match(rows$cytCombn, lab)]), got)
      }
    }
  })
})
