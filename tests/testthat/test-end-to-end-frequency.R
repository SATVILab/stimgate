# End-to-end checks of the default pipeline against simulated truth:
# negatives N(0, 1), responders N(5, 0.6), well separated, in tubes of
# 20,000 cells. Unstimulated tubes have no responders, so the true net
# frequency equals the stimulated response fraction.

.e2eSim <- function(seed, fracs, n = 20000L) {
  withr::with_seed(seed, {
    mk <- function(f) {
      nResp <- stats::rbinom(1L, n, f)
      x <- c(stats::rnorm(n - nResp, 0, 1), stats::rnorm(nResp, 5, 0.6))
      list(
        x = matrix(x, ncol = 1L, dimnames = list(NULL, "cyt")),
        resp = rep(c(FALSE, TRUE), c(n - nResp, nResp))
      )
    }
    tubes <- list()
    for (i in seq_along(fracs)) {
      tubes[[paste0("uns", i)]] <- mk(0)
      tubes[[paste0("stim", i)]] <- mk(fracs[[i]])
    }
    tubes
  })
}

# Gate with the default settings; return, per stimulated tube, the applied
# gate, the reported net frequency and label-based counts at that gate.
.e2eGate <- function(tubes) {
  pathProject <- withr::local_tempdir("e2e_freq")
  nBatch <- length(tubes) / 2
  withr::with_seed(1L, suppressMessages(gateStim(
    .data = lapply(tubes, `[[`, "x"), pathProject = pathProject,
    batchList = lapply(seq_len(nBatch), function(i) c(2 * i - 1, 2 * i)),
    chnl = "cyt", calcCytPosGates = FALSE
  )))
  gates <- getStimGates(pathProject)
  gates <- gates[gates$gateName == "loc_minClust", , drop = FALSE]
  stats <- getStimStats(pathProject)
  stats <- stats[
    stats$gateName == "loc_minClust" & stats$cytCombn == "cyt~+~", ,
    drop = FALSE
  ]
  indStim <- 2L * seq_len(nBatch)
  out <- lapply(indStim, function(ind) {
    gate <- gates$gate[as.numeric(gates$ind) == ind]
    stim <- tubes[[ind]]
    uns <- tubes[[ind - 1L]]
    pos <- stim$x[, 1] > gate
    data.frame(
      ind = ind,
      gate = gate,
      truth = mean(stim$resp),
      est = stats$propBs[as.numeric(stats$ind) == ind],
      sens = if (any(stim$resp)) mean(pos[stim$resp]) else NA_real_,
      fpStim = mean(pos[!stim$resp]),
      fpUns = mean(uns$x[, 1] > gate)
    )
  })
  do.call(rbind, out)
}

test_that("default gating recovers known response frequencies", {
  fracs <- c(0.01, 0.1, 0.002, 0.05, 0.005, 0.02)
  res <- .e2eGate(.e2eSim(1L, fracs))
  expect_identical(nrow(res), length(fracs))
  # The reported net frequency is the label-free estimate of the truth.
  expect_true(all(abs(res$est - res$truth) <= 0.3 * res$truth + 5e-4))
  # Gates sit between the populations: most responders are above the gate
  # and almost no negatives in either tube are.
  expect_true(all(res$sens >= 0.9))
  expect_true(all(res$fpStim <= 0.002))
  expect_true(all(res$fpUns <= 0.002))
})

test_that("null tubes gated beside responders give near-zero frequencies", {
  for (seed in 2:3) {
    fracs <- c(0, 0, 0, 0.01, 0.02, 0.05)
    res <- .e2eGate(.e2eSim(seed, fracs))
    null <- res[res$truth == 0, , drop = FALSE]
    resp <- res[res$truth > 0, , drop = FALSE]
    expect_identical(nrow(null), 3L)
    # Gate sharing may lower a null tube's gate, but not into the negatives.
    expect_true(all(abs(null$est) <= 1e-3), info = paste("seed", seed))
    expect_true(all(null$fpStim <= 0.002), info = paste("seed", seed))
    expect_true(
      all(abs(resp$est - resp$truth) <= 0.3 * resp$truth + 5e-4),
      info = paste("seed", seed)
    )
  }
})
