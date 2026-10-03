test_that("empirical cluster quantiles preserve the lower endpoint exactly", {
  withr::local_preserve_seed()
  set.seed(1)
  for (n in c(1L, 2L, 3L, 10L, 99L)) {
    x <- rnorm(n)
    taus <- c(0, 0.15, 0.6, 0.85, 1, seq_len(n) / n)
    expect_identical(
      .getCpClusterLocRqQuantile(x, taus),
      sort(x)[pmax(1L, ceiling(n * taus))]
    )
    for (tau in taus) {
      expected <- sort(x)[[max(1L, ceiling(n * tau))]]
      expect_identical(.getCpClusterLocRqQuantile(x, tau), expected)
    }
  }
  expect_identical(.getCpClusterLocRqQuantile(c(NA, Inf), 0.6), NA_real_)
  expect_identical(.getCpClusterLocRqQuantile(1:3, Inf), NA_real_)
  expect_identical(.getCpClusterLocRqQuantile(c("1", "2", NA), 0), 1)
  expect_identical(
    .getCpClusterLocRqQuantile(1:3, c(NA, 0.6, Inf)),
    c(NA_real_, 2, NA_real_)
  )
})

test_that("skip outputs preserve row order, types, gates and donor metadata", {
  gates <- .getCpClusterLocGateTblPrepare(tibble::tibble(
    ind = c("B", "A"), gate = c(Inf, 2),
    locGenerated = c(FALSE, TRUE), locGeneratedDirect = c(FALSE, TRUE),
    locSource = c("fallback", NA_character_),
    locReason = c("unavailable", NA_character_)
  ))
  out <- .getCpClusterLocSkipOut(gates, "no_donors")
  expect_identical(out$ind, gates$ind)
  expect_identical(out$cpJoinTgOrig, gates$gate)
  expect_identical(out$locSource, gates$locSource)
  expect_identical(out$locReason, gates$locReason)
  expect_identical(out$locGeneratedDirect, gates$locGeneratedDirect)
  expect_identical(out$locClusterReason, rep("no_donors", 2))
  expect_identical(out$locClusterNDirect, rep(NA_integer_, 2))
  expect_identical(
    .getCpClusterLocSkipOut(gates[0, ], "no_donors"), out[0, ]
  )
})

test_that("cluster quantiles preserve winsorisation and small donor groups", {
  gates <- .getCpClusterLocGateTblPrepare(tibble::tibble(
    ind = as.character(seq_len(15)),
    gate = c(seq_len(10), Inf, 20, 30, Inf, Inf),
    locGeneratedDirect = c(rep(TRUE, 10), FALSE, TRUE, TRUE, FALSE, FALSE),
    grp = c(rep("A", 11), rep("B", 3), "C")
  ))
  out <- .getCpClusterLocApplyQuantiles(
    gates,
    commonBw = 0.1,
    control = .getCpClusterControlUpdate(list()), nInitialClusters = 3L
  )
  expect_identical(
    out$cpJoinTgOrig,
    c(2, 2, 3, 4, 5, 6, 7, 8, 9, 9, 6, 20, 30, 30, Inf)
  )
  expect_identical(out$locClusterQ15, c(rep(2, 11), rep(NA_real_, 4)))
  expect_identical(out$locClusterQ60, c(rep(6, 11), rep(30, 3), NA_real_))
  expect_identical(out$locClusterQ85, c(rep(9, 11), rep(NA_real_, 4)))
  control <- .getCpClusterControlUpdate(list(winsorLower = c(0.15, 0.2)))
  expect_identical(
    .getCpClusterLocApplyQuantiles(gates, 0.1, control, 3L), out
  )
})
