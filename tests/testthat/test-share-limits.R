# Gate sharing: only responders donate, and lowered gates are limited.

.limit_ex <- function(x) {
  df <- data.frame(CD4 = x)
  attr(df, "chnlCut") <- "CD4"
  df
}

.limit_freq <- function(gate, x_stim, x_uns) {
  mean(x_stim > gate) - mean(x_uns > gate)
}

test_that("a tube with negative own frequency does not donate its low gate", {
  bg <- seq(0, 0.9, length.out = 60)
  ex_list <- list(
    uns = .limit_ex(c(bg, rep(1.5, 40))),
    stim1 = .limit_ex(c(bg, rep(1.5, 30), rep(4, 10))),
    stim2 = .limit_ex(c(bg, rep(1.5, 30), rep(2.7, 5), rep(4, 5))),
    # model-generated low gate, but fewer stim than unstim cells above it
    stim3 = .limit_ex(c(seq(0, 0.9, length.out = 70), rep(1.5, 30)))
  )
  cp <- c(stim1 = 3, stim2 = 2.5, stim3 = 1, uns = 2.2)
  attr(cp, "locGenerated") <- rep(TRUE, 4)
  attr(cp, "locGeneratedDirect") <- c(TRUE, TRUE, TRUE, FALSE)
  attr(cp, "locSource") <- c(rep("direct", 3), "unstim_summary")
  attr(cp, "propBsEst") <- c(0.1, 0.1, 0.05, NA)

  for (caps in list(c(Inf, Inf), c(1.5, 0.5))) {
    out <- .getCpUnsLocCombineCpWithMeta(
      cp, c("min", "no"),
      exListOrig = ex_list, shareCap = caps[[1]], cellCap = caps[[2]]
    )
    cp_min <- out[["min"]]
    expect_equal(unname(cp_min[c("stim1", "stim2", "stim3")]), rep(2.5, 3))
    meta <- .getCpUnsLocMetaFromCp(cp_min)
    expect_identical(meta$locResponder, c(TRUE, TRUE, FALSE, FALSE))
    expect_identical(meta$locShareLimit, rep("none", 4))
    expect_identical(meta$locShareProposed[1:3], rep(2.5, 3))
    expect_identical(meta$locGeneratedDirect, c(FALSE, TRUE, FALSE, FALSE))
    expect_identical(meta$locSource[3], "combined")
    # responder status is carried by every combination
    expect_identical(
      attr(out[["no"]], "locResponder"), c(TRUE, TRUE, FALSE, FALSE)
    )
  }

  # median of responders only (2.5 and 3), not of all generated gates
  out <- .getCpUnsLocCombineCpWithMeta(
    cp, "median",
    exListOrig = ex_list, shareCap = Inf, cellCap = Inf
  )
  expect_equal(unname(out[["median"]]["stim3"]), 2.75)
})

test_that("generated gates without a responder are kept", {
  ex_list <- list(
    uns = .limit_ex(c(1, 2, 3)),
    stim1 = .limit_ex(c(1, 1, 1)),
    stim2 = .limit_ex(c(1, 1, 3))
  )
  cp <- c(stim1 = 0.5, stim2 = 2.5, uns = 1.5)
  attr(cp, "locGenerated") <- rep(TRUE, 3)
  attr(cp, "locGeneratedDirect") <- c(TRUE, TRUE, FALSE)
  attr(cp, "locSource") <- c("direct", "direct", "unstim_summary")
  out <- .getCpUnsLocCombineCpWithMeta(cp, "min", exListOrig = ex_list)
  expect_equal(unname(as.numeric(out[["min"]])), c(0.5, 2.5, 1.5))
  expect_false(any(attr(out[["min"]], "locResponder")))
})

test_that("a responder accepts a lower gate only within locShareCap", {
  x_uns <- seq(0.01, 1, length.out = 100)
  x_stim <- c(
    seq(0.01, 1, length.out = 80), rep(1.3, 5), rep(1.6, 5), rep(3, 10)
  )
  args <- list(
    gs = 1.2, gc = 2.5, responder = TRUE, propBsEst = 0.1,
    freqDonor = 0.2, xStim = x_stim, xUns = x_uns, cellCap = 0.5
  )
  res <- do.call(.getCpShareApply, c(args, shareCap = 1.5))
  # cells at 1.6 and above give 0.15 = 1.5 * 0.1; the gate sits below 1.6
  expect_equal(res$gate, 1.45)
  expect_identical(res$limit, "responder_cap")
  expect_lte(.limit_freq(res$gate, x_stim, x_uns), 1.5 * 0.1 + 1e-12)
  expect_gt(.limit_freq(1.25, x_stim, x_uns), 1.5 * 0.1)

  res_inf <- do.call(.getCpShareApply, c(args, shareCap = Inf))
  expect_identical(res_inf, list(gate = 1.2, limit = "none"))

  # never above the current gate
  res_gc <- do.call(
    .getCpShareApply, c(utils::modifyList(args, list(gc = 1.4)), shareCap = 1.5)
  )
  expect_identical(res_gc, list(gate = 1.4, limit = "responder_cap"))

  # no estimate: keep the current gate
  res_na <- do.call(
    .getCpShareApply,
    c(utils::modifyList(args, list(propBsEst = NA_real_)), shareCap = 1.5)
  )
  expect_identical(res_na, list(gate = 2.5, limit = "responder_no_estimate"))

  # a higher shared gate is accepted
  res_up <- do.call(
    .getCpShareApply, c(utils::modifyList(args, list(gs = 2.8)), shareCap = 1.5)
  )
  expect_identical(res_up, list(gate = 2.8, limit = "none"))
})

test_that("a non-responder is limited by the cell and donor-median caps", {
  x_uns <- seq(0.01, 1, length.out = 100)
  x_stim <- c(
    seq(0.01, 1, length.out = 80), rep(1.3, 10), rep(1.6, 5), rep(3, 2),
    rep(3.5, 3)
  )
  args <- list(
    gs = 1.2, gc = Inf, responder = FALSE, propBsEst = NA_real_,
    xStim = x_stim, xUns = x_uns, shareCap = 1.5
  )
  # cell cap 10 / 100 = 0.1 binds below a donor median of 0.5
  res_cell <- do.call(
    .getCpShareApply, c(args, freqDonor = 0.5, cellCap = 10)
  )
  expect_equal(res_cell$gate, 1.45)
  expect_identical(res_cell$limit, "nonresponder_cell_cap")

  # donor median 0.03 binds below the cell cap of 0.1
  res_donor <- do.call(
    .getCpShareApply, c(args, freqDonor = 0.03, cellCap = 10)
  )
  expect_equal(res_donor$gate, 3.25)
  expect_identical(res_donor$limit, "nonresponder_donor_cap")
  expect_lte(.limit_freq(res_donor$gate, x_stim, x_uns), 0.03 + 1e-12)

  # Inf cell cap switches the rule off
  res_inf <- do.call(
    .getCpShareApply, c(args, freqDonor = 0.03, cellCap = Inf)
  )
  expect_identical(res_inf, list(gate = 1.2, limit = "none"))
})

test_that("default limits never give a non-responder half a net cell", {
  withr::local_seed(3)
  for (i in seq_len(20)) {
    n_stim <- sample(c(50, 500, 5000), 1)
    x_uns <- stats::rnorm(sample(c(50, 500, 5000), 1))
    x_stim <- c(
      stats::rnorm(n_stim),
      stats::rnorm(stats::rbinom(1, n_stim, 0.05), mean = 3)
    )
    gs <- stats::runif(1, -1, 3)
    res <- .getCpShareApply(
      gs = gs, gc = stats::runif(1, -1, 5), responder = FALSE,
      propBsEst = NA_real_, freqDonor = 0.2, xStim = x_stim, xUns = x_uns,
      shareCap = 1.5, cellCap = 0.5
    )
    expect_gte(res$gate, gs)
    expect_lte(
      length(x_stim) * .limit_freq(res$gate, x_stim, x_uns), 0.5 + 1e-9
    )
  }
})

test_that("cluster donors are responders and lowered gates are limited", {
  bg <- seq(0, 0.5, length.out = 900)
  x_uns <- seq(0, 0.5, length.out = 1000)
  gates <- c(1:7, Inf, 10)
  stim <- lapply(gates[1:6], function(g) c(bg, rep(g + 0.5, 100)))
  stim[[7]] <- c(bg[1:890], rep(6.5, 30), rep(6.8, 70), rep(9, 10))
  stim[[8]] <- stim[[9]] <- c(bg, rep(5.2, 100))
  ind <- LETTERS[seq_along(gates)]
  ex_lookup <- stats::setNames(lapply(seq_along(ind), function(i) {
    list(ind = ind[[i]], stim = .limit_ex(stim[[i]]), uns = .limit_ex(x_uns))
  }), ind)

  tbl <- .getCpClusterLocGateTblPrepare(tibble::tibble(
    ind = ind,
    gate = gates,
    locGenerated = c(rep(TRUE, 7), FALSE, TRUE),
    # B lost direct status in the batch step but is still a responder
    locGeneratedDirect = c(TRUE, FALSE, rep(TRUE, 5), FALSE, FALSE),
    locResponder = c(rep(TRUE, 7), FALSE, FALSE),
    propBsEst = c(rep(0.1, 6), 0.06, NA, NA),
    grp = "1"
  ))
  proposed <- .getCpClusterLocApplyQuantiles(
    tbl,
    commonBw = 0.1, control = .getCpClusterControlUpdate(list()),
    nInitialClusters = 1L
  )
  expect_identical(proposed$locClusterNDirect, rep(7L, 9))
  expect_identical(proposed$cpJoinTgOrig, c(2, 2:6, 6, 5, 5))
  expect_identical(proposed$locGeneratedDirect[1:2], c(TRUE, FALSE))
  expect_identical(proposed$locResponder, tbl$locResponder)

  # Inf limits leave the proposed gates
  unlimited <- .getCpClusterLocLimit(proposed, ex_lookup, Inf, Inf)
  expect_identical(unlimited, proposed)

  out <- .getCpClusterLocLimit(proposed, ex_lookup, 1.5, 0.5)
  # raising to q15 is unchanged
  expect_identical(out$cpJoinTgOrig[1:6], c(2, 2:6))
  # G is lowered to q85 = 6 only down to below 6.8 (frequency 0.08 <= 0.09)
  expect_equal(out$cpJoinTgOrig[[7]], 6.65)
  expect_identical(out$locShareLimit[[7]], "responder_cap")
  # non-donors: half a cell over 1000 cells admits no stimulated cell above
  expect_identical(out$cpJoinTgOrig[8:9], c(5.2, 5.2))
  expect_identical(out$locShareLimit[8:9], rep("nonresponder_cell_cap", 2))
  expect_identical(out$locShareProposed, proposed$cpJoinTgOrig)
  expect_identical(out$cpJoinLse, out$cpJoinTgOrig)
  expect_match(out$locReason[8:9], "_limited_by_nonresponder_cell_cap$")
  expect_match(out$locClusterAction[[7]], "_limited_by_responder_cap$")

  # The donor-median cap uses donors' frequencies at their own gates before
  # sharing (locOwnFreq), not at their current gates. With a loose cell cap,
  # non-donors at q60 = 5 (frequency 0.1) pass when the donors' own median
  # is 0.1 and are limited when it is 0.01.
  loose <- .getCpClusterLocLimit(proposed, ex_lookup, 1.5, 1e6)
  expect_identical(loose$cpJoinTgOrig[8:9], c(5, 5))
  low <- proposed
  low$locOwnFreq <- ifelse(low$locResponder, 0.01, NA_real_)
  lowOut <- .getCpClusterLocLimit(low, ex_lookup, 1.5, 1e6)
  expect_identical(lowOut$cpJoinTgOrig[8:9], c(5.2, 5.2))
  expect_identical(lowOut$locShareLimit[8:9], rep("nonresponder_donor_cap", 2))
})

test_that("stimControl validates the sharing limits", {
  ctrl <- stimControl()
  expect_identical(ctrl$locShareCap, 1.5)
  expect_identical(ctrl$locShareCellCap, 0.5)
  expect_no_error(stimControl(locShareCap = Inf, locShareCellCap = Inf))
  expect_no_error(stimControl(locShareCellCap = 0))
  expect_error(stimControl(locShareCap = 0.9), "at least 1")
  expect_error(stimControl(locShareCellCap = -1), "at least 0")
  expect_error(stimControl(locShareCap = "2"), "at least 1")
})

test_that("gateStim runs with the default sharing limits", {
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  path_project <- tempfile("stimgate_share_")
  withr::defer(unlink(path_project, recursive = TRUE))
  gateStim(
    path_project, gs, exampleData$batchList,
    marker = exampleData$marker
  )
  gates <- getStimGates(path_project)
  expect_true(nrow(gates) > 0L)
  expect_true(all(c("locResponder", "locShareLimit") %in% names(gates)))
  expect_true(all(gates$locShareLimit %in% c(
    "none", "responder_cap", "responder_no_estimate",
    "nonresponder_cell_cap", "nonresponder_donor_cap"
  )))
})


test_that("shared bandwidths widen only for samples below bwNcellMax", {
  settings <- list(bwNcellMax = 1e4)
  expect_equal(.getCpUnsLocBwSharedScale(0.03, 1e4, settings), 0.03)
  expect_equal(.getCpUnsLocBwSharedScale(0.03, 5e4, settings), 0.03)
  expect_equal(
    .getCpUnsLocBwSharedScale(0.03, 1e3, settings),
    0.03 * 10^(1 / 5)
  )
  settings$bwScaleNcell <- FALSE
  expect_equal(.getCpUnsLocBwSharedScale(0.03, 1e3, settings), 0.03)
  expect_true(stimControl()$bwScaleNcell)
  expect_error(stimControl(bwScaleNcell = NA), "bwScaleNcell")
})
