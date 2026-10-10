test_that("coexpression controls are validated and global-only", {
  expect_identical(stimControl()$cytPosMethod, "refine")
  expect_identical(stimControl(cytPosMethod = "coexpression")$coexNBin, 20L)
  for (bad in list(NULL, NA_character_, "unknown", c("refine", "coexpression"))) {
    expect_error(stimControl(cytPosMethod = bad), "cytPosMethod")
  }
  for (nm in c("coexNBin", "coexResidualMin", "coexPurityFrac", "coexZMin")) {
    for (bad in list(NULL, NA_real_, Inf, 0, -1, "2", c(1, 2))) {
      expect_error(do.call(stimControl, stats::setNames(list(bad), nm)), nm)
    }
    expect_error(.resolveMarkerControl(list(a = stats::setNames(list(2), nm)),
      chnl = "a", chnlLab = c(a = "a")), nm)
  }
  expect_error(stimControl(coexNBin = 1.5), "coexNBin")
  expect_error(stimControl(coexPurityFrac = 1.1), "coexPurityFrac")
  expect_error(.resolveMarkerControl(list(a = list(cytPosMethod = "coexpression")),
    chnl = "a", chnlLab = c(a = "a")), "cytPosMethod")
})

test_that("stimulation-specific diagonal lowers gates above the negative floor", {
  negative <- data.frame(a = seq(0.8, 1.2, length.out = 1000),
    b = seq(1.2, 0.8, length.out = 1000))
  response <- data.frame(a = seq(3.3, 5, length.out = 400),
    b = seq(2, 4, length.out = 400))
  dat <- list(stim = rbind(negative, response), uns = negative,
    gate = c(a = 3.2, b = 3.2))
  low <- .coexLowerGate(dat, "a", "b")
  expect_lt(low$cut, dat$gate[["b"]])
  expect_gte(low$cut, low$floorB)
  expect_gte(low$condCut, dat$gate[["a"]])
  cached <- dat
  cached$floor <- vapply(dat$uns, .coexNegFloor, numeric(1))
  cached$above <- lapply(dat[c("stim", "uns")], function(ex) {
    lapply(names(dat$gate), function(m) (ex[[m]] > dat$gate[[m]]) %in% TRUE) |>
      stats::setNames(names(dat$gate))
  })
  expect_identical(.coexLowerGate(cached, "a", "b"), low)
  dat$uns <- dat$stim
  expect_equal(.coexLowerGate(dat, "a", "b")$cut, dat$gate[["b"]])
  dat$uns <- negative
  dat$stim$b <- pmin(dat$stim$b, dat$gate[["b"]])
  expect_equal(.coexLowerGate(dat, "a", "b")$cut, dat$gate[["b"]])
})

test_that("pairwise positivity uses strict cuts without recursive conditioning", {
  gates <- tibble::tibble(chnl = c("a", "b", "c"), batch = "batch1",
    ind = "2", gateName = "loc_minClust", gate = c(3, 3, NA_real_),
    gateCyt = gate)
  attr(gates, "coexpression") <- tibble::tibble(ind = "2", batch = "batch1",
    gateName = "loc_minClust", chnlCond = "a", chnl = "b",
    cut = 1, condCut = 4, lowered = TRUE)
  ex <- tibble::tibble(a = c(4, 4.1, 4.1, 3, 0), b = c(2, 1, 1.1, 3, 3.1), c = 9)
  pos <- .getPosIndByChnl(ex, gates, gateTypeCytPos = "cyt")
  expect_identical(pos$b, c(FALSE, FALSE, TRUE, FALSE, TRUE))
  expect_identical(pos$c, rep(FALSE, 5))
  expect_identical(.getPosIndMult(ex, gates, gateTypeCytPos = "cyt"),
    c(FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_identical(.getPosInd(ex, gates, "b", gateTypeCytPos = "cyt"), pos$b)
  expect_identical(.getPosIndCytCombn(ex, gates, "b", c("a", "c"), "cyt"),
    c(FALSE, FALSE, FALSE, FALSE, TRUE))
})

local({
  withr::local_preserve_seed()
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "coexpression")
  set.seed(42)
  gateStim(pathProject, gs, exampleData$batchList, chnl = exampleData$chnl,
    control = stimControl(cytPosMethod = "coexpression"))

  test_that("example coexpression run saves pairwise gates and tuning", {
    low <- getStimGatesCoexpression(pathProject)
    expect_named(low, c("pop", "batch", "ind", "chnlCond", "markerCond",
      "chnl", "marker", "gate", "cut", "condCut", "floor", "floorCond",
      "z", "purityDp", "lowered"))
    nStim <- sum(lengths(exampleData$batchList) - 1L)
    nChnl <- length(exampleData$chnl)
    expect_equal(nrow(low), nStim * nChnl * (nChnl - 1L))
    expect_false(any(low$chnl == low$chnlCond))
    expect_type(low$lowered, "logical")
    gates <- getStimGates(pathProject)
    expect_identical(gates$gateCyt, gates$gate)
    expect_true(all(low$cut[low$lowered] >= low$floor[low$lowered]))
    settings <- stimgateMetaReadSettingsChnls(pathProject)
    expect_true(all(vapply(settings, function(x) {
      identical(x$cytPosMethod, "coexpression") && x$coexNBin == 20L &&
        x$coexResidualMin == 3.5 && x$coexPurityFrac == 0.75 && x$coexZMin == 2
    }, logical(1))))
  })

  # Recompute all exact combinations directly, without any positivity helper.
  test_that("cached expression reproduces stimulated and raw control combination counts", {
    low <- getStimGatesCoexpression(pathProject)
    gates <- getStimGates(pathProject)
    statistics <- getStimStats(pathProject)
    chnl <- exampleData$chnl
    for (indices in exampleData$batchList) {
      for (ind in indices[-1]) {
        sampleGates <- gates[gates$ind == as.character(ind), ]
        sampleLow <- low[low$ind == as.character(ind), ]
        for (tube in c(ind, indices[[1]])) {
          ex <- getStimExpr(pathProject, ind = tube, chnl = chnl)
          pos <- stats::setNames(lapply(chnl, function(m) {
            g <- sampleGates$gate[match(m, sampleGates$chnl)]
            if (!is.finite(g)) rep(FALSE, nrow(ex)) else (ex[[m]] > g) %in% TRUE
          }), chnl)
          for (i in which(sampleLow$lowered)) {
            a <- sampleLow$chnlCond[[i]]
            b <- sampleLow$chnl[[i]]
            pos[[b]] <- pos[[b]] | ((ex[[a]] > sampleLow$condCut[[i]]) %in% TRUE &
              (ex[[b]] > sampleLow$cut[[i]]) %in% TRUE)
          }
          labels <- vapply(seq_len(nrow(ex)), function(i) {
            paste0(chnl, ifelse(vapply(pos, `[`, logical(1), i), "~+~", "~-~"), collapse = "")
          }, character(1))
          rows <- statistics[statistics$ind == as.character(ind), ]
          counts <- vapply(rows$cytCombn, function(label) sum(labels == label), integer(1))
          expect_equal(unname(counts), if (tube == ind) rows$countStim else rows$countUns)
          expect_equal(sum(counts), nrow(ex))
        }
      }
    }
  })

  test_that("FCS and expression selections use saved coexpression rules", {
    out <- file.path(dirname(exampleData$pathGs), "positive_fcs")
    manifest <- writeStimFCS(pathProject, gs, indBatchList = exampleData$batchList,
      pathDirSave = out)
    expect_s3_class(manifest, "tbl_df")
    expect_true(file.exists(file.path(out, "coexpression.csv")))
    for (indices in exampleData$batchList) {
      for (ind in indices[-1]) {
        selected <- getStimExpr(pathProject, ind = ind, chnl = exampleData$chnl,
          chnlGate = exampleData$chnl)
        expect_equal(manifest$nCellPos[manifest$ind == as.character(ind)], nrow(selected))
      }
    }
  })

  test_that("disabled coexpression preserves ordinary results and getter errors", {
    disabled <- file.path(dirname(exampleData$pathGs), "disabled")
    baseline <- file.path(dirname(exampleData$pathGs), "baseline")
    set.seed(42)
    gateStim(disabled, gs, exampleData$batchList, chnl = exampleData$chnl,
      control = stimControl(cytPosMethod = "coexpression", calcCytPosGates = FALSE))
    set.seed(42)
    gateStim(baseline, gs, exampleData$batchList, chnl = exampleData$chnl,
      control = stimControl(calcCytPosGates = FALSE))
    expect_equal(getStimGates(disabled), getStimGates(baseline))
    expect_equal(getStimStats(disabled), getStimStats(baseline))
    expect_error(getStimGatesCoexpression(disabled), "not gated")
    expect_error(getStimGatesCoexpression(baseline), "not gated")
  })
})
