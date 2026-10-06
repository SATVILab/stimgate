test_that("getStimExpr reads saved channel data and filters correctly", {
  tmp <- tempfile("stimgate_ex_")
  dir.create(
    file.path(tmp, "sampleData", "pop_POP1", "ind_1"),
    recursive = TRUE
  )
  dir.create(
    file.path(tmp, "sampleData", "pop_POP1", "ind_2"),
    recursive = TRUE
  )

  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )
  saveRDS(
    c(4, 5, 6),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC2.rds")
  )
  saveRDS(
    c(7, 8),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_2", "chnl_BC1.rds")
  )
  saveRDS(
    c(9, 10),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_2", "chnl_BC2.rds")
  )
  res <- getStimExpr(tmp)
  expect_equal(nrow(res), 5)
  expect_true(all(c("pop", "ind", "BC1", "BC2") %in% names(res)))
  expect_equal(unique(res$pop), "POP1")
  expect_equal(sum(res$ind == "1"), 3)
  expect_equal(sum(res$ind == "2"), 2)

  resBc1 <- getStimExpr(tmp, chnl = "BC1")
  expect_true("BC1" %in% names(resBc1))
  expect_false("BC2" %in% names(resBc1))

  resInd1 <- getStimExpr(tmp, ind = "1")
  expect_equal(nrow(resInd1), 3)

  expect_error(getStimExpr(""))
})

test_that("getStimExpr applies bias only to unstim sample", {
  tmp <- tempfile("stimgate_ex_bias_")
  dir.create(
    file.path(tmp, "sampleData", "pop_POP1", "ind_1"),
    recursive = TRUE
  )
  dir.create(
    file.path(tmp, "sampleData", "pop_POP1", "ind_2"),
    recursive = TRUE
  )

  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )
  saveRDS(
    c(4, 5, 6),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC2.rds")
  )
  saveRDS(
    c(7, 8),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_2", "chnl_BC1.rds")
  )
  saveRDS(
    c(9, 10),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_2", "chnl_BC2.rds")
  )

  # Create metaData with matching names so chnl lookup works
  chnlList <- list(BC1 = list(biasUns = 10), BC2 = list(biasUns = -2))
  chnlLab <- c(BC1 = "BC1", BC2 = "BC2")
  batchList <- list(batch1 = c("1", "2"))

  dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
  saveRDS(chnlList, file.path(tmp, "metaData", "chnlSettings.rds"))
  saveRDS(chnlLab, file.path(tmp, "metaData", "chnlLab.rds"))
  saveRDS(batchList, file.path(tmp, "metaData", "batchList.rds"))

  # unstim is ind 1 -> expect bias added
  resUns <- getStimExpr(tmp, ind = "1", bias = TRUE)
  expect_equal(resUns$BC1, c(1 + 10, 2 + 10, 3 + 10))
  expect_equal(resUns$BC2, c(4 - 2, 5 - 2, 6 - 2))

  # stim is ind 2 -> bias should not be applied
  resStim <- getStimExpr(tmp, ind = "2", bias = TRUE)
  expect_equal(resStim$BC1, c(7, 8))
  expect_equal(resStim$BC2, c(9, 10))
})

test_that("getStimExpr uses completed marker settings and tolerates absent bias", {
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- withr::local_tempdir()
  invisible(gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  ))

  settings <- stimgateMetaReadSettingsChnls(pathProject)
  chnlLab <- stimgateMetaReadChnlLab(pathProject)
  expect_setequal(names(settings), unname(chnlLab[exampleData$chnl]))
  expect_false(any(exampleData$chnl %in% names(settings)))
  batch <- exampleData$batchList[[1]]
  exUns <- getStimExpr(pathProject, ind = as.character(batch[[1]]))
  exBias <- getStimExpr(pathProject, ind = as.character(batch[[1]]), bias = TRUE)
  for (ch in exampleData$chnl) {
    savedBias <- settings[[chnlLab[[ch]]]]$biasUns
    expect_length(savedBias, 1L)
    expect_equal(exBias[[ch]], exUns[[ch]] + savedBias)
  }
  indStim <- as.character(batch[[2]])
  expect_equal(
    getStimExpr(pathProject, ind = indStim, bias = TRUE),
    getStimExpr(pathProject, ind = indStim)
  )

  # Older or incomplete settings can lack a bias or a channel entirely.
  settings[[chnlLab[[exampleData$chnl[[1]]]]]]$biasUns <- NULL
  settings[[chnlLab[[exampleData$chnl[[2]]]]]] <- NULL
  saveRDS(settings, file.path(pathProject, "metaData", "chnlSettings.rds"))
  expect_equal(
    getStimExpr(pathProject, ind = as.character(batch[[1]]), bias = TRUE),
    exUns
  )
})

test_that("getStimExpr excludes minimum observed values when excMin = TRUE", {
  tmp <- tempfile("stimgate_ex_excmin_")
  dir.create(
    file.path(tmp, "sampleData", "pop_POP1", "ind_1"),
    recursive = TRUE
  )

  # BC1 min is 1, BC2 min is 4; after exclusion only the row (3,6) should remain
  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )
  saveRDS(
    c(4, 4, 6),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC2.rds")
  )

  resNoexc <- getStimExpr(tmp, ind = "1", excMin = FALSE)
  expect_equal(nrow(resNoexc), 3)

  resExc <- getStimExpr(tmp, ind = "1", excMin = TRUE)
  expect_equal(nrow(resExc), 1)
  expect_equal(resExc$BC1, 3)
  expect_equal(resExc$BC2, 6)
})

test_that("getStimExpr uses marker parameter to rename channels", {
  tmp <- tempfile("stimgate_ex_marker_")
  dir.create(
    file.path(tmp, "sampleData", "pop_POP1", "ind_1"),
    recursive = TRUE
  )

  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )
  saveRDS(
    c(4, 5, 6),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC2.rds")
  )

  # Create chnlLab mapping (channel -> marker name)
  chnlLab <- c(BC1 = "IFNg", BC2 = "IL2")
  dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
  saveRDS(chnlLab, file.path(tmp, "metaData", "chnlLab.rds"))

  res <- getStimExpr(tmp, marker = "IFNg")
  expect_true("IFNg" %in% names(res))
  expect_false("BC1" %in% names(res))
  expect_equal(res$IFNg, c(1, 2, 3))
})

test_that("getStimExpr errors when both marker and chnl specified", {
  tmp <- tempfile("stimgate_ex_marker_chnl_conflict_")
  dir.create(
    file.path(tmp, "sampleData", "pop_POP1", "ind_1"),
    recursive = TRUE
  )

  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )

  expect_error(
    getStimExpr(tmp, marker = "IFNg", chnl = "BC1"),
    "Must not specify both marker and chnl"
  )
})

test_that("getStimExpr applies transFn to specified channels", {
  tmp <- tempfile("stimgate_ex_trans_")
  dir.create(
    file.path(tmp, "sampleData", "pop_POP1", "ind_1"),
    recursive = TRUE
  )

  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )
  saveRDS(
    c(4, 5, 6),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC2.rds")
  )

  # Test transforming all channels
  transFnDouble <- function(x) x * 2
  resAll <- getStimExpr(tmp, transFn = transFnDouble)
  expect_equal(resAll$BC1, c(2, 4, 6))
  expect_equal(resAll$BC2, c(8, 10, 12))

  # Test transforming only BC1
  resBc1 <- getStimExpr(
    tmp,
    transFn = transFnDouble,
    transChnl = "BC1"
  )
  expect_equal(resBc1$BC1, c(2, 4, 6))
  expect_equal(resBc1$BC2, c(4, 5, 6))
})

test_that("getStimExpr applies transFn to markers when using marker parameter", {
  tmp <- tempfile("stimgate_ex_trans_marker_")
  dir.create(
    file.path(tmp, "sampleData", "pop_POP1", "ind_1"),
    recursive = TRUE
  )

  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )
  saveRDS(
    c(4, 5, 6),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC2.rds")
  )

  # Create chnlLab mapping (channel -> marker name)
  chnlLab <- c(BC1 = "IFNg", BC2 = "IL2")
  dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
  saveRDS(chnlLab, file.path(tmp, "metaData", "chnlLab.rds"))

  transFnDouble <- function(x) x * 2
  res <- getStimExpr(
    tmp,
    marker = c("IFNg", "IL2"),
    transFn = transFnDouble,
    transMarker = "IFNg"
  )
  expect_equal(res$IFNg, c(2, 4, 6))
  expect_equal(res$IL2, c(4, 5, 6))
})

test_that("getStimExpr errors when both chnlGate and markerGate specified", {
  tmp <- tempfile("stimgate_ex_gate_conflict_")
  dir.create(
    file.path(tmp, "sampleData", "pop_POP1", "ind_1"),
    recursive = TRUE
  )

  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sampleData", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )

  # Create chnlLab mapping for markerGate
  chnlLab <- c(BC1 = "IFNg")
  dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
  saveRDS(chnlLab, file.path(tmp, "metaData", "chnlLab.rds"))

  expect_error(
    getStimExpr(tmp, chnlGate = "BC1", markerGate = "IFNg"),
    "Must not specify both chnlGate and markerGate"
  )
})

test_that("getStimExpr and plotStim filter using saved stimulation gates", {
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- withr::local_tempdir()
  invisible(gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker
  ))

  gateTbl <- getStimGates(pathProject)
  chnlGated <- unique(gateTbl$chnl)
  chnlLab <- stimgateMetaReadChnlLab(pathProject)
  expect_gt(nrow(gateTbl), 0L)

  # Include unstimulated indices: they have expression but no saved gates.
  for (chnlCurr in chnlGated) {
    markerCurr <- unname(chnlLab[chnlCurr])
    exAll <- getStimExpr(pathProject, pop = "root", chnl = chnlCurr)
    for (gateType in c("base", "cyt")) {
      expected <- rep(FALSE, nrow(exAll))
      for (indCurr in unique(exAll$ind)) {
        gateInd <- gateTbl[
          gateTbl$chnl == chnlCurr & gateTbl$ind == indCurr,
        ]
        if (nrow(gateInd) == 0L) {
          next
        }
        rows <- exAll$ind == indCurr
        # One requested channel has no other cytokine context, so the cyt+
        # rule reduces to base positivity as well.
        expected[rows] <- exAll[[chnlCurr]][rows] > gateInd$gate[[1]]
      }
      res <- getStimExpr(
        pathProject,
        pop = "root",
        chnl = chnlCurr,
        chnlGate = chnlCurr,
        gateTypeCytPos = gateType
      )
      expect_equal(nrow(res), sum(expected))
      expect_equal(res$ind, exAll$ind[expected])
      expect_equal(res[[chnlCurr]], exAll[[chnlCurr]][expected])

      resMarker <- getStimExpr(
        pathProject,
        pop = "root",
        marker = markerCurr,
        markerGate = markerCurr,
        gateTypeCytPos = gateType
      )
      expect_equal(resMarker$ind, res$ind)
      expect_equal(resMarker[[markerCurr]], res[[chnlCurr]])
    }

    # The filter must keep a meaningful number of cells for some stimulated
    # sample (the counts per sample depend on the bandwidth defaults).
    expect_gt(max(table(res$ind)), 10L)
    for (useMarker in c(FALSE, TRUE)) {
      args <- if (useMarker) {
        list(marker = markerCurr, markerGate = markerCurr)
      } else {
        list(chnl = chnlCurr, chnlGate = chnlCurr)
      }
      p <- do.call(plotStim, c(list(
        pathProject = pathProject,
        .data = gs,
        ind = unique(exAll$ind),
        excMin = FALSE,
        grid = FALSE
      ), args))
      expect_length(p, 1L)
      expect_s3_class(p[[1]], "ggplot")
    }
  }

  exAll <- getStimExpr(pathProject, pop = "root", chnl = chnlGated)
  expected <- rep(FALSE, nrow(exAll))
  for (indCurr in unique(exAll$ind)) {
    gates <- gateTbl[gateTbl$ind == indCurr, ]
    if (nrow(gates) == 0L) {
      next
    }
    rows <- exAll$ind == indCurr
    base <- lapply(chnlGated, function(ch) {
      exAll[[ch]][rows] > gates$gate[gates$chnl == ch]
    })
    cyt <- lapply(chnlGated, function(ch) {
      exAll[[ch]][rows] > gates$gateCyt[gates$chnl == ch]
    })
    # Multifunctional cyt+ positivity requires base positivity on one
    # channel and cyt+ positivity on a different channel.
    pos <- rep(FALSE, sum(rows))
    for (i in seq_along(chnlGated)) {
      for (j in setdiff(seq_along(chnlGated), i)) {
        pos <- pos | (base[[i]] & cyt[[j]])
      }
    }
    expected[rows] <- pos
  }
  resMult <- getStimExpr(
    pathProject,
    pop = "root",
    chnl = chnlGated,
    chnlGate = chnlGated,
    mult = TRUE
  )
  expect_equal(nrow(resMult), sum(expected))
  expect_equal(resMult$ind, exAll$ind[expected])
  for (ch in chnlGated) {
    expect_equal(resMult[[ch]], exAll[[ch]][expected])
  }
})

test_that("cytokine-positive filtering helpers respect gate context and exclusions", {
  ex <- tibble::tibble(
    IFNg = c(0, 5, 3, 3, 5, 5, 3),
    IL2 = c(0, 0, 5, 0, 5, 0, 0),
    TNFa = c(0, 0, 0, 5, 0, 5, 0)
  )
  gateTbl <- tibble::tibble(
    chnl = c("IFNg", "IL2", "TNFa"),
    ind = "1",
    gate = c(4, 4, 4),
    gateCyt = c(2, 2, 2)
  )

  exBase <- .dataGetExCytPosInc(
    ex,
    gateTbl,
    mult = FALSE,
    chnl = c("IFNg", "IL2", "TNFa"),
    gateTypeCytPos = "base"
  )
  expect_equal(
    exBase[, c("IFNg", "IL2", "TNFa")],
    tibble::tibble(
      IFNg = c(5, 3, 3, 5, 5),
      IL2 = c(0, 5, 0, 5, 0),
      TNFa = c(0, 0, 5, 0, 5)
    )
  )

  exCyt <- .dataGetExCytPosInc(
    ex,
    gateTbl,
    mult = FALSE,
    chnl = c("IFNg", "IL2", "TNFa"),
    gateTypeCytPos = "cyt"
  )
  expect_equal(
    exCyt[, c("IFNg", "IL2", "TNFa")],
    tibble::tibble(
      IFNg = c(5, 3, 3, 5, 5),
      IL2 = c(0, 5, 0, 5, 0),
      TNFa = c(0, 0, 5, 0, 5)
    )
  )

  exCytMult <- .dataGetExCytPosInc(
    ex,
    gateTbl,
    mult = TRUE,
    chnl = c("IFNg", "IL2", "TNFa"),
    gateTypeCytPos = "cyt"
  )
  expect_equal(
    exCytMult[, c("IFNg", "IL2", "TNFa")],
    tibble::tibble(
      IFNg = c(3, 3, 5, 5),
      IL2 = c(5, 0, 5, 0),
      TNFa = c(0, 5, 0, 5)
    )
  )

  exSubset <- .dataGetExCytPosInc(
    ex,
    gateTbl,
    mult = FALSE,
    chnl = "IFNg",
    gateTypeCytPos = "cyt"
  )
  expect_equal(
    nrow(exSubset),
    5L
  )
  expect_equal(
    exSubset$IFNg,
    c(5, 3, 3, 5, 5)
  )
  expect_false(any(exSubset$IFNg == 3 & exSubset$IL2 == 0 & exSubset$TNFa == 0))

  exGateCytRequired <- tibble::tibble(
    IFNg = c(1, 3, 5),
    IL2 = c(5, 5, 0),
    TNFa = c(0, 0, 0)
  )
  exGateCytRequiredPos <- .dataGetExCytPosInc(
    exGateCytRequired,
    gateTbl,
    mult = FALSE,
    chnl = "IFNg",
    gateTypeCytPos = "cyt"
  )
  expect_true(any(exGateCytRequiredPos$IFNg == 3 & exGateCytRequiredPos$IL2 == 5))
  expect_false(any(exGateCytRequiredPos$IFNg == 1 & exGateCytRequiredPos$IL2 == 5))

  exExc <- .dataGetExCytPosExc(
    exCyt,
    combnExc = list(c("IFNg", "IL2")),
    gateTblInd = gateTbl,
    chnlGate = c("IFNg", "IL2", "TNFa"),
    gateTypeCytPos = "cyt"
  )
  expect_equal(
    exExc[, c("IFNg", "IL2", "TNFa")],
    tibble::tibble(
      IFNg = c(5, 3, 5),
      IL2 = c(0, 0, 0),
      TNFa = c(0, 5, 5)
    )
  )

  tmp <- tempfile("stimgate_ex_cyt_filter_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  dir.create(file.path(tmp, "sampleData", "pop_root", "ind_1"), recursive = TRUE)
  saveRDS(
    c(0, 5, 3, 3, 5, 5, 3),
    file.path(tmp, "sampleData", "pop_root", "ind_1", "chnl_IFNg.rds")
  )
  saveRDS(
    c(0, 0, 5, 0, 5, 0, 0),
    file.path(tmp, "sampleData", "pop_root", "ind_1", "chnl_IL2.rds")
  )
  saveRDS(
    c(0, 0, 0, 5, 0, 5, 0),
    file.path(tmp, "sampleData", "pop_root", "ind_1", "chnl_TNFa.rds")
  )
  for (nm in c("IFNg", "IL2", "TNFa")) {
    dir.create(
      file.path(tmp, "gates", "poproot", paste0("chnl", nm), "all"),
      recursive = TRUE
    )
    gateTblCurr <- tibble::tibble(
      chnl = nm,
      ind = "1",
      gate = 4,
      gateCyt = 2
    )
    saveRDS(
      gateTblCurr,
      file.path(tmp, "gates", "poproot", paste0("chnl", nm), "all", "gateTbl.rds")
    )
  }

  res <- getStimExpr(
    tmp,
    pop = "root",
    chnl = c("IFNg", "IL2", "TNFa"),
    chnlGate = c("IFNg", "IL2", "TNFa")
  )
  expect_equal(
    res[, c("IFNg", "IL2", "TNFa")],
    tibble::tibble(
      IFNg = c(5, 3, 3, 5, 5),
      IL2 = c(0, 5, 0, 5, 0),
      TNFa = c(0, 0, 5, 0, 5)
    ),
    check.attributes = FALSE
  )
})

test_that("getStimExpr returns zero rows when no cells are stimulation-positive", {
  tmp <- tempfile("stimgate_ex_zero_pos_")
  withr::defer(unlink(tmp, recursive = TRUE))
  dir.create(file.path(tmp, "sampleData", "pop_root", "ind_1"), recursive = TRUE)
  saveRDS(
    c(0, 5, 3, 3, 5, 5, 3),
    file.path(tmp, "sampleData", "pop_root", "ind_1", "chnl_IFNg.rds")
  )
  dir.create(
    file.path(tmp, "gates", "poproot", "chnlIFNg", "all"),
    recursive = TRUE
  )
  saveRDS(
    tibble::tibble(chnl = "IFNg", ind = "1", gate = 100, gateCyt = 100),
    file.path(tmp, "gates", "poproot", "chnlIFNg", "all", "gateTbl.rds")
  )

  res <- suppressMessages(getStimExpr(
    tmp,
    pop = "root",
    chnl = "IFNg",
    chnlGate = "IFNg",
    excMin = TRUE
  ))
  expect_identical(nrow(res), 0L)
  expect_identical(colnames(res), c("pop", "ind", "IFNg"))
  expect_equal(attr(res, "probGMin")[["root"]][["1"]][["IFNg"]], 6 / 7)
  expect_equal(
    attr(res, "nCellPos"),
    tibble::tibble(pop = "root", ind = "1", nCellPos = 0L)
  )
})

test_that(
  "getStimExpr attaches nCellPos attribute for every requested pop and ind",
  {
    tmp <- tempfile("stimgate_ex_ncell_")
    withr::defer(unlink(tmp, recursive = TRUE))

    dir.create(
      file.path(tmp, "sampleData", "pop_root", "ind_1"),
      recursive = TRUE
    )
    dir.create(
      file.path(tmp, "sampleData", "pop_root", "ind_2"),
      recursive = TRUE
    )
    dir.create(
      file.path(tmp, "sampleData", "pop_sub", "ind_1"),
      recursive = TRUE
    )
    dir.create(
      file.path(tmp, "sampleData", "pop_sub", "ind_2"),
      recursive = TRUE
    )

    saveRDS(
      c(1, 2, 10, 20),
      file.path(tmp, "sampleData", "pop_root", "ind_1", "chnl_IFNg.rds")
    )
    saveRDS(
      c(1, 2, 3),
      file.path(tmp, "sampleData", "pop_root", "ind_2", "chnl_IFNg.rds")
    )
    saveRDS(
      c(10, 20, 30),
      file.path(tmp, "sampleData", "pop_sub", "ind_1", "chnl_IFNg.rds")
    )
    saveRDS(
      c(1, 2),
      file.path(tmp, "sampleData", "pop_sub", "ind_2", "chnl_IFNg.rds")
    )

    dir.create(
      file.path(tmp, "gates", "poproot", "chnlIFNg", "all"),
      recursive = TRUE
    )
    saveRDS(
      tibble::tibble(
        chnl = "IFNg",
        ind = c("1", "2"),
        gate = c(5, 5),
        gateCyt = c(5, 5)
      ),
      file.path(tmp, "gates", "poproot", "chnlIFNg", "all", "gateTbl.rds")
    )

    dir.create(
      file.path(tmp, "gates", "popsub", "chnlIFNg", "all"),
      recursive = TRUE
    )
    saveRDS(
      tibble::tibble(
        chnl = "IFNg",
        ind = c("1", "2"),
        gate = c(5, 5),
        gateCyt = c(5, 5)
      ),
      file.path(tmp, "gates", "popsub", "chnlIFNg", "all", "gateTbl.rds")
    )

    res <- suppressMessages(getStimExpr(
      tmp,
      pop = c("root", "sub"),
      ind = c("1", "2"),
      chnl = "IFNg",
      chnlGate = "IFNg"
    ))

    # Samples with zero positive cells return zero rows
    expect_equal(nrow(res[res$pop == "root" & res$ind == "2", ]), 0L)
    expect_equal(nrow(res[res$pop == "sub" & res$ind == "2", ]), 0L)
    expect_equal(nrow(res[res$pop == "root" & res$ind == "1", ]), 2L)
    expect_equal(nrow(res[res$pop == "sub" & res$ind == "1", ]), 3L)

    # Attribute nCellPos includes EVERY requested pop/ind combination
    n_cell_pos <- attr(res, "nCellPos")
    expect_s3_class(n_cell_pos, "tbl_df")
    expect_named(n_cell_pos, c("pop", "ind", "nCellPos"))
    expect_type(n_cell_pos$pop, "character")
    expect_type(n_cell_pos$ind, "character")
    expect_type(n_cell_pos$nCellPos, "integer")

    expected_n_cell_pos <- tibble::tibble(
      pop = c("root", "root", "sub", "sub"),
      ind = c("1", "2", "1", "2"),
      nCellPos = c(2L, 0L, 3L, 0L)
    )
    expect_equal(n_cell_pos, expected_n_cell_pos)
  }
)
