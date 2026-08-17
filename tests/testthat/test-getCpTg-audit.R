test_that(".get_cp_tg_call_audit enumerates live, optional and dead paths", {
  audit_tbl <- stimgate:::.get_cp_tg_call_audit()

  expect_s3_class(audit_tbl, "tbl_df")
  expect_equal(nrow(audit_tbl), 5)
  expect_true(all(c(
    "callSite",
    "activation",
    "downstreamConsumer",
    "observableOutput",
    "classification",
    "intermediateArtifacts"
  ) %in% names(audit_tbl)))

  expect_true(any(audit_tbl$classification == "dead for current output"))
  expect_true(any(
    audit_tbl$classification == "backwards-compatible but optional"
  ))
  expect_true(any(grepl("tolClust", audit_tbl$callSite, fixed = TRUE)))
  expect_true(any(grepl("tgType = 'tg'", audit_tbl$callSite, fixed = TRUE)))
})

test_that(".get_cp_tg_migration_note_157 summarises migration work", {
  note <- stimgate:::.get_cp_tg_migration_note_157()

  expect_true(is.character(note))
  expect_length(note, 1)
  expect_match(note, "Issue #157 migration note", fixed = TRUE)
  expect_match(note, "default init tgClust tailgate path has been removed", fixed = TRUE)
  expect_match(note, "legacy single-positive", fixed = TRUE)
})

test_that("default init flow no longer generates tgClust gate rows", {
  chnlSettings <- list(
    tolCtrl = 0.01,
    tolClust = 0.02,
    chnlCut = "marker",
    gateCombn = "prejoin",
    excMin = FALSE,
    minCell = 5,
    cpMin = 0,
    bw = 1,
    tol = 1e-2
  )

  exList <- list(
    uns = data.frame(marker = c(1, 2, 3)),
    stim = data.frame(marker = c(2, 3, 4))
  )
  attr(exList[[1]], "batch") <- "batch1"
  attr(exList[[2]], "batch") <- "batch1"

  local_mocked_bindings(
    .getCpUnsLoc = function(...) list(list(cp = c(uns = 1))),
    .getCpTg = function(exList, chnlSettings, tgType, stage, pathProject) {
      if (identical(tgType, "tolClust")) {
        stop("should not evaluate dead tgClust tailgate path")
      }
      list(cp = c(uns = 1))
    },
    .gateBatchTbl = function(gateList, batch) gateList
  )

  gateList <- stimgate:::.gateBatchAll(
    indBatch = 1,
    batch = "batch1",
    exList = exList,
    .data = NULL,
    chnlSettings = chnlSettings,
    stage = "init",
    pathProject = tempdir()
  )

  expect_false(any(grepl("tgClust", names(gateList), fixed = TRUE)))
  expect_equal(stimgate:::.gateBatchTblUse("tgCtrl_0.01"), "ctrl")
  expect_equal(stimgate:::.gateBatchTblUse("tgClust"), "gate")
})
