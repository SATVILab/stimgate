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
  expect_match(note, "tgClust dead plumbing", fixed = TRUE)
  expect_match(note, "legacy single-positive", fixed = TRUE)
})
