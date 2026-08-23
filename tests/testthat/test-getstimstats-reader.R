test_that("getStimStats returns saved RDS table unchanged", {
  tmp_dir <- file.path(tempdir(), "stimgate_reader_rds_test")
  dir.create(tmp_dir, showWarnings = FALSE)
  withr::defer(unlink(tmp_dir, recursive = TRUE))

  # Fixture uses non-character ind (integer) and factor batch to verify the
  # reader contract: getStimStats() must return exactly what is stored.
  fixture <- data.frame(
    ind = 1:2,
    batch = factor(c("batch1", "batch1")),
    chnl = c("APC-A", "APC-A"),
    marker = c("IFNg", "IFNg"),
    count_pos = c(42L, 7L),
    freq_pos = c(0.042, 0.007),
    stringsAsFactors = FALSE
  )
  saveRDS(fixture, file.path(tmp_dir, "gateStats.rds"))

  result <- getStimStats(tmp_dir)

  expect_equal(result, fixture)
})

test_that("getStimStats reads CSV fallback when no RDS present", {
  tmp_dir <- file.path(tempdir(), "stimgate_reader_csv_test")
  dir.create(tmp_dir, showWarnings = FALSE)
  withr::defer(unlink(tmp_dir, recursive = TRUE))

  fixture <- data.frame(
    ind = c("sampleX", "sampleY"),
    count_pos = c(5L, 99L),
    stringsAsFactors = FALSE
  )
  utils::write.csv(fixture, file.path(tmp_dir, "gateStats.csv"), row.names = FALSE)

  result <- getStimStats(tmp_dir)

  expect_true(is.data.frame(result))
  expect_equal(result[["ind"]], fixture[["ind"]])
  expect_equal(result[["count_pos"]], fixture[["count_pos"]])
})

test_that("getStimStats prefers RDS over CSV when both exist", {
  tmp_dir <- file.path(tempdir(), "stimgate_reader_precedence_test")
  dir.create(tmp_dir, showWarnings = FALSE)
  withr::defer(unlink(tmp_dir, recursive = TRUE))

  rds_fixture <- data.frame(
    ind = "from_rds",
    count_pos = 1L,
    stringsAsFactors = FALSE
  )
  csv_fixture <- data.frame(
    ind = "from_csv",
    count_pos = 999L,
    stringsAsFactors = FALSE
  )
  saveRDS(rds_fixture, file.path(tmp_dir, "gateStats.rds"))
  utils::write.csv(csv_fixture, file.path(tmp_dir, "gateStats.csv"), row.names = FALSE)

  result <- getStimStats(tmp_dir)

  expect_equal(result[["ind"]], "from_rds")
})

test_that("getStimStats errors with 'No stats file found' when project has no stats", {
  tmp_dir <- file.path(tempdir(), "stimgate_reader_nostats_test")
  dir.create(tmp_dir, showWarnings = FALSE)
  withr::defer(unlink(tmp_dir, recursive = TRUE))

  expect_error(getStimStats(tmp_dir), "No stats file found")
})
