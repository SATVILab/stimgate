if (!Sys.getenv("STIMGATE_CLEAR_TEST_DATA") == "FALSE") {
  pathDirCacheTestData <- testthat::test_path("cache", "test_data")
  if (dir.exists(pathDirCacheTestData)) {
    unlink(pathDirCacheTestData, recursive = TRUE, force = TRUE)
  }
}
