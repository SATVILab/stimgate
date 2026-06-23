if (!Sys.getenv("STIMGATE_CLEAR_TEST_DATA") == "FALSE") {
  path_dir_cache_test_data <- testthat::test_path("cache", "test_data")
  if (dir.exists(path_dir_cache_test_data)) {
    unlink(path_dir_cache_test_data, recursive = TRUE, force = TRUE)
  }
}
