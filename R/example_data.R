#' Get example GatingSet
#'
#' Load the canonical packaged example dataset shipped under
#' \code{inst/extdata/stimgate_example_data/}. This keeps the regular package
#' examples and tests on a deterministic, realistic dataset without requiring
#' the simulation machinery to be installed with the package.
#'
#' @return A list with the saved example-data path, channel labels, marker labels,
#'   and sample-to-condition mapping.
#' @export
getExampleData <- function() {
  example_dir <- .getExampleDataDir()
  meta <- readRDS(file.path(example_dir, "metadata.rds"))

  fcs_paths <- file.path(example_dir, meta$fcsNames)
  ff_list <- lapply(fcs_paths, flowCore::read.FCS, transformation = FALSE)
  fs <- flowCore::flowSet(ff_list)
  flowCore::sampleNames(fs) <- meta$sampleNames
  gs <- flowWorkspace::GatingSet(fs)

  tmp_dir <- tempfile(pattern = "stimgate_example_data_")
  dir.create(tmp_dir)
  path_gs <- file.path(tmp_dir, "gs")
  flowWorkspace::save_gs(gs, path = path_gs)

  list(
    pathGs = path_gs,
    batchList = meta$batchList,
    chnl = meta$chnl,
    marker = meta$marker
  )
}

#' @rdname getExampleData
#' @export
getTestData <- getExampleData

#' @keywords internal
.getExampleDataDir <- function() {
  example_dir <- system.file(
    "extdata", "stimgate_example_data",
    package = "stimgate"
  )
  if (!nzchar(example_dir) || !dir.exists(example_dir)) {
    stop(
      "stimgate example data not found. ",
      "Run data-raw/create_test_fixture.R to regenerate it."
    )
  }
  example_dir
}
