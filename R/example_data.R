#' Get example GatingSet
#'
#' Load the canonical package fixture shipped under
#' \code{inst/extdata/stimgate_test_fixture/}. This keeps the regular package
#' examples and tests on a deterministic, realistic dataset without requiring
#' the simulation machinery to be installed with the package.
#'
#' @param scenario Character ignored for backwards compatibility. Default: "default".
#' @param dirCache Directory ignored for backwards compatibility. Default: NULL.
#' @param clear Logical ignored for backwards compatibility. Default: FALSE.
#' @param nCell Integer ignored for backwards compatibility. Default: 10000.
#' @param nInd Integer ignored for backwards compatibility. Default: 8.
#' @return A list with the saved fixture path, channel labels, marker labels and
#'   sample-to-condition mapping.
#' @export
getExampleData <- function(
  scenario = "default",
  dirCache = NULL,
  clear = FALSE,
  nCell = 1e4,
  nInd = 8
) {
  if (
    !identical(scenario, "default") ||
      !is.null(dirCache) ||
      isTRUE(clear) ||
      !isTRUE(all.equal(as.numeric(nCell), 1e4)) ||
      !isTRUE(all.equal(as.numeric(nInd), 8))
  ) {
    warning(
      "`getExampleData()` now loads the canonical saved package fixture; `scenario`, `dirCache`, `clear`, `nCell`, and `nInd` are deprecated and ignored.",
      call. = FALSE
    )
  }
  .getTestFixture()
}

#' @rdname getExampleData
#' @export
getTestData <- getExampleData

#' Load the package-shipped tiny cytometry fixture
#'
#' Reads the small deterministic FCS files stored in
#' \code{inst/extdata/stimgate_test_fixture/} (or the installed equivalent)
#' and assembles them into a \code{GatingSet} in a temporary directory.
#' Returns the same list structure as \code{\link{getExampleData}} so that
#' tests and examples that do not depend on simulation can substitute this
#' fixture for a much faster run.
#'
#' The fixture contains 2 biological samples, 2 conditions each (unstimulated
#' and stimulated), 2 markers, and roughly 10,000 cells per flow frame with a
#' deliberately rare but non-trivial stimulated response. It is generated
#' deterministically by \code{data-raw/create_test_fixture.R}.
#'
#' @return A named list with elements \code{pathGs}, \code{batchList},
#'   \code{chnl}, and \code{marker}.
#' @keywords internal
.getTestFixture <- function() {
  fixture_dir <- system.file(
    "extdata", "stimgate_test_fixture",
    package = "stimgate"
  )
  if (!nzchar(fixture_dir) || !dir.exists(fixture_dir)) {
    stop(
      "stimgate test fixture not found. ",
      "Run data-raw/create_test_fixture.R to regenerate it."
    )
  }

  meta <- readRDS(file.path(fixture_dir, "metadata.rds"))

  fcs_paths <- file.path(fixture_dir, meta$fcsNames)
  ff_list <- lapply(fcs_paths, flowCore::read.FCS, transformation = FALSE)
  fs <- flowCore::flowSet(ff_list)
  flowCore::sampleNames(fs) <- meta$sampleNames
  gs <- flowWorkspace::GatingSet(fs)

  tmp_dir <- tempfile(pattern = "stimgate_fixture_")
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
