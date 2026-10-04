#' Load example cytometry data
#'
#' Load the packaged dataset and save a GatingSet in a temporary directory.
#' Use the returned paths and labels to try [gateStim()].
#'
#' @return A list with `pathGs` (saved GatingSet path), `batchList` (sample
#'   indices by batch, control first), `chnl` (channels) and `marker` (labels).
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' exampleData$batchList
#' @export
getExampleData <- function() {
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
