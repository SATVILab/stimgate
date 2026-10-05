.acs_path_calls <- function(x) {
  if (!is.call(x) && !is.expression(x) && !is.pairlist(x)) return(list())
  if (is.call(x) && identical(x[[1]], quote(projr::projr_path_get))) return(list(x))
  unlist(lapply(as.list(x), .acs_path_calls), recursive = FALSE)
}

.acs_path_document_calls <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (grepl("\\.qmd$", path)) {
    open <- grep("^```\\{r", lines)
    code <- unlist(lapply(open, function(start) {
      end <- which(lines == "```" & seq_along(lines) > start)[1L]
      lines[seq.int(start + 1L, end - 1L)]
    }))
  } else code <- lines
  .acs_path_calls(parse(text = code))
}

test_that("every ACS projr path remains stable across Quarto working-directory changes", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  files <- c("analysis/9-real-compare-acs-cytof.qmd",
    "analysis/10-real-compare-acs-cytof-validation.qmd",
    "scripts/r/acs_cytof-manual.R", "scripts/r/acs_cytof-gate.R")
  calls <- unlist(lapply(file.path(root, files), .acs_path_document_calls), recursive = FALSE)
  expect_length(calls, 7L)
  for (call in calls) expect_identical(as.list(call)[["format"]], "absolute")
  testthat::skip_if_not_installed("projr")
  withr::local_envvar(PROJR_PROFILE = NA)
  fixture <- withr::local_tempdir()
  project <- file.path(fixture, "stimgate")
  analysis <- file.path(project, "analysis")
  dir.create(analysis, recursive = TRUE)
  dir.create(file.path(project, ".git"))
  writeLines(c("directories:", "  raw-data-large:",
    "    path: ../stimgate-store/_raw/data/large", "    ignore-git: no",
    "  raw-data-small:", "    path: _raw/data/small", "    ignore-git: no",
    "  cache:", "    path: ../stimgate-store/_tmp", "    ignore-git: no"),
    file.path(project, "_projr.yml"))
  fcs <- file.path(fixture, "stimgate-store", "_raw", "data", "large",
    "comparison_data", "acscytof", "fcs", "tcrgd")
  dir.create(fcs, recursive = TRUE)
  file.create(file.path(fcs, paste0(1:10, ".fcs")))
  # Evaluate the real QMD/helper calls where setup initially runs. Do not create
  # cache outputs: those absent targets must also resolve absolutely.
  env <- new.env(parent = baseenv())
  env$fn <- "manual.csv"
  paths <- withr::with_dir(analysis, lapply(calls, function(call) {
    call <- as.call(c(as.list(call), list(create = FALSE)))
    eval(call, env)
  }))
  withr::with_dir(project, {
    gate <- new.env(parent = baseenv())
    source(file.path(root, "scripts/r/acs_cytof-gate.R"), local = gate)
    expect_length(gate$.acsCytofFcsFiles(file.path(paths[[1]], "tcrgd")), 10L)
    for (path in paths) expect_true(grepl("^(/|[A-Za-z]:[/\\\\])", path))
    expect_equal(normalizePath(file.path(paths[[1]], "tcrgd"), winslash = "/"),
      normalizePath(fcs, winslash = "/"))
    expect_false(dir.exists(paths[[2]]))
    expect_equal(paths[[4]], paths[[5]])
  })
})
