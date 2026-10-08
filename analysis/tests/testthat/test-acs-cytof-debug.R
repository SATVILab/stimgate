root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
qmd13_path <- file.path(root_dir, "analysis", "13-real-debug-acs-cytof.qmd")
qmd9_path <- file.path(root_dir, "analysis", "9-real-compare-acs-cytof.qmd")

.load_acs_debug_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "analysis-plot-style.R", "acs_cytof-helper.R",
    "acs_cytof-gate.R", "acs_cytof-methods.R", "sim-debug-loc.R",
    "acs_cytof-debug.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  env
}

# The example data as a two-donor ACS-like population, with ACS helpers
# mocked only where the example's layout (two tubes per donor, two markers)
# differs from ACS.
.mock_acs_example <- function(env, ex) {
  sampleMap <- tibble::tibble(
    ind = as.character(unlist(ex$batchList)),
    SampleID = rep(c("d1", "d2"), lengths(ex$batchList)),
    stim = unlist(lapply(ex$batchList, function(b) c("uns", rep("p1", length(b) - 1L))))
  )
  env$.acsCytofEnsureCurrentCheckout <- function(...) invisible(TRUE)
  env$.acsCytofReadPreprocessing <- function(...) list(sampleMap = sampleMap)
  env$.acsCytofBatchList <- function(...) unname(ex$batchList)
  env$.acsCytofChannelMap <- function() stats::setNames(ex$marker, ex$chnl)
  env$.acsCytofGateStim <- function(gs, pathProject, batchList, ...) {
    stimgate::gateStim(
      pathProject = pathProject, .data = gs, batchList = batchList,
      chnl = ex$chnl, control = stimgate::stimControl(clusterGates = TRUE)
    )
  }
  sampleMap
}

.quiet <- function(code) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(code)))
  out
}

test_that(".simDebugLoc keys every tube by channel and maps real donors", {
  env <- .load_acs_debug_env()
  ex <- getExampleData()
  withr::defer(unlink(dirname(ex$pathGs), recursive = TRUE))
  sampleMap <- .mock_acs_example(env, ex)
  gs <- flowWorkspace::load_gs(ex$pathGs)
  pathRef <- withr::local_tempdir()
  pathDbg <- withr::local_tempdir()

  withr::local_seed(11)
  .quiet(env$.acsCytofGateStim(gs, pathRef, unname(ex$batchList)))
  withr::local_seed(11)
  info <- env$.acsCytofDebugTubeInfo(sampleMap, "tcrgd")
  out <- .quiet(env$.simDebugLoc(
    env$.acsCytofGateStim(gs, pathDbg, unname(ex$batchList)),
    sample = NULL,
    tubeInfo = info,
    onRecord = function(rec) rec
  ))
  # Tracing changes nothing: the same gates as an untraced run.
  expect_equal(getStimGates(pathDbg), getStimGates(pathRef))

  expect_s3_class(out, "simDebugLocList")
  indStim <- as.character(unlist(lapply(ex$batchList, `[`, -1L)))
  expect_setequal(
    names(out),
    paste0("dataset1_ind", rep(indStim, each = 2L), "_", ex$chnl)
  )
  for (rec in out) {
    expect_identical(rec$dataset, 1L)
    expect_identical(rec$chnl, attr(rec$inputs$exTblStimNoMin, "chnlCut"))
    expect_identical(rec$sample, sampleMap$SampleID[sampleMap$ind == rec$ind])
    expect_identical(rec$tube$stim, "p1")
    expect_null(rec$truth)
  }
  expect_identical(info("1")$stim, sampleMap$stim[sampleMap$ind == "1"])
  expect_error(info("99"), "not in the ACS sample map")

  kept <- .quiet(env$.simDebugLoc(
    env$.acsCytofGateStim(gs, withr::local_tempdir(), unname(ex$batchList)),
    sample = NULL,
    onRecord = function(rec) rec$chnl
  ))
  expect_setequal(unlist(kept), ex$chnl)
})

test_that("pages draw without simulated truth and record the final gates", {
  env <- .load_acs_debug_env()
  ex <- getExampleData()
  withr::defer(unlink(dirname(ex$pathGs), recursive = TRUE))
  .mock_acs_example(env, ex)
  base <- withr::local_tempdir()
  paths9 <- list(
    gs = ex$pathGs, stimgate = file.path(base, "a9", "stimgate"),
    tailgate = file.path(base, "a9", "tailgate", "result.rds"),
    fbeta = file.path(base, "a9", "fbeta", "result.rds")
  )
  gs <- flowWorkspace::load_gs(ex$pathGs)
  withr::local_seed(5)
  .quiet(env$.acsCytofGateStim(gs, paths9$stimgate, unname(ex$batchList)))
  hash9 <- tools::md5sum(list.files(
    c(ex$pathGs, paths9$stimgate),
    recursive = TRUE, full.names = TRUE
  ))

  pathsDebug <- env$.acsCytofDebugPaths("tcrgd", file.path(base, "debug"))
  comparison <- tibble::tibble(
    pop = "tcrgd", SampleID = "d1", stim = "p1", cyt = ex$marker[[1]],
    freqStimManual = 5, freqUnsManual = 1, freqBsManual = 4,
    freqBsAnalysis9 = 6, diffAnalysis9 = 2, absDiffAnalysis9 = 2
  )
  htmlKeys <- tibble::tibble(
    pop = "tcrgd", SampleID = "d2", stim = "p1", cyt = ex$marker[[2]],
    htmlReason = "random"
  )
  withr::local_seed(5)
  res <- .quiet(env$.acsCytofDebugRunPopulation(
    pop = "tcrgd", paths9 = paths9, pathsDebug = pathsDebug,
    settings = list(biasUns = c(tcrgd = NA)),
    select = list(cyt = ex$marker[[1]]),
    comparison = comparison, htmlKeys = htmlKeys,
    manifest = list(run_id = "test")
  ))
  expect_true(res$success)
  # Analysis 9's GatingSet and project are untouched.
  expect_identical(
    tools::md5sum(list.files(
      c(ex$pathGs, paths9$stimgate),
      recursive = TRUE, full.names = TRUE
    )),
    hash9
  )
  sm <- readRDS(pathsDebug$summary)
  expect_identical(sm$manifest$run_id, "test")
  tbl <- sm$table
  expect_true(all(tbl$matchesAnalysis9))
  expect_true(all(is.na(tbl$tailgateGate)))
  # Pages for the selected cytokine plus the report's sample.
  paged <- tbl[!is.na(tbl$page), ]
  expect_setequal(
    paste(paged$SampleID, paged$cyt),
    c(paste(c("d1", "d2"), ex$marker[[1]]), paste("d2", ex$marker[[2]]))
  )
  expect_true(all(is.na(paged$pageError)))
  expect_true(all(file.exists(file.path(pathsDebug$pages, unique(paged$pdf)))))
  expect_equal(sort(paged$page[paged$cyt == ex$marker[[1]]]), 1:2)
  row <- tbl[tbl$SampleID == "d1" & tbl$cyt == ex$marker[[1]], ]
  expect_equal(row$freqBsManual, 4)
  expect_true(is.finite(row$propBsFinal))
  # Records are deleted, except the report's, and the project is removed.
  expect_false(dir.exists(pathsDebug$records))
  expect_false(dir.exists(pathsDebug$stimgate))
  html <- list.files(pathsDebug$html, full.names = TRUE)
  expect_length(html, 1L)

  x <- readRDS(html)
  expect_null(x$rec$truth)
  expect_null(x$rec$inputs$exTblUnsBias)
  expect_false("fdp" %in% names(env$.simDebugLocSummary(x$rec)))
  plots <- env$.simDebugLocPlots(x$rec, extraLines = env$.acsCytofDebugLines(x$row))
  expect_false("truth" %in% names(plots))
  expect_true("final gate" %in% as.character(attr(plots, "lines")$line))
  info <- env$.acsCytofDebugInfo(x$rec, x$row)
  expect_true(all(names(env$.acsCytofDebugHeadings) %in% names(info)))
  expect_identical(
    info$final$value[info$final$name == "same gates as Analysis 9"], "yes"
  )
  expect_s3_class(env$.acsCytofDebugFigure(x$rec, x$row), "ggplot")
})

test_that("HTML samples cover the largest, closest and random errors", {
  env <- .load_acs_debug_env()
  cand <- tidyr::expand_grid(
    pop = "b", stim = c("p1", "mtb"), cyt = "IL2",
    SampleID = paste0("s", 1:6)
  )
  comparison <- dplyr::mutate(cand, absDiffAnalysis9 = rep(c(5, 1, 3, 0.1, 2, 4), 2))
  comparison$absDiffAnalysis9[comparison$stim == "mtb"][6] <- NA
  withr::local_seed(3)
  seed_before <- .Random.seed
  sel <- env$.acsCytofDebugSelectHtml(cand, comparison, 2L, 1L, 1L, seed = 9L)
  expect_identical(.Random.seed, seed_before)
  p1 <- sel[sel$stim == "p1", ]
  expect_identical(p1$SampleID[p1$htmlReason == "largest error"], c("s1", "s6"))
  expect_identical(p1$SampleID[p1$htmlReason == "close agreement"], "s4")
  expect_length(p1$SampleID[p1$htmlReason == "random"], 1L)
  expect_false(any(duplicated(sel[c("stim", "SampleID")])))
  expect_identical(
    sel, env$.acsCytofDebugSelectHtml(cand, comparison, 2L, 1L, 1L, seed = 9L)
  )
  none <- env$.acsCytofDebugSelectHtml(cand, NULL, 2L, 1L, 1L, seed = 9L)
  expect_true(all(none$htmlReason == "random (no manual comparison)"))
  expect_equal(nrow(none), 8L)
  req <- env$.acsCytofDebugSelectHtml(cand, comparison, samples = "s2")
  expect_true(all(req$SampleID == "s2" & req$htmlReason == "requested"))
  expect_equal(nrow(env$.acsCytofDebugSelectHtml(cand[0, ], comparison)), 0L)
})

test_that("the manual-matched value reproduces the manual frequency", {
  env <- .load_acs_debug_env()
  xUns <- c(0, 1, 2, 3)
  xStim <- c(0, 1, 2, 3, 4, 5, 6, 7)
  x <- env$.acsCytofDebugManualMatch(xStim, xUns, 0.25)
  expect_equal(mean(xStim > x) - mean(xUns > x), 0.25)
  expect_true(is.na(env$.acsCytofDebugManualMatch(xStim, xUns, 0)))
  expect_true(is.na(env$.acsCytofDebugManualMatch(xStim, xUns, 0.9)))
  expect_identical(env$.acsCytofDebugList("all"), NULL)
  expect_identical(env$.acsCytofDebugList("none"), character(0L))
  expect_identical(env$.acsCytofDebugList("IL2, TNF"), c("IL2", "TNF"))
})

test_that("Analysis 13 shares Analysis 9's settings and never writes its caches", {
  env <- .load_acs_debug_env()
  settings <- env$.acsCytofDebugSettings(qmd9_path)
  lines9 <- readLines(qmd9_path, warn = FALSE)
  start <- which(lines9 == "#| label: scientific-settings")
  end <- start + which(lines9[-seq_len(start)] == "```")[1L]
  ref <- new.env()
  ref$run_preprocessing <- ref$run_stimgate <- ref$run_comparators <- FALSE
  eval(parse(text = lines9[(start + 1L):(end - 1L)]), ref)
  expect_identical(settings$popVec, ref$pop_vec)
  expect_identical(settings$bwMtd, ref$stimgate_bw_mtd)
  expect_identical(settings$bwScope, ref$stimgate_bw_scope)
  expect_identical(settings$locThresholdMethod, ref$stimgate_loc_threshold_method)
  expect_identical(settings$biasUns, ref$bias_uns_vec_by_pop)
  expect_identical(settings$seed, ref$analysis_seed)

  # Separate cache folders.
  p9 <- env$.acsCytofPopulationPaths("b", "fcs", "cache/acs_cytof/gs", "cache/acs_cytof/scratch")
  p13 <- env$.acsCytofDebugPaths("b", "cache/acs_cytof_debug")
  for (path in unlist(p13)) {
    expect_false(startsWith(path, "cache/acs_cytof/"))
    expect_false(any(startsWith(path, unlist(p9))))
  }
  body13 <- paste(deparse(body(env$.acsCytofDebugRunPopulation)), collapse = "\n")
  expect_true(grepl("backend_readonly = TRUE", body13, fixed = TRUE))
  expect_true(grepl(".acsCytofGateStim(", body13, fixed = TRUE))
  expect_true(grepl("pathsDebug$dir", body13, fixed = TRUE))

  lines <- readLines(qmd13_path, warn = FALSE)
  starts <- which(startsWith(lines, "```{r"))
  code <- unlist(lapply(starts, function(s) {
    e <- s + which(lines[(s + 1L):length(lines)] == "```")[1L]
    lines[(s + 1L):(e - 1L)]
  }))
  expect_no_error(parse(text = code))
  code <- paste(code, collapse = "\n")
  # No Analysis 9 writer, and Analysis 9's paths are resolved without creating.
  for (writer in c(
    ".acsCytofRunPopulation", "create_gatingset", ".acsCytofManualSave",
    ".acsCytofRunComparisonMethods", "comp_against_manual_cyt", "save_gs",
    ".write_rds_atomic"
  )) {
    expect_false(grepl(writer, code, fixed = TRUE), info = writer)
  }
  calls <- regmatches(code, gregexpr(
    "projr::projr_path_get\\([^)]*\\)", code
  ))[[1L]]
  expect_length(calls, 3L)
  expect_true(all(grepl("\"acs_cytof\"", calls, fixed = TRUE)))
  expect_true(all(grepl("create = FALSE", calls, fixed = TRUE)))
  expect_true(grepl("\"acs_cytof_debug\"", code, fixed = TRUE))
  expect_true(grepl(".acsCytofDebugSettings(", code, fixed = TRUE))
  expect_true(grepl("seed = settings$seed", code, fixed = TRUE))
  expect_true(grepl(".acsCytofMapPopulations(", code, fixed = TRUE))
})

test_that("Analysis 13 chunks re-gate, publish pages and print chosen figures", {
  env <- .load_acs_debug_env()
  ex <- getExampleData()
  withr::defer(unlink(dirname(ex$pathGs), recursive = TRUE))
  .mock_acs_example(env, ex)
  root <- .local_projr_root()
  base <- withr::local_tempdir()
  lines <- readLines(qmd13_path, warn = FALSE)
  chunk <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1L]
    parse(text = lines[(start + 1L):(end - 1L)])
  }
  env$root_dir <- root
  env$fig_key <- c("13-real-debug-acs-cytof", "quick")
  env$analysis_qmd <- "analysis/13-real-debug-acs-cytof.qmd"
  env$analysis_quick <- FALSE
  env$output_group <- NULL
  env$settings <- list(
    popVec = "tcrgd", bwMtd = "nrd0", bwScope = "cytokine",
    locThresholdMethod = "region", biasUns = c(tcrgd = NA), seed = 4L
  )
  env$run_pops <- "tcrgd"
  env$paths9 <- list(tcrgd = list(
    gs = ex$pathGs, stimgate = file.path(base, "a9", "stimgate"),
    tailgate = file.path(base, "none"), fbeta = file.path(base, "none")
  ))
  env$paths_debug <- list(
    tcrgd = env$.acsCytofDebugPaths("tcrgd", file.path(base, "debug"))
  )
  env$path_comparison <- file.path(base, "missing.rds")
  env$sel_stims <- env$sel_samples <- env$html_stims <- env$html_samples <- NULL
  env$sel_cyts <- ex$marker[[1]]
  env$html_pops <- "tcrgd"
  env$html_cyts <- ex$marker[[1]]
  env$html_n_largest <- 1L
  env$html_n_random <- 0L
  env$html_n_closest <- 0L
  env$html_seed <- 1L
  env$n_workers <- 1L

  env$run_simulations <- TRUE
  .quiet(eval(chunk("run-regating"), envir = env))
  expect_true(file.exists(env$paths_debug$tcrgd$summary))
  env$run_plots <- TRUE
  .quiet(eval(chunk("read-summaries"), envir = env))
  expect_equal(nrow(env$summary_tbl), 2L * length(ex$chnl))
  # As in a render, so figures are embedded at their saved size.
  withr::local_options(knitr.in.progress = TRUE)
  out <- paste(utils::capture.output(suppressMessages(
    for (label in c(
      "population-summary", "combination-summary", "publish-pages",
      "html-figures"
    )) {
      eval(chunk(label), envir = env)
    }
  )), collapse = "\n")
  expect_match(out, "tcrgd**: 4 tubes x cytokines", fixed = TRUE)
  expect_match(out, "random (no manual comparison)", fixed = TRUE)
  expect_match(out, "data:image/png;base64", fixed = TRUE)
  expect_true(file.exists(.projr_output_path(
    root, "table", "13-real-debug-acs-cytof", "quick", "per-sample-summary.csv"
  )))
  expect_true(file.exists(.projr_output_path(
    root, "fig", "13-real-debug-acs-cytof", "quick", "per_sample", "tcrgd",
    "p1", paste0(ex$marker[[1]], ".pdf")
  )))

  # With plots off, nothing is read.
  env$run_plots <- FALSE
  env$summary_tbl <- NULL
  .quiet(eval(chunk("read-summaries"), envir = env))
  expect_equal(nrow(env$summary_tbl), 0L)
})
