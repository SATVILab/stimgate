root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
qmd9a_path <- file.path(root_dir, "analysis", "9a-real-tune-acs-cytof-tailgate.qmd")

.load_tailgate_select_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (f in c(
    "analysis-runtime.R", "analysis-plot-style.R", "acs_cytof-helper.R",
    "acs_cytof-gate.R", "acs_cytof-methods.R", "acs_cytof-manual.R",
    "acs_cytof-plot_cyt.R", "acs_cytof-tailgate-select.R"
  )) {
    source(file.path(root_dir, "scripts", "r", f), local = env)
  }
  env
}

# Eight configurations whose estimates degrade with their position, over two
# eligible strata, one ineligible stratum and a p4 stratum.
.tailgate_select_table <- function(methods) {
  withr::with_seed(1, {
    grid <- expand.grid(
      SampleID = paste0("d", 1:8), pop = "CD4 T cells",
      cyt = c("IFNg", "TNF", "IL17"), stim = c("mtb", "p4"),
      stringsAsFactors = FALSE
    )
    grid$freq_uns_man <- 0.005
    grid$freq_bs_man <- ifelse(grid$cyt == "IL17", 0.001, seq(0.1, 1.5, length.out = 8))
    grid$freq_stim_man <- grid$freq_bs_man + grid$freq_uns_man
    do.call(rbind, lapply(seq_along(methods), function(k) {
      g <- grid
      g$method <- methods[[k]]
      g$freq_bs_auto <- g$freq_bs_man + stats::rnorm(nrow(g), 0, 0.02 * k)
      g$thresholdFailed <- FALSE
      g
    }))
  })
}

test_that("Analysis 9a compares four Tailgate parameter sets with and without zeros", {
  env <- .load_tailgate_select_env()
  configs <- env$.acsCytofTailgateConfigs()
  expect_identical(nrow(configs), 8L)
  expect_false(anyDuplicated(configs$method) > 0L)
  expect_setequal(paste(configs$config, configs$removeZero),
    as.vector(outer(c("simulation", "autotol", "acs", "package"), c(TRUE, FALSE), paste)))
  expect_identical(names(env$.acsCytofTailgateConfigLabels()), configs$method)
  for (i in seq_len(nrow(configs))) {
    p <- env$.acsCytofComparatorSettings(configs$method[[i]])$params
    expect_identical(p[c("autoTol", "tol", "bias", "removeZero")],
      as.list(configs[i, c("autoTol", "tol", "bias", "removeZero")]))
  }
  sim <- configs[configs$config == "simulation", ]
  expect_true(all(sim$autoTol) && all(sim$tol == 0.01) && all(sim$bias == 0.1))
  pkg <- configs[configs$config == "package", ]
  expect_true(all(!pkg$autoTol) && all(pkg$tol == 0.01) && all(pkg$bias == 0))
  acs <- configs[configs$config == "acs", ]
  expect_true(all(!acs$autoTol) && all(acs$tol == 2e-5) && all(acs$bias == 0.2))
  skip_if_not_installed("cytoUtils")
  # The "package" set is cytoUtils' own default.
  defaults <- formals(cytoUtils:::.cytokine_cutpoint)
  expect_equal(defaults$tol, 0.01)
  expect_false(defaults$auto_tol)
})

test_that("Scoring uses eligible non-p4 strata and donors gated by every configuration", {
  env <- .load_tailgate_select_env()
  methods <- env$.acsCytofTailgateSelectionMethods()
  tbl <- .tailgate_select_table(methods)
  # One failed gate for the last configuration removes that donor for all.
  failed <- tbl$method == methods[[8]] & tbl$cyt == "IFNg" & tbl$stim == "mtb" & tbl$SampleID == "d3"
  tbl$thresholdFailed[failed] <- TRUE
  strata <- env$.acsCytofTailgateScoreStrata(tbl)
  expect_false(any(strata$stim == "p4"))
  expect_setequal(unique(strata$cyt[!strata$eligible]), "IL17")
  ifng <- strata[strata$cyt == "IFNg", ]
  expect_true(all(ifng$n_common == 7L))
  expect_identical(ifng$n_failed, c(rep(0L, 7), 1L))
  rows <- tbl[tbl$method == methods[[1]] & tbl$cyt == "IFNg" & tbl$stim == "mtb" & tbl$SampleID != "d3", ]
  expect_equal(ifng$pcc[[1]], stats::cor(rows$freq_bs_auto, rows$freq_bs_man))
  expect_true(ifng$pcc_own[[1]] != ifng$pcc[[1]])

  ranking <- env$.acsCytofTailgateRanking(strata)
  expect_identical(ranking$method[[1]], methods[[1]])
  expect_true(all(ranking$n_scored == 2L) && all(ranking$n_eligible == 2L))
  expect_identical(env$.acsCytofTailgateSelectionStatus(ranking), "selected")
  expect_equal(ranking$mean_pcc[[1]], mean(strata$pcc[strata$eligible & strata$method == methods[[1]]]))
})

test_that("A stratum with an undefined PCC is left out for every configuration", {
  env <- .load_tailgate_select_env()
  methods <- env$.acsCytofTailgateSelectionMethods()
  tbl <- .tailgate_select_table(methods)
  # A constant estimate gives one configuration an undefined PCC for TNF.
  constant <- tbl$method == methods[[2]] & tbl$cyt == "TNF"
  tbl$freq_bs_auto[constant] <- 0.5
  ranking <- env$.acsCytofTailgateRanking(env$.acsCytofTailgateScoreStrata(tbl))
  expect_true(all(ranking$n_scored == 1L))
  expect_identical(ranking$n_undefined[ranking$method == methods[[2]]], 1L)
  expect_true(all(ranking$n_undefined[ranking$method != methods[[2]]] == 0L))
})

test_that("Ties and missing scores select nothing, and Analysis 9b refuses them", {
  env <- .load_tailgate_select_env()
  methods <- env$.acsCytofTailgateSelectionMethods()
  tbl <- .tailgate_select_table(methods)
  tbl$freq_bs_auto[tbl$method == methods[[2]]] <- tbl$freq_bs_auto[tbl$method == methods[[1]]]
  ranking <- env$.acsCytofTailgateRanking(env$.acsCytofTailgateScoreStrata(tbl))
  expect_identical(env$.acsCytofTailgateSelectionStatus(ranking), "tie")
  expect_identical(sum(ranking$best), 2L)
  record <- local({
    attr(tbl, "manifest") <- list(methods = list(), manualInputHash = "x")
    withr::with_dir(root_dir, env$.acsCytofTailgateSelectionRecord(ranking, tbl))
  })
  expect_identical(record$status, "tie")
  expect_setequal(record$tied, methods[1:2])
  expect_error(env$.acsCytofComparatorSettings("tailgate", record), "did not select")
  expect_error(env$.acsCytofComparatorSettings("tailgate"), "Analysis 9a")
  ranking$mean_pcc <- NA_real_
  ranking$best <- FALSE
  expect_identical(env$.acsCytofTailgateSelectionStatus(ranking), "unavailable")
})

test_that("Analysis 9b's Tailgate settings come from the recorded selection", {
  env <- .load_tailgate_select_env()
  configs <- env$.acsCytofTailgateConfigs()
  pick <- configs[configs$method == "tailgate_simulation_rmzero", ]
  selection <- list(status = "selected", method = pick$method, config = as.list(pick))
  settings <- env$.acsCytofComparatorSettings("tailgate", selection)
  expect_identical(settings$tailgateConfig, pick$method)
  expect_identical(settings$params, env$.acsCytofComparatorSettings(pick$method)$params)
  expect_equal(settings$cacheVersion, 3L)
  # A cached result from another configuration fails the settings check.
  other <- utils::modifyList(selection, list(
    method = "tailgate_acs_rmzero", config = as.list(configs[configs$method == "tailgate_acs_rmzero", ])
  ))
  expect_false(isTRUE(all.equal(settings, env$.acsCytofComparatorSettings("tailgate", other))))
  # A changed definition invalidates the saved selection.
  stale <- selection
  stale$config$bias <- 0.3
  expect_error(env$.acsCytofComparatorSettings("tailgate", stale), "no longer matches")

  manifest <- list(context = list(gitSha = "a", preprocessing = list(sampleMap = 1)))
  tbl <- data.frame(x = 1)
  attr(tbl, "manifest") <- list(methods = list(cd4 = list(stimgate = manifest)))
  selection$preprocessing <- list(cd4 = list(preprocessing = list(sampleMap = 1)))
  expect_true(env$.acsCytofCheckSelectionPreprocessing(selection, tbl))
  selection$preprocessing$cd4$preprocessing$sampleMap <- 2
  expect_error(env$.acsCytofCheckSelectionPreprocessing(selection, tbl), "differ")
})

test_that("The Typst table matches the HTML table and compiles", {
  env <- .load_tailgate_select_env()
  methods <- env$.acsCytofTailgateSelectionMethods()
  ranking <- env$.acsCytofTailgateRanking(
    env$.acsCytofTailgateScoreStrata(.tailgate_select_table(methods))
  )
  html <- env$.acsCytofTailgateHtmlTable(ranking)
  typst <- env$.acsCytofTailgateTypst(ranking)
  expect_match(typst, "^#figure\\(")
  for (value in c(html$`Mean PCC`, html$`Strata scored`, html$`Failed gates`)) {
    expect_true(grepl(value, typst, fixed = TRUE))
  }
  expect_match(typst, "$2 times 10^(-5)$", fixed = TRUE)
  expect_match(html$Configuration[[1]], "(selected)", fixed = TRUE)
  expect_identical(lengths(regmatches(typst, gregexpr("\n    \\[\\*", typst))), 1L)

  quarto <- Sys.which("quarto")
  skip_if(!nzchar(quarto), "Quarto (with Typst) is unavailable.")
  dir <- withr::local_tempdir()
  file <- file.path(dir, "table.typ")
  writeLines(c('#set page(paper: "a4", margin: (left: 1.5in, rest: 1in))', typst), file)
  status <- system2(quarto, c("typst", "compile", file, file.path(dir, "table.pdf")),
    stdout = TRUE, stderr = TRUE)
  expect_null(attr(status, "status"))
  expect_true(file.exists(file.path(dir, "table.pdf")))
})

test_that("Thesis heatmaps are drawn at the thesis text width", {
  env <- .load_tailgate_select_env()
  expect_identical(env$.acsCytofThesisFigureSize(TRUE)$width, 146.5)
  expect_identical(env$.acsCytofThesisFigureSize(FALSE)$units, "mm")
  methods <- env$.acsCytofTailgateSelectionMethods()
  strata <- env$.acsCytofTailgateScoreStrata(.tailgate_select_table(methods))
  p <- env$.acsCytofValidationPlotCorrelation(
    env$.acsCytofTailgateCorrelationTable(strata), methods[[1]], "pcc", thesis = TRUE
  )
  expect_null(p$labels$caption)
  dir <- withr::local_tempdir()
  paths <- env$.acsCytofSaveThesisFigure(p, file.path(dir, "fig"), env$.acsCytofThesisFigureSize(TRUE))
  expect_true(all(file.exists(paths)))
  png <- png::readPNG(paths[[2]])
  expect_equal(ncol(png) / 300 * 25.4, 146.5, tolerance = 0.01)
})

test_that("Analysis 9a saves thesis outputs to output-thesis and selects by PCC", {
  qmd <- paste(readLines(qmd9a_path), collapse = "\n")
  expect_match(qmd, '"output-thesis"', fixed = TRUE)
  expect_match(qmd, 'thesis_parts <- c("appendix", fig_key[-1L])', fixed = TRUE)
  expect_match(qmd, "acs-tailgate-ranking.typ", fixed = TRUE)
  expect_match(qmd, "```typst", fixed = TRUE)
  expect_match(qmd, "code-copy: true", fixed = TRUE)
  expect_match(qmd, "in-sample", fixed = TRUE)
  projr <- yaml::read_yaml(file.path(root_dir, "_projr.yml"))
  expect_identical(projr$directories[["output-thesis"]]$path, "_output/thesis")
  expect_identical(projr$directories$output$path, "_output/general")
  expect_identical(projr$build$github[["output-thesis"]]$content, "output-thesis")
})
