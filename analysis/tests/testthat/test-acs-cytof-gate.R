root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
script_helper <- file.path(root_dir, "scripts", "r", "acs_cytof-helper.R")
script_gate <- file.path(root_dir, "scripts", "r", "acs_cytof-gate.R")
qmd_path <- file.path(root_dir, "analysis", "9-real-compare-acs-cytof.qmd")
launcher_path <- file.path(
  root_dir,
  "scripts",
  "slurm",
  "dev-9-real-compare-acs-cytof.sh"
)

.load_acs_gate_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_helper, local = env)
  source(script_gate, local = env)
  env
}

test_that("ACS batches use mapped donors and stimuli regardless of file order", {
  env <- .load_acs_gate_env()
  mapped <- data.frame(
    SampleID = rep(c("a", "b"), each = 5),
    stim = rep(c("uns", "p1", "mtb", "ebv", "p4"), 2),
    ind = 1:10
  )
  expect_equal(env$.acsCytofBatchList(mapped), list(a = 1:5, b = 6:10))
  shuffled <- mapped[c(8, 5, 1, 9, 2, 7, 4, 10, 3, 6), ]
  shuffled$ind <- seq_len(nrow(shuffled))
  batches <- env$.acsCytofBatchList(shuffled)
  expect_equal(shuffled$stim[batches$a], c("uns", "p1", "mtb", "ebv", "p4"))
  expect_true(all(shuffled$SampleID[batches$a] == "a"))
  expect_error(env$.acsCytofBatchList(mapped[-3, ]), "Missing or duplicate")
  mapped$stim[3] <- "uns"
  expect_error(env$.acsCytofBatchList(mapped), "Missing or duplicate")
  expect_error(env$.acsCytofValidateSampleCount(19L), "multiple of five")
  expect_error(env$.acsCytofValidateSampleCount(0L), "at least 5")
})

test_that("the TCRgd tester and full population use separate output paths", {
  env <- .load_acs_gate_env()

  full <- env$.acsCytofPopulationPaths(
    pop = "tcrgd",
    pathFcsBase = "raw",
    pathGsBase = "gs",
    pathScratchBase = "scratch"
  )
  tester <- env$.acsCytofPopulationPaths(
    pop = "tcrgd",
    pathFcsBase = "raw",
    pathGsBase = "gs",
    pathScratchBase = "scratch",
    outputGroup = "tester"
  )

  expect_equal(full$fcs, tester$fcs)
  expect_equal(full$gs, file.path("gs", "tcrgd"))
  expect_equal(tester$gs, file.path("gs", "tester", "tcrgd"))
  expect_equal(
    tester$stimgate,
    file.path("scratch", "tester", "tcrgd", "stimgate")
  )
  expect_equal(
    tester$tailgate,
    file.path("scratch", "tester", "tcrgd", "tailgate", "result.rds")
  )
  expect_equal(
    tester$fbeta,
    file.path("scratch", "tester", "tcrgd", "fbeta", "result.rds")
  )
  expect_false(identical(full$stimgate, tester$stimgate))
})

test_that("the reusable population runner preserves the ACS StimGate contract", {
  env <- .load_acs_gate_env()
  runner_formals <- names(formals(env$.acsCytofRunPopulation))
  runner_body <- paste(deparse(body(env$.acsCytofRunPopulation)), collapse = "\n")

  expect_true(all(c(
    "pop", "runPreprocessing", "runMethods", "runPlots",
    "nSample", "biasUns", "outputGroup"
  ) %in% runner_formals))
  expect_true(grepl('gateCombn = "min"', runner_body, fixed = TRUE))
  expect_true(grepl('bwMtd = "nrd0"', runner_body, fixed = TRUE))
  expect_true(grepl("biasUnsFactor = biasUnsFactor", runner_body, fixed = TRUE))
  expect_null(formals(env$.acsCytofRunPopulation)$biasUns)
  expect_identical(formals(env$.acsCytofRunPopulation)$biasUnsFactor, 4)
  expect_true(grepl("calcCytPosGates = TRUE", runner_body, fixed = TRUE))
  expect_true(grepl("clusterGates = TRUE", runner_body, fixed = TRUE))
})

test_that("analysis 9 uses one runner for the tester and configured populations", {
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl(
    'pop_vec <- c("tcrgd", "cd4", "cd8", "nk", "nk_pre", "b")',
    content,
    fixed = TRUE
  ))
  expect_gte(
    lengths(regmatches(
      content,
      gregexpr(".acsCytofRunPopulation(", content, fixed = TRUE)
    )),
    1L
  )
  expect_true(grepl(".acsCytofRunPopulationSafe(", content, fixed = TRUE))
  expect_true(grepl('outputGroup = "tester"', content, fixed = TRUE))
  expect_true(grepl(".acsCytofMapPopulations(", content, fixed = TRUE))
  expect_true(grepl(
    ".acsCytofRunComparisonMethods",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "run-comparison-methods-in-parallel",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    'methods = c("stimgate", "tailgate", "fbeta")',
    content,
    fixed = TRUE
  ))
})


test_that("analysis 9 validates execution controls before running", {
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl(
    'tester_n_sample <- as.integer(.get_qmd_param_env(',
    content,
    fixed = TRUE
  ))
  expect_true(grepl(".acsCytofValidateSampleCount(tester_n_sample)", content, fixed = TRUE))
  expect_true(grepl(
    'n_workers <- as.integer(.get_qmd_param_env(',
    content,
    fixed = TRUE
  ))
  expect_true(grepl("is.na(n_workers) || n_workers < 1L", content, fixed = TRUE))
})

test_that("analysis 9 preprocessing reaches every configured population", {
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl(
    "if (isTRUE(run_preprocessing) || isTRUE(run_stimgate))",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "runPreprocessing = run_preprocessing_vec_by_pop[[pop]]",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "runPlots = run_stimgate_plots_vec_by_pop[[pop]]",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "runPreprocessingPlots = run_preprocessing_plots_vec_by_pop[[pop]]",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("analysis_seed <- 20260823L", content, fixed = TRUE))
  expect_true(grepl("set.seed(analysis_seed)", content, fixed = TRUE))
  # Both the StimGate and comparator population maps use the analysis seed.
  expect_equal(
    lengths(regmatches(
      content,
      gregexpr("seed = analysis_seed", content, fixed = TRUE)
    )),
    2L
  )
  expect_false(grepl("seed = TRUE", content, fixed = TRUE))
  expect_false(grepl(
    'Sys.setenv("STIMGATE_DEBUG" = "TRUE")',
    content,
    fixed = TRUE
  ))
})

test_that("analysis 9 does not continue after a population-stage failure", {
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl(
    'label = "ACS population preprocessing/StimGate runs"',
    content,
    fixed = TRUE
  ))
  helper_body <- paste(
    deparse(body(.load_acs_gate_env()$.acsCytofMapPopulations)),
    collapse = "\n"
  )
  expect_true(grepl("stop(", helper_body, fixed = TRUE))
  expect_false(grepl(
    '"ACS CyTOF populations failed: "',
    content,
    fixed = TRUE
  ))
})

test_that("analysis 9 Slurm launcher exports the controls the QMD reads", {
  content <- paste(readLines(launcher_path, warn = FALSE), collapse = "\n")

  expect_true(grepl(
    'run_methods_default="${RUN_METHODS:-true}"',
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    'export RUN_STIMGATE="${RUN_STIMGATE:-$run_methods_default}"',
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    'export RUN_COMPARATORS="${RUN_COMPARATORS:-$run_methods_default}"',
    content,
    fixed = TRUE
  ))
  expect_true(grepl('echo "RUN_STIMGATE: $RUN_STIMGATE"', content, fixed = TRUE))
  expect_true(grepl(
    'echo "RUN_COMPARATORS: $RUN_COMPARATORS"',
    content,
    fixed = TRUE
  ))
  expect_true(grepl("#SBATCH --ntasks=6", content, fixed = TRUE))
})

test_that("analysis 9 stages inherit simulation controls and allow overrides", {
  lines <- readLines(qmd_path, warn = FALSE)
  yaml_end <- which(lines == "---")[2L]
  yaml_params <- yaml::yaml.load(paste(lines[2:(yaml_end - 1L)], collapse = "\n"))
  expect_true(yaml_params$params$run_simulations)
  expect_false(yaml_params$params$run_plots)
  expect_false(yaml_params$execute$warning)
  expect_false(yaml_params$execute$message)
  for (flag in c("run_preprocessing", "run_stimgate", "run_comparators")) {
    expect_null(yaml_params$params[[flag]])
  }

  env <- .load_acs_gate_env()
  source(file.path(root_dir, "scripts", "r", "analysis-runtime.R"), local = env)
  withr::local_envvar(c(
    RUN_SIMULATIONS = NA, RUN_PREPROCESSING = NA, RUN_STIMGATE = NA,
    RUN_COMPARATORS = NA, RUN_PLOTS = NA
  ))
  start <- which(lines == "#| label: setup")
  end <- which(lines == "```" & seq_along(lines) > start)[1L]
  code <- parse(text = lines[seq.int(start + 1L, end - 1L)])
  controls <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      as.character(expr[[2]]) %in% c(
        "run_simulations", "run_preprocessing", "run_stimgate", "run_comparators"
      )
  }, as.list(code))
  expect_length(controls, 4L)
  stages <- c("run_preprocessing", "run_stimgate", "run_comparators")
  for (enabled in c(TRUE, FALSE)) {
    env$params <- yaml_params$params
    env$params$run_simulations <- enabled
    for (expr in controls) eval(expr, env)
    expect_identical(vapply(stages, get, logical(1), envir = env),
                     stats::setNames(rep(enabled, 3L), stages))
  }
  Sys.setenv(RUN_SIMULATIONS = "true", RUN_PREPROCESSING = "false")
  for (expr in controls) eval(expr, env)
  expect_false(env$run_preprocessing)
  expect_true(env$run_stimgate)
  expect_true(env$run_comparators)
})

test_that("a replaced directory keeps its last good content on failure", {
  env <- .load_acs_gate_env()
  path_dir <- tempfile("acs-replace-")
  withr::defer(unlink(
    list.files(dirname(path_dir), pattern = basename(path_dir), full.names = TRUE),
    recursive = TRUE
  ))

  env$.acsCytofReplaceDir(path_dir, function(path_tmp) {
    writeLines("v1", file.path(path_tmp, "out.txt"))
  })
  expect_equal(readLines(file.path(path_dir, "out.txt")), "v1")

  expect_error(
    env$.acsCytofReplaceDir(path_dir, function(path_tmp) {
      writeLines("partial", file.path(path_tmp, "out.txt"))
      stop("gating failed")
    }),
    "gating failed"
  )
  expect_equal(readLines(file.path(path_dir, "out.txt")), "v1")
  expect_equal(
    list.files(dirname(path_dir), pattern = paste0("^", basename(path_dir))),
    basename(path_dir)
  )

  env$.acsCytofReplaceDir(path_dir, function(path_tmp) {
    writeLines("v2", file.path(path_tmp, "out.txt"))
  })
  expect_equal(readLines(file.path(path_dir, "out.txt")), "v2")
  expect_equal(
    list.files(dirname(path_dir), pattern = paste0("^", basename(path_dir))),
    basename(path_dir)
  )
})
