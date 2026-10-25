root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")

test_that(".simBandwidthBwOne agrees numerically with stimgate:::.bwCalcOne", {
  withr::local_preserve_seed()
  for (f in c(script_misc, script_bw)) {
    if (!file.exists(f)) stop("Expected analysis helper not found: ", f)
  }

  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)

  set.seed(42)
  x_sample <- rnorm(200, mean = 5, sd = 1.5)

  bw_params <- list(
    x = x_sample,
    bwMtd = "nrd0",
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L
  )

  bw_wrapper_params <- c(
    bw_params,
    list(bwMin = 1e-10, bwMax = 1e10, bwFallback = NULL)
  )

  set.seed(99)
  bw_direct <- do.call(env$.simBandwidthBwOne, bw_wrapper_params)

  set.seed(99)
  bw_pkg <- do.call(stimgate:::.bwCalcOne, bw_params)

  expect_equal(as.numeric(bw_direct), as.numeric(bw_pkg))
})

test_that(".simBandwidthBwOne resolves the current stimgate namespace in future workers", {
  # Fixed expected values need R's default generators; an earlier test
  # can leave the parallel (L'Ecuyer) generator selected.
  withr::local_seed(1L, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  withr::local_preserve_seed()
  skip_if_not_installed("future")
  skip_if_not_installed("pkgload")

  root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

  x_sample <- rnorm(200, mean = 5, sd = 1.5)
  bw_params <- list(
    x = x_sample,
    bwMtd = "nrd0",
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L
  )
  bw_wrapper_params <- c(
    bw_params,
    list(bwMin = 1e-10, bwMax = 1e10, bwFallback = NULL)
  )

  set.seed(99)
  bw_pkg <- do.call(stimgate:::.bwCalcOne, bw_params)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 2L)

  result <- future::future(
    {
      env <- new.env(parent = globalenv())
      source(script_misc, local = env)
      source(script_bw, local = env)
      env$.simBandwidthEnsureCurrentCheckout(root_dir)
      ns <- asNamespace("stimgate")

      suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
      set.seed(99)
      bw_worker <- do.call(env$.simBandwidthBwOne, bw_wrapper_params)

      set.seed(99)
      bw_worker_pkg <- do.call(
        get(".bwCalcOne", mode = "function", envir = ns),
        bw_params
      )

      list(
        ns_path = normalizePath(getNamespaceInfo(ns, "path"), winslash = "/", mustWork = FALSE),
        bw_worker = bw_worker,
        bw_pkg = bw_worker_pkg
      )
    },
    seed = TRUE
  )

  worker_out <- future::value(result)

  expect_equal(
    normalizePath(worker_out$ns_path, winslash = "/", mustWork = FALSE),
    normalizePath(root_dir, winslash = "/", mustWork = FALSE)
  )
  expect_true(is.finite(as.numeric(worker_out$bw_worker)))
  expect_equal(as.numeric(worker_out$bw_worker), as.numeric(worker_out$bw_pkg))
  expect_equal(as.numeric(worker_out$bw_worker), as.numeric(bw_pkg))
})


test_that("bandwidth simulation wrappers reject retired methods before fallback", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)
  for (method in paste0(c("hpi0", "hpi1", "hpi2", "hpi3", "nrd0", "sj"), "Norm")) {
    expect_error(env$.simBandwidthBwOne(
      x = numeric(), bwMtd = method, bwMin = 0, bwMax = Inf,
      bwAdj = 1, bwNcellMin = NULL, bwNcellMax = NULL, bwFallback = 0.2
    ), "Unsupported bandwidth method")
    expect_error(env$.simBandwidthEstBwDirect(bwMtd = method),
      "Unsupported bandwidth method")
    expect_error(env$.simBandwidthBsFreq(bwMtd = method),
      "Unsupported bandwidth method")
  }
})

test_that("retired bandwidth options fail as unused arguments", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_bw, local = env)
  for (option in c(
    "bwAdaptive", "bwAdaptiveCore", "bwAdaptiveExtra", "bwAdaptiveCrossover",
    "bwAdaptiveTransitionWidth", "bwAdaptiveDensityN", "bwAdaptivePadFrac"
  )) {
    expect_error(do.call(env$.simBandwidthBsFreq, setNames(list(TRUE), option)),
      "unused argument")
  }
  for (wrapper in c(".simBandwidthBwOne", ".simBandwidthEstBwDirect")) {
    for (option in c(
      "normPeakMinRel", "normExtraFrac", "normExtraMax", "normLambda",
      "normDensityN", "normExcessBwMtd", "normExcessNcell",
      "normAdaptiveNcell", "normMtd", "adaptive"
    )) {
      expect_error(do.call(env[[wrapper]], setNames(list(TRUE), option)),
        "unused argument", info = paste(wrapper, option))
    }
  }
})
