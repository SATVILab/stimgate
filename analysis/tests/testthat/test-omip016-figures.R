root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

.omip016FiguresTestEnv <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (f in c(
    "analysis-runtime.R", "analysis-plot-style.R", "acs_cytof-helper.R", "acs_cytof-gate.R",
    "sim-misc.R", "sim-compare-freq_bs.R", "acs_cytof-methods.R", "sim-debug-loc.R",
    "acs_cytof-debug.R", "omip016-prepare.R", "omip016-methods.R", "omip016-figures.R"
  )) {
    source(file.path(root_dir, "scripts", "r", f), local = env)
  }
  env
}

# Two tubes of CD4 T cells on the GatingSet scale (asinh fluorescence,
# linear SSC-A), manual labels from the rectangles and one gate per method.
.omip016FiguresFixture <- function(env) {
  channels <- c("FITC-A", "PE-A", "Horizon V450-A", "Alexa Fluor 700-A", "APC-A")
  markers <- c("CD154", "MIP1b", "IFNg", "TNFa", "IL2")
  thresholds <- data.frame(
    marker = markers, channel = channels,
    y_param = c("<Horizon V450-A>", "<Horizon V450-A>", "SSC-A", "SSC-A", "SSC-A"),
    x_max = 1e5, y_min = c(-100, -100, 1000, 1000, 1000),
    y_max = c(1e5, 1e5, 1e5, 1e5, 1e5),
    threshold_raw = 300, stringsAsFactors = FALSE
  )
  thresholds$threshold_trans <- env$.omip016Transform(thresholds$threshold_raw)
  withr::with_seed(1, {
    frames <- lapply(c(uns = 200L, seb = 300L), function(n) {
      ex <- cbind(matrix(stats::rnorm(n * 5L, 1.5, 1), n, 5L), stats::runif(n, 2e4, 9e4))
      colnames(ex) <- c(channels, "SSC-A")
      flowCore::flowFrame(ex)
    })
  })
  cs <- flowWorkspace::flowSet_to_cytoset(flowCore::flowSet(frames))
  flowWorkspace::sampleNames(cs) <- c("uns.fcs", "seb.fcs")
  gs <- flowWorkspace::GatingSet(cs)
  labels <- lapply(1:2, function(ind) {
    ex <- env$.omip016Expr(gs, ind, c(channels, "SSC-A"))
    out <- data.frame(event = seq_len(nrow(ex)))
    for (j in seq_along(markers)) {
      out[[markers[[j]]]] <- ex[, channels[[j]]] > thresholds$threshold_trans[[j]]
    }
    out
  })
  names(labels) <- c("uns.fcs", "seb.fcs")
  prep <- list(
    gs = gs, labels = labels,
    pre = list(
      sampleMap = data.frame(ind = c("1", "2"), file = names(labels), stim = c("uns", "seb")),
      batchList = list(1:2), responseChannels = stats::setNames(markers, channels)
    )
  )
  gates <- do.call(rbind, lapply(c("stimgate", "manual"), function(method) {
    data.frame(
      ind = 2L, chnl = channels, method = method,
      gate = if (method == "stimgate") 1.8 else thresholds$threshold_trans,
      gateCyt = NA_real_, thresholdOrigin = NA_character_, thresholdReason = NA_character_,
      thresholdFallbackUsed = FALSE, stringsAsFactors = FALSE
    )
  }))
  list(prep = prep, thresholds = thresholds, scored = env$.omip016Score(prep, gates))
}

test_that("Figure 1E gate bounds follow the analysis scale and edge extension", {
  env <- .omip016FiguresTestEnv()
  fx <- .omip016FiguresFixture(env)
  g <- env$.omip016FigureGates(fx$thresholds)
  expect_identical(g$marker, names(env$.omip016FigureMarkers))
  expect_equal(g$xmin, rep(asinh(300 / 150), 5L))
  # Compensated lower bounds at or below zero reach the panel edge;
  # SSC-A bounds stay linear.
  expect_equal(g$ymin[g$y_channel == "Horizon V450-A"], c(-Inf, -Inf))
  expect_equal(g$ymin[g$y_channel == "SSC-A"], rep(1000, 3L))
  expect_equal(g$ymax[g$y_channel == "Horizon V450-A"], rep(asinh(1e5 / 150), 2L))
})

test_that("Figure 1E panels count the manual labels and methods reproduce scored counts", {
  env <- .omip016FiguresTestEnv()
  fx <- .omip016FiguresFixture(env)
  fd <- env$.omip016FigurePanels(fx$prep, fx$thresholds, c("uns", "seb"))
  expect_equal(
    fd$reference$count,
    unlist(lapply(fx$prep$labels, function(l) colSums(l[names(env$.omip016FigureMarkers)])),
      use.names = FALSE)
  )
  expect_identical(nrow(fd$cells), 5L * (200L + 300L))
  md <- env$.omip016MethodFigureData(fx$prep, fx$thresholds, fx$scored,
    c(stimgate = "StimGate"))
  expect_identical(levels(md$reference$tube), "SEB")
  sg <- fx$scored[fx$scored$method == "stimgate", ]
  expect_equal(md$methods$frequency, 100 * sg$countStim[match(md$methods$marker, sg$cyt)] / 300)
  bad <- fx$scored
  bad$countStim[bad$method == "stimgate"][[1]] <- -1L
  expect_error(
    env$.omip016MethodFigureData(fx$prep, fx$thresholds, bad, c(stimgate = "StimGate")),
    "does not reproduce"
  )
  colours <- c(StimGate = "#0072B2")
  markerMap <- data.frame(
    channel = c(fx$thresholds$channel, "SSC-A"),
    marker = c(fx$thresholds$marker, "SSC-A"),
    fluorochrome = c("FITC", "PE", "V450", "Alexa Fluor 700", "APC", "")
  )
  expect_s3_class(env$.omip016AuthorFigure(fd, markerMap), "ggplot")
  expect_s3_class(env$.omip016MethodsFigure(md, markerMap, colours), "ggplot")
  expect_s3_class(env$.omip016UnivariateFigure(md, colours), "ggplot")
})
