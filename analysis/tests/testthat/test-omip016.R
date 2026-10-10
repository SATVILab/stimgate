root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

.omip016TestEnv <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (f in c(
    "analysis-runtime.R", "acs_cytof-helper.R", "acs_cytof-gate.R", "sim-misc.R",
    "sim-compare-freq_bs.R", "acs_cytof-methods.R", "omip016-prepare.R", "omip016-methods.R"
  )) {
    source(file.path(root_dir, "scripts", "r", f), local = env)
  }
  env
}

test_that("OMIP-016 helpers source in dependency order", {
  expect_no_error(.omip016TestEnv())
})

test_that("Boolean gates evaluate left to right or with & before |", {
  env <- .omip016TestEnv()
  g <- list(
    c(TRUE, TRUE, FALSE, FALSE),
    c(TRUE, FALSE, TRUE, FALSE),
    c(FALSE, FALSE, FALSE, TRUE)
  )
  # (G1 | G2) & !G0 versus G1 | (G2 & !G0).
  expect_identical(
    env$.omip016EvalBoolean("G1 | G2 & ! G0", g, "left"),
    c(FALSE, FALSE, TRUE, TRUE)
  )
  expect_identical(
    env$.omip016EvalBoolean("G1 | G2 & ! G0", g, "precedence"),
    c(TRUE, FALSE, TRUE, TRUE)
  )
  expect_identical(env$.omip016EvalBoolean("G0 | G2", g, "left"), c(TRUE, TRUE, FALSE, TRUE))
  expect_error(env$.omip016EvalBoolean("G0 | G5", g), "missing")
  expect_error(env$.omip016EvalBoolean("G0 ^ G1", g), "Unsupported")
})

test_that("polygon membership and FlowJo edge extension", {
  env <- .omip016TestEnv()
  v <- cbind(x = c(-50, 100, 100, -50), y = c(0, 0, 10, 10))
  x <- c(0, 50, 150, -200, -200)
  y <- c(5, 5, 5, 5, 20)
  expect_identical(env$.omip016InPolygon(x, y, v[, 1], v[, 2]), c(TRUE, TRUE, FALSE, FALSE, FALSE))
  # Only the compensated (x) axis is extended, and only at or below zero.
  ext <- env$.omip016ExtendVertices(v, axisCompensated = c(TRUE, FALSE), dataMin = c(-200, 0))
  expect_identical(ext[, 1], c(-201, 100, 100, -201))
  expect_identical(ext[, 2], v[, 2])
  expect_identical(env$.omip016InPolygon(x, y, ext[, 1], ext[, 2]), c(TRUE, TRUE, FALSE, TRUE, FALSE))
})

test_that("gate tree applies parents, polygons and sibling Booleans", {
  env <- .omip016TestEnv()
  ex <- cbind(`FSC-A` = c(1, 2, 3, 4, 50), `FITC-A` = c(1, 5, 1, 5, 5))
  ff <- flowCore::flowFrame(ex)
  pops <- data.frame(
    pop_id = 1:4, parent_id = c(0L, 1L, 1L, 1L),
    name = c("Cells", "Pos", "Low", "Both"),
    path = c("Cells", "Cells/Pos", "Cells/Low", "Cells/Both"),
    type = c("polygon", "polygon", "polygon", "boolean"),
    x_param = c("FSC-A", "<FITC-A>", "FSC-A", ""),
    y_param = c("FITC-A", "FSC-A", "FITC-A", ""),
    bool_expr = c("", "", "", "G0 & G1"),
    bool_refs = c("", "", "", "/Pos;/Low"),
    stringsAsFactors = FALSE
  )
  sq <- function(path, x0, x1, y0, y1) {
    data.frame(path = path, vertex = 1:4, x = c(x0, x1, x1, x0), y = c(y0, y0, y1, y1))
  }
  vertices <- rbind(
    sq("Cells", 0, 10, 0, 10), sq("Cells/Pos", 3, 100, 0, 10),
    sq("Cells/Low", 0, 2.5, 0, 10)
  )
  m <- env$.omip016ApplyGates(ff, pops, vertices)
  expect_identical(unname(m[, "Cells"]), c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(unname(m[, "Cells/Pos"]), c(FALSE, TRUE, FALSE, TRUE, FALSE))
  expect_identical(unname(m[, "Cells/Both"]), c(FALSE, TRUE, FALSE, FALSE, FALSE))
})

test_that("classification applies StimGate's cytokine-positive gates", {
  env <- .omip016TestEnv()
  x <- cbind(a = c(2, 1.5, 1.5, 0), b = c(2, 2, 0, 2))
  # Ordinary gates at 1.8; a's cytokine-positive gate 1 applies only to
  # cells ordinarily positive for b.
  out <- env$.omip016Classify(x, gate = c(1.8, 1.8), gateCyt = c(1, NA))
  expect_identical(unname(out[, "a"]), c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(unname(out[, "b"]), c(TRUE, TRUE, FALSE, TRUE))
  plain <- env$.omip016Classify(x, gate = c(1.8, Inf), gateCyt = c(NA, NA))
  expect_identical(unname(plain[, "a"]), c(TRUE, FALSE, FALSE, FALSE))
  expect_true(all(is.na(plain[, "b"])))
})

test_that("batch list puts the unstimulated tube first", {
  env <- .omip016TestEnv()
  sm <- data.frame(ind = c("1", "2", "3"), stim = c("gag", "uns", "seb"))
  expect_identical(env$.omip016BatchList(sm), list(c(2L, 1L, 3L)))
  expect_error(env$.omip016BatchList(sm[-2, ]), "unstimulated")
})

test_that("curated OMIP-016 maps are consistent", {
  env <- .omip016TestEnv()
  small <- file.path(root_dir, "_raw", "data", "small", "comparison_data", "omip016")
  sm <- env$.omip016SampleMap(small)
  mm <- env$.omip016MarkerMap(small)
  raw <- utils::read.csv(file.path(small, "raw_manifest.csv"))
  expect_setequal(sm$file, raw$file)
  expect_identical(unname(env$.omip016ResponseChannels(mm)), c("CD154", "MIP1b", "IFNg", "IL2", "TNFa"))
  beads <- mm$bead_file[nzchar(mm$bead_file)]
  expect_true(all(beads %in% sm$file[sm$role == "bead_single_stain"]))
})

test_that("the OMIP-016 workspace decodes to the expected gate tree", {
  path_jo <- tryCatch(
    withr::with_dir(root_dir, projr::projr_path_get(
      "raw-data-large", "OMIP-016", "PHY1002_20110118_ICS_SG_OMIP.jo",
      format = "absolute", create = FALSE
    )),
    error = function(e) ""
  )
  skip_if_not(file.exists(path_jo), "OMIP-016 workspace not available")
  skip_if(!nzchar(Sys.which("python3")), "python3 not available")
  env <- .omip016TestEnv()
  out <- withr::local_tempdir()
  samples <- c("13523 17012011_NS_F01.fcs", "13523 17012011_SEB_F06.fcs")
  ws <- env$.omip016ParseWorkspace(path_jo, out, samples, root_dir)
  expect_identical(dim(ws$spill), c(10L, 10L))
  expect_equal(ws$spill["FITC-A", "PE-A"], 0.23000403)
  for (s in samples) {
    pops <- ws$pops[ws$pops$sample == s, ]
    expect_identical(nrow(pops), 22L)
    cd4 <- pops[pops$path == "Single cells/Lymphocytes/Lives/CD3+/CD4+", ]
    expect_identical(cd4$type, "boolean")
    expect_identical(cd4$bool_expr, "G1 | G4 | G5 | G2 | G3 & ! G0")
    ifng <- ws$vertices[ws$vertices$sample == s &
      ws$vertices$path == "Single cells/Lymphocytes/Lives/CD3+/CD4+/IFNg+", ]
    expect_equal(min(ifng$x), 580.966674804688)
  }
})
