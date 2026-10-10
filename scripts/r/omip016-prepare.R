# OMIP-016 (FlowRepository FR-FCM-ZZ2T) preparation for StimGate and the
# Tailgate/F-beta comparators.
#
# The raw FCS files are uncompensated (their $SPILL is the identity). The
# authors' FlowJo for Mac workspace holds the compensation matrix and the
# manual gates; scripts/python/omip016_flowjo_jo.py decodes them. Here they
# are re-applied to the raw events to give the CD4 T-cell population and the
# manual cytokine-positive labels used as the reference. The raw directory is
# only read.
#
# Requires scripts/r/acs_cytof-helper.R (for .acsCytofReplaceDir()).

# Fluorescence channels are stored as asinh(x / cofactor) in the GatingSet;
# manual gate vertices stay on the compensated linear scale.
.omip016Cofactor <- 150

.omip016Transform <- function(x, cofactor = .omip016Cofactor) {
  asinh(x / cofactor)
}

# FlowJo extends gate vertices drawn at or below this value on a compensated
# axis to the bottom of the data (CytoML's extend_val default), so events
# piled below the axis minimum stay inside gates touching that edge.
.omip016ExtendVal <- 0

.omip016Version <- 1L

.omip016Paths <- function(pathOut) {
  list(
    root = pathOut,
    workspace = file.path(pathOut, "workspace"),
    labels = file.path(pathOut, "labels"),
    fcs = file.path(pathOut, "fcs", "cd4"),
    gs = file.path(pathOut, "gs", "cd4"),
    reference = file.path(pathOut, "reference"),
    qc = file.path(pathOut, "qc"),
    manifest = file.path(pathOut, "manifest.rds")
  )
}

# Stored beside, not inside, the GatingSet folder: flowWorkspace::load_gs()
# rejects any extra file or folder in it.
.omip016PreprocessingFile <- function(pathGs) {
  paste0(pathGs, ".omip016-preprocessing.rds")
}

.omip016ReadCsv <- function(path) {
  if (!file.exists(path)) stop("OMIP-016 input not found: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.omip016SampleMap <- function(pathSmall) {
  x <- .omip016ReadCsv(file.path(pathSmall, "sample_map.csv"))
  need <- c("file", "role", "SampleID", "stim", "tube_name", "well_id")
  if (!all(need %in% names(x))) {
    stop("sample_map.csv needs columns: ", paste(need, collapse = ", "))
  }
  if (anyDuplicated(x$file)) stop("sample_map.csv repeats a file.")
  cells <- x[x$role == "cells", , drop = FALSE]
  if (sum(cells$stim == "uns") != 1L || nrow(cells) < 2L ||
      anyDuplicated(cells$stim) || length(unique(cells$SampleID)) != 1L) {
    stop("OMIP-016 cell samples must be one donor with one unstimulated tube and distinct stimuli.")
  }
  x
}

.omip016MarkerMap <- function(pathSmall) {
  x <- .omip016ReadCsv(file.path(pathSmall, "marker_map.csv"))
  need <- c("channel", "marker", "fluorochrome", "role", "response", "manual_gate", "bead_file")
  if (!all(need %in% names(x))) {
    stop("marker_map.csv needs columns: ", paste(need, collapse = ", "))
  }
  x$response <- as.logical(x$response)
  if (anyNA(x$response) || anyDuplicated(x$channel)) {
    stop("marker_map.csv needs one row per channel and a logical 'response'.")
  }
  x
}

# Response markers in gating order, named by channel.
.omip016ResponseChannels <- function(markerMap) {
  x <- markerMap[markerMap$response, , drop = FALSE]
  stats::setNames(x$marker, x$channel)
}

# The unstimulated tube first, as gateStim() requires; then the other tubes
# in sample-map order. Indices refer to the GatingSet's sample order.
.omip016BatchList <- function(sampleMapGs) {
  if (!all(c("ind", "stim") %in% names(sampleMapGs))) {
    stop("The GatingSet sample map needs 'ind' and 'stim'.")
  }
  ind <- as.integer(sampleMapGs$ind)
  isUns <- sampleMapGs$stim == "uns"
  if (sum(isUns) != 1L) stop("Expected exactly one unstimulated tube.")
  list(c(ind[isUns], ind[!isUns]))
}

.omip016ValidateRaw <- function(pathRaw, pathSmall) {
  expected <- .omip016ReadCsv(file.path(pathSmall, "raw_manifest.csv"))
  path <- file.path(pathRaw, expected$file)
  exists <- file.exists(path)
  md5 <- rep(NA_character_, length(path))
  md5[exists] <- unname(tools::md5sum(path[exists]))
  size <- rep(NA_real_, length(path))
  size[exists] <- file.size(path[exists])
  out <- data.frame(
    file = expected$file,
    exists = exists,
    size = size,
    size_ok = exists & size == expected$bytes,
    md5 = md5,
    md5_ok = exists & md5 == expected$md5,
    stringsAsFactors = FALSE
  )
  extra <- setdiff(list.files(pathRaw), expected$file)
  if (!all(out$md5_ok & out$size_ok)) {
    bad <- out$file[!(out$md5_ok & out$size_ok)]
    stop("OMIP-016 raw files missing or changed: ", paste(bad, collapse = ", "))
  }
  attr(out, "extra") <- extra
  out
}

.omip016Python <- function() {
  python <- Sys.getenv("OMIP016_PYTHON", unset = "")
  if (!nzchar(python)) python <- Sys.which("python3")
  if (!nzchar(python)) stop("python3 is required to decode the OMIP-016 workspace.")
  unname(python)
}

.omip016ParseWorkspace <- function(pathJo, pathOut, samples, pathRoot = ".") {
  script <- file.path(pathRoot, "scripts", "python", "omip016_flowjo_jo.py")
  if (!file.exists(script)) stop("Workspace parser not found: ", script)
  out <- system2(
    .omip016Python(),
    c("-I", shQuote(script), shQuote(pathJo), shQuote(pathOut), shQuote(samples)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("Decoding the OMIP-016 workspace failed:\n", paste(out, collapse = "\n"))
  }
  .omip016ReadWorkspace(pathOut)
}

.omip016ReadWorkspace <- function(path) {
  comp <- .omip016ReadCsv(file.path(path, "compensation.csv"))
  spill <- as.matrix(comp[, -1, drop = FALSE])
  rownames(spill) <- comp$channel
  if (!identical(rownames(spill), colnames(spill))) {
    stop("Compensation matrix rows and columns differ.")
  }
  list(
    spill = spill,
    pops = .omip016ReadCsv(file.path(path, "populations.csv")),
    vertices = .omip016ReadCsv(file.path(path, "vertices.csv")),
    meta = jsonlite::read_json(file.path(path, "workspace.json"))
  )
}

.omip016ReadFcs <- function(path) {
  flowCore::read.FCS(
    path,
    transformation = FALSE,
    truncate_max_range = FALSE,
    emptyValue = FALSE
  )
}

.omip016Compensate <- function(ff, spill) {
  if (!all(colnames(spill) %in% flowCore::colnames(ff))) {
    stop("Compensation channels are missing from the FCS file.")
  }
  flowCore::compensate(ff, spill)
}

# "<FITC-A>" denotes the compensated FITC-A parameter in the workspace.
.omip016ParamChannel <- function(param) {
  sub("^<(.*)>$", "\\1", param)
}

# Move vertices at or below `extendVal` on compensated axes to below the
# data minimum, mirroring FlowJo/CytoML handling of gates drawn to the edge.
.omip016ExtendVertices <- function(vertices, axisCompensated, dataMin,
                                   extendVal = .omip016ExtendVal) {
  for (j in seq_len(2L)) {
    if (isTRUE(axisCompensated[[j]])) {
      low <- vertices[, j] <= extendVal
      vertices[low, j] <- min(dataMin[[j]], min(vertices[, j])) - 1
    }
  }
  vertices
}

# Point-in-polygon by ray casting (even-odd rule), vectorised over events.
.omip016InPolygon <- function(x, y, vx, vy) {
  n <- length(vx)
  inside <- logical(length(x))
  j <- n
  for (i in seq_len(n)) {
    crosses <- ((vy[[i]] > y) != (vy[[j]] > y)) &
      (x < (vx[[j]] - vx[[i]]) * (y - vy[[i]]) / (vy[[j]] - vy[[i]]) + vx[[i]])
    inside <- xor(inside, crosses)
    j <- i
  }
  inside
}

# Evaluate a FlowJo Boolean expression over operands G0, G1, ... (logical
# vectors). `!` binds to the next operand. "left" applies & and | from left
# to right; "precedence" binds & before |.
.omip016EvalBoolean <- function(expr, operands, mode = c("left", "precedence")) {
  mode <- match.arg(mode)
  tokens <- regmatches(expr, gregexpr("G[0-9]+|[|&!]", expr))[[1]]
  if (!identical(paste(tokens, collapse = ""), gsub("[[:space:]]", "", expr))) {
    stop("Unsupported Boolean expression: ", expr)
  }
  terms <- list()
  ops <- character()
  negate <- FALSE
  for (tok in tokens) {
    if (tok == "!") {
      negate <- !negate
    } else if (tok %in% c("|", "&")) {
      ops <- c(ops, tok)
    } else {
      k <- as.integer(sub("G", "", tok)) + 1L
      if (k > length(operands)) stop("Boolean operand ", tok, " is missing.")
      val <- operands[[k]]
      terms[[length(terms) + 1L]] <- if (negate) !val else val
      negate <- FALSE
    }
  }
  if (length(ops) != length(terms) - 1L) stop("Malformed Boolean expression: ", expr)
  if (mode == "left") {
    out <- terms[[1]]
    for (i in seq_along(ops)) {
      out <- if (ops[[i]] == "&") out & terms[[i + 1L]] else out | terms[[i + 1L]]
    }
    return(out)
  }
  # & first: split into |-separated groups of &-joined terms.
  groups <- split(seq_along(terms), cumsum(c(TRUE, ops == "|")))
  Reduce(`|`, lapply(groups, function(g) Reduce(`&`, terms[g])))
}

# Membership of every workspace population for one compensated frame:
# a logical matrix with one column per population path.
.omip016ApplyGates <- function(ffComp, pops, vertices,
                               booleanMode = "left", extend = TRUE) {
  ex <- flowCore::exprs(ffComp)
  pops <- pops[order(pops$pop_id), , drop = FALSE]
  member <- matrix(
    FALSE, nrow(ex), nrow(pops),
    dimnames = list(NULL, pops$path)
  )
  byId <- stats::setNames(seq_len(nrow(pops)), pops$pop_id)
  # Polygons first: a Boolean gate refers to siblings stored after it.
  inside <- list()
  for (i in which(pops$type == "polygon")) {
    v <- vertices[vertices$path == pops$path[[i]], , drop = FALSE]
    v <- as.matrix(v[order(v$vertex), c("x", "y")])
    params <- c(pops$x_param[[i]], pops$y_param[[i]])
    chnl <- .omip016ParamChannel(params)
    if (!all(chnl %in% colnames(ex))) {
      stop("Gate '", pops$path[[i]], "' uses missing channels.")
    }
    if (isTRUE(extend)) {
      v <- .omip016ExtendVertices(
        v, axisCompensated = grepl("^<.*>$", params),
        dataMin = c(min(ex[, chnl[[1]]]), min(ex[, chnl[[2]]]))
      )
    }
    inside[[pops$path[[i]]]] <- .omip016InPolygon(
      ex[, chnl[[1]]], ex[, chnl[[2]]], v[, 1], v[, 2]
    )
  }
  resolve <- function(i) {
    parent <- pops$parent_id[[i]]
    parentIn <- if (parent == 0L) rep(TRUE, nrow(ex)) else member[, byId[[as.character(parent)]]]
    if (pops$type[[i]] == "polygon") {
      return(parentIn & inside[[pops$path[[i]]]])
    }
    refs <- strsplit(pops$bool_refs[[i]], ";", fixed = TRUE)[[1]]
    prefix <- if (parent == 0L) "" else paste0(pops$path[[byId[[as.character(parent)]]]], "/")
    operands <- lapply(paste0(prefix, sub("^/", "", refs)), function(p) {
      k <- match(p, pops$path)
      if (is.na(k) || pops$type[[k]] != "polygon") {
        stop("Boolean reference '", p, "' must be a sibling polygon gate.")
      }
      parentIn & inside[[p]]
    })
    parentIn & .omip016EvalBoolean(pops$bool_expr[[i]], operands, booleanMode)
  }
  # Parents precede children in the pre-order table.
  for (i in seq_len(nrow(pops))) member[, i] <- resolve(i)
  member
}

# Manual cytokine gates as 1D thresholds: the lower x bound of each
# rectangle drawn on the response channel, and whether the other bounds can
# exclude any CD4 T cell above it (then the gate is not a pure threshold).
.omip016ManualThresholds <- function(pops, vertices, markerMap, pathCd4) {
  resp <- markerMap[markerMap$response, , drop = FALSE]
  out <- lapply(seq_len(nrow(resp)), function(k) {
    path <- paste0(pathCd4, "/", resp$manual_gate[[k]])
    i <- match(path, pops$path)
    if (is.na(i) || pops$type[[i]] != "polygon" ||
        .omip016ParamChannel(pops$x_param[[i]]) != resp$channel[[k]]) {
      stop("No manual gate on ", resp$channel[[k]], " at ", path, ".")
    }
    v <- vertices[vertices$path == path, , drop = FALSE]
    data.frame(
      marker = resp$marker[[k]],
      channel = resp$channel[[k]],
      gate_path = path,
      n_vertices = nrow(v),
      x_min = min(v$x), x_max = max(v$x),
      y_param = pops$y_param[[i]], y_min = min(v$y), y_max = max(v$y),
      threshold_raw = min(v$x),
      threshold_trans = .omip016Transform(min(v$x)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# Single-stain bead estimate of the spillover matrix, for comparison with
# the workspace matrix only: the median of beads above the 99th percentile of
# the unstained beads in the primary channel, minus the unstained median.
.omip016BeadSpillover <- function(pathRaw, markerMap, sampleMap, channels) {
  negFile <- sampleMap$file[sampleMap$role == "bead_unstained"]
  if (length(negFile) != 1L) stop("Expected one unstained bead file.")
  neg <- flowCore::exprs(.omip016ReadFcs(file.path(pathRaw, negFile)))[, channels]
  negMed <- apply(neg, 2, stats::median)
  beads <- markerMap[nzchar(markerMap$bead_file) & markerMap$channel %in% channels, ]
  rows <- lapply(seq_len(nrow(beads)), function(k) {
    ex <- flowCore::exprs(.omip016ReadFcs(file.path(pathRaw, beads$bead_file[[k]])))[, channels]
    primary <- beads$channel[[k]]
    pos <- ex[ex[, primary] > stats::quantile(neg[, primary], 0.99), , drop = FALSE]
    signal <- apply(pos, 2, stats::median) - negMed
    data.frame(
      bead_file = beads$bead_file[[k]],
      fluorochrome_channel = primary,
      detector = channels,
      n_positive_beads = nrow(pos),
      spill_beads = unname(signal / signal[[primary]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.omip016WriteCsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

.omip016WriteCd4Fcs <- function(ffComp, keep, path, file, wsName) {
  ff <- ffComp[keep, ]
  kw <- flowCore::keyword(ff)
  kw[["OMIP016_SOURCE_FILE"]] <- file
  kw[["OMIP016_POPULATION"]] <- "Single cells/Lymphocytes/Lives/CD3+/CD4+"
  kw[["OMIP016_COMPENSATION"]] <- paste0("applied: ", wsName)
  kw[["OMIP016_TRANSFORM"]] <- "none (compensated linear values)"
  flowCore::keyword(ff) <- kw
  flowCore::write.FCS(ff, path)
  invisible(path)
}

.omip016PopulationCounts <- function(member, file) {
  n <- colSums(member)
  pathVec <- colnames(member)
  parent <- sub("/[^/]*$", "", pathVec)
  parentN <- ifelse(grepl("/", pathVec), n[match(parent, pathVec)], nrow(member))
  data.frame(
    file = file, path = pathVec, count = unname(n),
    parent_count = unname(parentN),
    pct_parent = 100 * unname(n) / unname(parentN),
    stringsAsFactors = FALSE
  )
}

# Prepare everything into `pathOut` (built in a temporary sibling and
# swapped in only on success).
.omip016Prepare <- function(pathRaw, pathSmall, pathOut, pathRoot = ".",
                            writePlots = TRUE) {
  pathRaw <- normalizePath(pathRaw, winslash = "/", mustWork = TRUE)
  pathSmall <- normalizePath(pathSmall, winslash = "/", mustWork = TRUE)
  pathRoot <- normalizePath(pathRoot, winslash = "/", mustWork = TRUE)
  rawCheck <- .omip016ValidateRaw(pathRaw, pathSmall)
  sampleMap <- .omip016SampleMap(pathSmall)
  markerMap <- .omip016MarkerMap(pathSmall)
  cells <- sampleMap[sampleMap$role == "cells", , drop = FALSE]
  pathJo <- file.path(pathRaw, sampleMap$file[sampleMap$role == "workspace"])
  if (length(pathJo) != 1L) stop("Expected one workspace in sample_map.csv.")
  dir.create(dirname(pathOut), recursive = TRUE, showWarnings = FALSE)

  .acsCytofReplaceDir(pathOut, function(pathTmp) {
    paths <- .omip016Paths(pathTmp)
    for (p in paths[c("labels", "fcs", "reference", "qc")]) {
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
    }
    .omip016WriteCsv(rawCheck, file.path(paths$qc, "raw_validation.csv"))

    ws <- .omip016ParseWorkspace(pathJo, paths$workspace, cells$file, pathRoot)
    pathCd4 <- "Single cells/Lymphocytes/Lives/CD3+/CD4+"
    thresholds <- .omip016ManualThresholds(
      ws$pops[ws$pops$sample == cells$file[[1]], ],
      ws$vertices[ws$vertices$sample == cells$file[[1]], ],
      markerMap, pathCd4
    )
    .omip016WriteCsv(thresholds, file.path(paths$reference, "manual_thresholds.csv"))

    resp <- .omip016ResponseChannels(markerMap)
    countList <- list()
    sensList <- list()
    oneDimList <- list()
    freqList <- list()
    cd4Frames <- list()
    for (k in seq_len(nrow(cells))) {
      file <- cells$file[[k]]
      pops <- ws$pops[ws$pops$sample == file, , drop = FALSE]
      verts <- ws$vertices[ws$vertices$sample == file, , drop = FALSE]
      ffComp <- .omip016Compensate(.omip016ReadFcs(file.path(pathRaw, file)), ws$spill)
      member <- .omip016ApplyGates(ffComp, pops, verts)
      # Sensitivity of the reconstruction to the two FlowJo conventions that
      # cannot be checked against the workspace (it stores no statistics).
      alt <- list(
        boolean_precedence = .omip016ApplyGates(ffComp, pops, verts, booleanMode = "precedence"),
        no_edge_extension = .omip016ApplyGates(ffComp, pops, verts, extend = FALSE)
      )
      countList[[k]] <- .omip016PopulationCounts(member, file)
      sensList[[k]] <- do.call(rbind, lapply(names(alt), function(nm) {
        data.frame(
          file = file, variant = nm, path = colnames(member),
          count_primary = colSums(member), count_variant = colSums(alt[[nm]]),
          n_events_differing = colSums(member != alt[[nm]]),
          stringsAsFactors = FALSE
        )
      }))

      cd4 <- member[, pathCd4]
      ex <- flowCore::exprs(ffComp)
      labels <- data.frame(event = which(cd4))
      for (j in seq_len(nrow(thresholds))) {
        labels[[thresholds$marker[[j]]]] <- member[cd4, thresholds$gate_path[[j]]]
        above <- ex[cd4, thresholds$channel[[j]]] > thresholds$threshold_raw[[j]]
        oneDimList[[length(oneDimList) + 1L]] <- data.frame(
          file = file, marker = thresholds$marker[[j]],
          n_cd4 = sum(cd4),
          n_manual_pos = sum(labels[[thresholds$marker[[j]]]]),
          n_above_threshold = sum(above),
          n_disagree = sum(above != labels[[thresholds$marker[[j]]]]),
          stringsAsFactors = FALSE
        )
        freqList[[length(freqList) + 1L]] <- data.frame(
          file = file, SampleID = cells$SampleID[[k]], stim = cells$stim[[k]],
          pop = "CD4 T cells", cyt = thresholds$marker[[j]],
          chnl = thresholds$channel[[j]],
          n_cell = sum(cd4), count_man = sum(labels[[thresholds$marker[[j]]]]),
          stringsAsFactors = FALSE
        )
      }
      saveRDS(
        list(file = file, populations = member, primaryConventions = list(
          booleanMode = "left", extendVal = .omip016ExtendVal
        )),
        file.path(paths$labels, paste0(tools::file_path_sans_ext(file), "-populations.rds"))
      )
      saveRDS(labels, file.path(paths$labels, paste0(tools::file_path_sans_ext(file), "-cd4-manual.rds")))
      .omip016WriteCd4Fcs(
        ffComp, cd4, file.path(paths$fcs, file), file, ws$meta$compensation$name
      )
      cd4Frames[[file]] <- ffComp[cd4, ]
      rm(member, alt, ffComp, ex)
    }
    counts <- do.call(rbind, countList)
    .omip016WriteCsv(counts, file.path(paths$qc, "population_counts.csv"))
    .omip016WriteCsv(do.call(rbind, sensList), file.path(paths$qc, "convention_sensitivity.csv"))
    .omip016WriteCsv(do.call(rbind, oneDimList), file.path(paths$qc, "manual_gate_vs_threshold.csv"))

    freq <- do.call(rbind, freqList)
    freq$freq_stim_man <- 100 * freq$count_man / freq$n_cell
    uns <- freq[freq$stim == "uns", c("cyt", "freq_stim_man", "count_man", "n_cell")]
    names(uns) <- c("cyt", "freq_uns_man", "count_uns_man", "n_cell_uns")
    freq <- merge(freq, uns, by = "cyt", sort = FALSE)
    freq$freq_bs_man_raw <- freq$freq_stim_man - freq$freq_uns_man
    freq$freq_bs_man <- pmax(freq$freq_bs_man_raw, 0)
    freq <- freq[order(match(freq$stim, cells$stim), match(freq$cyt, thresholds$marker)), ]
    .omip016WriteCsv(freq, file.path(paths$reference, "manual_frequencies.csv"))

    spillBeads <- .omip016BeadSpillover(pathRaw, markerMap, sampleMap, colnames(ws$spill))
    spillBeads$spill_workspace <- ws$spill[cbind(
      spillBeads$fluorochrome_channel, spillBeads$detector
    )]
    .omip016WriteCsv(spillBeads, file.path(paths$qc, "compensation_beads_vs_workspace.csv"))

    # CD4 T cells, compensated then asinh-transformed, for gateStim().
    sampleMapGs <- data.frame(
      ind = as.character(seq_len(nrow(cells))), file = cells$file,
      SampleID = cells$SampleID, stim = cells$stim, stringsAsFactors = FALSE
    )
    batchList <- .omip016BatchList(sampleMapGs)
    cs <- flowWorkspace::flowSet_to_cytoset(flowCore::flowSet(cd4Frames))
    flowWorkspace::sampleNames(cs) <- cells$file
    gs <- flowWorkspace::GatingSet(cs)
    transChannels <- colnames(ws$spill)
    trans <- flowWorkspace::transformerList(
      transChannels,
      flowWorkspace::flow_trans(
        "asinh150",
        trans.fun = function(x) asinh(x / 150),
        inverse.fun = function(x) 150 * sinh(x)
      )
    )
    gs <- flowWorkspace::transform(gs, trans)
    dir.create(dirname(paths$gs), recursive = TRUE, showWarnings = FALSE)
    flowWorkspace::save_gs(gs, paths$gs)

    preprocessing <- list(
      version = .omip016Version,
      dataset = "OMIP-016 (FlowRepository FR-FCM-ZZ2T)",
      population = pathCd4,
      sampleMap = sampleMapGs,
      batchList = batchList,
      responseChannels = resp,
      settings = list(
        compensation = ws$meta$compensation$name,
        transform = "asinh(x / 150) on compensated fluorescence channels",
        transformChannels = transChannels,
        booleanMode = "left",
        extendVal = .omip016ExtendVal
      ),
      inputContentHash = unname(stats::setNames(rawCheck$md5, rawCheck$file)),
      workspaceMd5 = ws$meta$workspaceMd5,
      gitSha = if (exists(".git_sha", mode = "function")) .git_sha(pathRoot) else NA_character_,
      createdAt = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    saveRDS(preprocessing, .omip016PreprocessingFile(paths$gs))

    if (isTRUE(writePlots)) {
      .omip016PlotQc(cd4Frames, thresholds, cells, file.path(paths$qc, "cd4_response_channels.png"))
    }
    saveRDS(preprocessing, paths$manifest)
  })
  invisible(.omip016Paths(pathOut))
}

# Histograms of each response channel in CD4 T cells per tube, with the
# manual threshold, on the GatingSet's asinh scale.
.omip016PlotQc <- function(cd4Frames, thresholds, cells, path) {
  df <- do.call(rbind, lapply(names(cd4Frames), function(file) {
    ex <- flowCore::exprs(cd4Frames[[file]])
    do.call(rbind, lapply(seq_len(nrow(thresholds)), function(j) {
      data.frame(
        stim = cells$stim[match(file, cells$file)],
        marker = thresholds$marker[[j]],
        x = .omip016Transform(ex[, thresholds$channel[[j]]])
      )
    }))
  }))
  thr <- data.frame(marker = thresholds$marker, x = thresholds$threshold_trans)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_histogram(bins = 120) +
    ggplot2::geom_vline(data = thr, ggplot2::aes(xintercept = .data$x), colour = "red") +
    ggplot2::facet_grid(stim ~ marker, scales = "free") +
    ggplot2::scale_y_sqrt() +
    ggplot2::labs(x = "asinh(x / 150), compensated", y = "CD4 T cells (sqrt scale)") +
    ggplot2::theme_bw()
  ggplot2::ggsave(path, p, width = 30, height = 14, units = "cm", dpi = 120)
  invisible(path)
}
