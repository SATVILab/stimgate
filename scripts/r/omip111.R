# OMIP-111 real-data comparison. Source after analysis-runtime.R,
# acs_cytof-helper.R, sim-compare-freq_bs.R and acs_cytof-methods.R.
.omip111Markers <- function(population) {
  switch(population,
    CD4 = c(
      IFNg = "AF488-A", IL2 = "PE-A", TNF = "PE-Cy7-A",
      IL4_5 = "BV421-A", IL17A = "R718-A"
    ),
    CD8 = c(IFNg = "AF488-A", IL2 = "PE-A", TNF = "PE-Cy7-A"),
    stop("Unknown OMIP-111 population: ", population)
  )
}

.omip111SampleMap <- function(files) {
  pattern <- "^[EF][0-9]{2} (C57|BALB)_M([1-5])_(BFA|Stim)_TS WLSM[.]fcs$"
  if (length(files) != 20L || any(!grepl(pattern, files))) {
    stop("OMIP-111 requires the 20 explicitly identified full-panel mouse tubes.")
  }
  strain <- sub(pattern, "\\1", files)
  mouse <- paste0(strain, "_M", sub(pattern, "\\2", files))
  condition <- ifelse(sub(pattern, "\\3", files) == "BFA", "uns", "stim")
  out <- data.frame(
    sample = files, strain = strain, mouse = mouse,
    condition = condition, stringsAsFactors = FALSE
  )
  pairs <- split(out, out$mouse)
  if (length(pairs) != 10L || any(vapply(pairs, function(x) {
    nrow(x) != 2L || !setequal(x$condition, c("uns", "stim"))
  }, logical(1)))) {
    stop("Every mouse must have exactly one unstimulated and one stimulated tube.")
  }
  out[order(out$strain, out$mouse, match(out$condition, c("uns", "stim"))), ]
}

.omip111Preprocess <- function(pathRaw, pathPrepared, preprocessing, pathPython, pathImporter) {
  if (!is.finite(preprocessing$cofactor) || preprocessing$cofactor <= 0 ||
    !is.finite(preprocessing$parentRelativeTolerance) || preprocessing$parentRelativeTolerance < 0 ||
    !is.finite(preprocessing$frequencyTolerancePp) || preprocessing$frequencyTolerancePp < 0) {
    stop("Invalid preprocessing cofactor or import tolerance.")
  }
  files <- list.files(pathRaw, pattern = "^[EF][0-9]{2} .* WLSM[.]fcs$")
  sampleMap <- .omip111SampleMap(files)
  sources <- file.path(pathRaw, c(
    files, "OMIP-ICS-C57BL_6.wsp",
    "OMIP-ICS-BALB_c.wsp", "Supplementary_Note_1.xlsx"
  ))
  if (any(!file.exists(sources))) stop("Missing OMIP-111 source files.")
  fingerprint <- stats::setNames(vapply(sources, function(path) {
    digest::digest(file = path, algo = "sha256")
  }, character(1)), basename(sources))
  .acsCytofReplaceDir(pathPrepared, function(pathTmp) {
    masksDir <- file.path(pathTmp, "membership")
    status <- system2(pathPython, c(shQuote(pathImporter), shQuote(pathRaw), shQuote(masksDir)),
      stdout = file.path(pathTmp, "import.log"), stderr = file.path(pathTmp, "import-errors.log")
    )
    if (status != 0L) {
      stop(
        "FlowKit import failed: ",
        paste(readLines(file.path(pathTmp, "import-errors.log"), warn = FALSE), collapse = "\n")
      )
    }
    populations <- utils::read.csv(file.path(masksDir, "populations.csv"))
    validation <- utils::read.csv(file.path(masksDir, "gate-validation.csv"))
    references <- utils::read.csv(file.path(masksDir, "references.csv"))
    importer <- jsonlite::read_json(file.path(masksDir, "importer.json"))
    if (!identical(importer$FlowKit, "1.3.2")) stop("This analysis requires FlowKit 1.3.2.")
    if (nrow(populations) != 40L || nrow(references) != 160L ||
      nrow(validation) != 540L) {
      stop("Incomplete imported population/reference cohort.")
    }
    if (any(abs(references$importedLineage / references$flowjoLineage - 1) >
      preprocessing$parentRelativeTolerance) ||
      any(abs(references$importedNonNaive / references$flowjoNonNaive - 1) >
        preprocessing$parentRelativeTolerance)) {
      stop("Parent count import exceeds recorded tolerance.")
    }
    if (anyDuplicated(unique(references[c("population", "marker", "cd44LowerRaw", "cd44UpperRaw")])[
      c("population", "marker")])) {
      stop("A cytokine gate's CD44 interval differs between tubes.")
    }
    # Import agreement is judged on the authors' own frequency scale (count of
    # the non-naive gate, as printed in Figures S11/S12).
    references$frequencyDifferencePp <- 100 * (references$importedCount / references$importedNonNaive -
      references$flowjoCount / references$flowjoNonNaive)
    references$referenceAccepted <- abs(references$frequencyDifferencePp) <=
      preprocessing$frequencyTolerancePp
    input <- list()
    for (sid in sampleMap$sample) {
      frame <- flowCore::read.FCS(file.path(pathRaw, sid),
        transformation = FALSE,
        truncate_max_range = FALSE
      )
      rawExpression <- flowCore::exprs(frame)
      sampleInput <- list()
      for (population in c("CD4", "CD8")) {
        markers <- .omip111Markers(population)
        row <- populations[populations$sample == sid & populations$population == population, ]
        if (nrow(row) != 1L) stop("Missing/duplicate parent population.")
        membership <- utils::read.csv(file.path(masksDir, row$file))
        indices <- membership$eventIndex
        if (length(indices) != row$nCell || anyDuplicated(indices) ||
          any(indices < 1L | indices > nrow(rawExpression)) || length(indices) < 100L) {
          stop("Invalid or too-small parent membership: ", sid, " ", population)
        }
        raw <- rawExpression[indices, unname(markers), drop = FALSE]
        expression <- asinh(raw / preprocessing$cofactor)
        colnames(expression) <- names(markers)
        expression[] <- .omip111Float32(as.numeric(expression))
        # Rows are the authors' non-naive cells; `parent` marks each
        # cytokine's CD44+ parent (the non-naive cells within that tube's
        # author gate's CD44 interval, as the authors gated). The reference is
        # the authors' 2-D cytokine-by-CD44 gate membership.
        readMask <- function(prefix) {
          mask <- as.matrix(membership[paste0(prefix, names(markers))]) == 1L
          colnames(mask) <- names(markers)
          mask
        }
        parent <- readMask("parent_")
        manual <- readMask("")
        projected <- readMask("projected_")
        if (!all(is.finite(expression))) stop("Non-finite processed expression.")
        refs <- references[references$sample == sid & references$population == population, ]
        refs <- refs[match(names(markers), refs$marker), ]
        if (anyNA(refs$marker) || any(colSums(parent) != refs$importedParent) ||
          any(colSums(manual) != refs$importedCount) || any(manual & !parent) ||
          any(colSums(projected) != refs$projectedCount) || any(colSums(parent) < 100L) ||
          any(sum(membership$nonNaive) != refs$importedNonNaive)) {
          stop("Manual membership/reference counts disagree.")
        }
        # Check the exported masks against the raw FCS values read by R: the
        # cytokine-only cut exactly, and the CD44 interval and the authors'
        # rectangle on raw units, which differ from FlowJo's display-scale
        # gates only at their edges.
        lower <- matrix(refs$cytokineLowerRaw, nrow(raw), ncol(raw), byrow = TRUE)
        upper <- matrix(refs$cytokineUpperRaw, nrow(raw), ncol(raw), byrow = TRUE)
        cd44 <- rawExpression[indices, "[RealBlue 744]-A"]
        rawParent <- (membership$nonNaive == 1L) &
          outer(cd44, refs$cd44LowerRaw, ">=") & outer(cd44, refs$cd44UpperRaw, "<")
        rectangle <- (membership$nonNaive == 1L) & raw >= lower & raw < upper & rawParent
        if (!identical(unname(parent & raw > lower), unname(projected)) ||
          any(colSums(rawParent != parent) != refs$rawParentMismatch) ||
          any(colSums(rectangle != manual) != refs$rawRectangleMismatch)) {
          stop("Python reference masks do not match the raw FCS values read by R.")
        }
        cytokineCut <- raw > lower
        colnames(cytokineCut) <- names(markers)
        sampleInput[[population]] <- list(
          expression = expression, parent = parent, manual = manual,
          projected = projected, cytokineCut = cytokineCut, eventIndex = indices,
          cd44Raw = cd44, nonNaive = membership$nonNaive == 1L
        )
      }
      input[[sid]] <- sampleInput
    }
    manifest <- list(
      semantics = "omip111-v10", sourceSha256 = fingerprint,
      importer = importer, importerSourceSha256 = digest::digest(file = pathImporter, algo = "sha256"),
      preprocessing = preprocessing,
      scale = paste0("asinh(raw unmixed fluorescence / ", preprocessing$cofactor, "); float32"),
      compensation = "none; source FCS already spectrally unmixed",
      populations = paste(
        "reconstructed author-gated non-naive conventional CD4/CD8 T cells within each",
        "cytokine gate's CD44 interval (FlowJo display scale), one parent per cytokine"
      ),
      cd44Intervals = unique(references[c("population", "marker", "cd44LowerRaw", "cd44UpperRaw")])
    )
    saveRDS(
      list(
        samples = sampleMap, input = input, manifest = manifest,
        validation = validation, references = references, populations = populations
      ),
      file.path(pathTmp, "inputs.rds")
    )
    writeLines("Audited parent counts and labelled reference-frequency discrepancies", file.path(pathTmp, "COMPLETE"))
  })
  invisible(pathPrepared)
}

# One cytokine's parent cells: its author gate's CD44 interval.
.omip111MarkerInput <- function(input, marker) {
  keep <- input$parent[, marker]
  list(
    expression = input$expression[keep, marker, drop = FALSE],
    manual = input$manual[keep, marker], projected = input$projected[keep, marker]
  )
}

.omip111Float32 <- function(x) {
  readBin(writeBin(as.numeric(x), raw(), size = 4L), "numeric", n = length(x), size = 4L)
}

.omip111Counts <- function(countStim, nStim, countUns, nUns, manualStim, manualUns) {
  if (nStim <= 0 || nUns <= 0) stop("Frequency denominators must be positive.")
  data.frame(
    countStim = countStim, nCellStim = nStim, countUns = countUns,
    nCellUns = nUns, manualCountStim = manualStim, manualCountUns = manualUns,
    propStim = countStim / nStim, propUns = countUns / nUns,
    propBs = countStim / nStim - countUns / nUns,
    manualStim = manualStim / nStim, manualUns = manualUns / nUns,
    manualBs = manualStim / nStim - manualUns / nUns
  )
}

.omip111Run <- function(prepared, pathResults, settings, pathFbeta) {
  if (!identical(prepared$manifest$preprocessing, settings$preprocessing)) {
    stop("Prepared inputs use different preprocessing settings.")
  }
  # Missing competitor infrastructure is fatal, before running StimGate.
  if (!requireNamespace("cytoUtils", quietly = TRUE) ||
    !requireNamespace("reticulate", quietly = TRUE) || !file.exists(pathFbeta) ||
    !reticulate::py_module_available("numpy")) {
    stop("OMIP-111 requires cytoUtils, reticulate, numpy and the existing fbeta.py.")
  }
  fbetaEnv <- .simCompareFbetaEnvironment(pathFbeta = pathFbeta)
  .acsCytofReplaceDir(pathResults, function(pathTmp) {
    resultRows <- list()
    for (strain in c("C57", "BALB")) {
      sampleMap <- prepared$samples[prepared$samples$strain == strain, ]
      sampleMap$ind <- seq_len(nrow(sampleMap))
      batchList <- lapply(split(sampleMap, sampleMap$mouse), function(pair) {
        pair$ind[match(c("uns", "stim"), pair$condition)]
      })
      for (population in c("CD4", "CD8")) {
        # Each cytokine has its own parent (its author gate's CD44 interval),
        # so StimGate runs once per cytokine; cytokine-positive refinement is
        # off, so no joint gating across cytokines is lost.
        for (marker in names(.omip111Markers(population))) {
          tubes <- stats::setNames(lapply(sampleMap$sample, function(sid) {
            .omip111MarkerInput(prepared$input[[sid]][[population]], marker)
          }), sampleMap$sample)
          data <- lapply(tubes, `[[`, "expression")
          pathStim <- file.path(pathTmp, strain, population, marker, "stimgate")
          message("Running StimGate: ", strain, " ", population, " ", marker)
          # calcCytPosGates is a gateStim() argument; the other settings are
          # stimControl() settings.
          control <- do.call(
            stimgate::stimControl,
            settings$control[setdiff(names(settings$control), "calcCytPosGates")]
          )
          .analysis_with_seed(settings$seed, stimgate::gateStim(
            pathProject = pathStim, .data = data, batchList = batchList,
            marker = marker, biasUns = NULL, calcCytPosGates = isTRUE(settings$control$calcCytPosGates),
            control = control
          ))
          gates <- stimgate::getStimGates(pathStim)
          gates <- gates[gates$gateName == "loc_minClust", ]
          stats <- stimgate::getStimStats(pathStim)
          stats <- stats[stats$gateName == "loc_minClust", ]
          if (nrow(gates) != 5L) stop("Incomplete StimGate gate table.")
          for (batch in batchList) {
            iUns <- batch[[1]]
            iStim <- batch[[2]]
            sidUns <- sampleMap$sample[[iUns]]
            sidStim <- sampleMap$sample[[iStim]]
            xUns <- data[[iUns]][, marker]
            xStim <- data[[iStim]][, marker]
            refUns <- tubes[[iUns]]$manual
            refStim <- tubes[[iStim]]$manual
            info <- data.frame(
              strain = strain, mouse = sampleMap$mouse[[iStim]],
              population = population, marker = marker, sampleStim = sidStim, sampleUns = sidUns,
              projectedBs = mean(tubes[[iStim]]$projected) - mean(tubes[[iUns]]$projected)
            )
            gate <- gates[as.integer(as.character(gates$ind)) == iStim & gates$marker == marker, ]
            if (nrow(gate) != 1L) stop("Missing/duplicate final StimGate gate.")
            marginal <- stats[as.integer(as.character(stats$ind)) == iStim &
              grepl(paste0(marker, "~+~"), stats$cytCombn, fixed = TRUE), ]
            if (nrow(marginal) != 1L || marginal$nCellStim != length(xStim) ||
              marginal$nCellUns != length(xUns)) {
              stop("Incomplete or inconsistent StimGate combination counts.")
            }
            counts <- .omip111Counts(
              marginal$countStim, length(xStim),
              marginal$countUns, length(xUns), sum(refStim), sum(refUns)
            )
            provenance <- gate[, intersect(c(
              "locGenerated", "locGeneratedDirect", "locSource",
              "locReason", "gate", "gateCyt", "locShareLimit", "locShareProposed"
            ), names(gate)), drop = FALSE]
            resultRows[[length(resultRows) + 1L]] <- dplyr::bind_cols(
              info,
              data.frame(
                method = "StimGate", threshold = gate$gate,
                thresholdOrigin = as.character(gate$locSource),
                thresholdFallbackUsed = !isTRUE(gate$locGeneratedDirect), error = NA_character_
              ),
              counts, provenance
            )
            for (method in c("F-beta", "Tailgate")) {
              comparator <- .acsCytofComparatorSettings(
                if (method == "F-beta") "fbeta_default" else "tailgate_simulation_allcells"
              )
              threshold <- .acsCytofThresholdOne(
                method = method, xUns = xUns, xStim = xStim,
                settings = comparator, pathFbeta = pathFbeta, fbetaEnv = fbetaEnv
              )
              counts <- .omip111Counts(
                if (is.finite(threshold$threshold)) sum(xStim > threshold$threshold) else NA_real_,
                length(xStim),
                if (is.finite(threshold$threshold)) sum(xUns > threshold$threshold) else NA_real_,
                length(xUns), sum(refStim), sum(refUns)
              )
              resultRows[[length(resultRows) + 1L]] <- dplyr::bind_cols(
                info,
                data.frame(method = method, error = if (threshold$thresholdOrigin == "runtime_error") {
                  "Comparator estimation failed; see render log"
                } else {
                  NA_character_
                }), threshold, counts
              )
            }
          }
          saveRDS(
            list(
              context = prepared$manifest, settings = settings,
              channelSettings = stimgate::stimgateMetaReadSettingsChnls(pathStim)
            ),
            file.path(dirname(pathStim), "manifest.rds")
          )
          saveRDS(dplyr::bind_rows(resultRows), file.path(pathTmp, "progress.rds"))
        }
      }
    }
    results <- dplyr::bind_rows(resultRows)
    writeLines(capture.output(utils::sessionInfo()), file.path(pathTmp, "session-info.txt"))
    # Paper-scale context: the FlowJo count of the non-naive gate, as printed
    # in Figures S11/S12 (stimulated tubes only, no background subtraction).
    originals <- prepared$references |>
      dplyr::mutate(flowjoProp = .data$flowjoCount / .data$flowjoNonNaive) |>
      dplyr::select(.data$sample, .data$population, .data$marker, .data$flowjoProp, .data$referenceAccepted)
    results <- results |>
      dplyr::left_join(dplyr::rename(originals, sampleStim = "sample", flowjoStim = "flowjoProp", referenceAcceptedStim = "referenceAccepted"),
        by = c("sampleStim", "population", "marker")
      ) |>
      dplyr::left_join(dplyr::rename(originals, sampleUns = "sample", flowjoUns = "flowjoProp", referenceAcceptedUns = "referenceAccepted"),
        by = c("sampleUns", "population", "marker")
      ) |>
      dplyr::mutate(
        flowjoBs = .data$flowjoStim - .data$flowjoUns,
        projectedVsAuthorDifferencePp = 100 * (.data$projectedBs - .data$manualBs)
      )
    results$errorPp <- 100 * (results$propBs - results$manualBs)
    # The reference is the reconstructed author 2-D gate; tubes whose import
    # misses FlowJo's frequency by more than the tolerance are not scored.
    results$referenceAccepted <- results$referenceAcceptedStim & results$referenceAcceptedUns
    results$unscreenedErrorPp <- results$errorPp
    results$errorPp[!results$referenceAccepted] <- NA_real_
    results$absoluteErrorPp <- abs(results$errorPp)
    key <- results[c("strain", "mouse", "population", "marker", "method")]
    if (nrow(results) != 240L || anyDuplicated(key)) stop("Incomplete comparison cohort.")
    resultsNonNaive <- .omip111RunNonNaive(prepared, pathTmp, settings, pathFbeta, fbetaEnv)
    saveRDS(
      list(
        results = results, resultsNonNaive = resultsNonNaive,
        settings = settings, preprocessing = prepared$manifest
      ),
      file.path(pathTmp, "result.rds")
    )
    utils::write.csv(results, file.path(pathTmp, "comparison.csv"), row.names = FALSE)
    utils::write.csv(resultsNonNaive, file.path(pathTmp, "comparison-non-naive-fit.csv"), row.names = FALSE)
    writeLines("Complete: 240 method/sample/population/marker outcomes.", file.path(pathTmp, "COMPLETE"))
  })
  invisible(pathResults)
}

# The same methods fitted on every non-naive cell (the authors' parent gate)
# instead of each cytokine's CD44 parent. All cytokines share these cells, so
# StimGate runs once per strain and population. Counts are over the non-naive
# cells of each tube; .omip111FitEval() rescores them on either population.
.omip111RunNonNaive <- function(prepared, pathTmp, settings, pathFbeta, fbetaEnv) {
  rows <- list()
  for (strain in c("C57", "BALB")) {
    sampleMap <- prepared$samples[prepared$samples$strain == strain, ]
    sampleMap$ind <- seq_len(nrow(sampleMap))
    batchList <- lapply(split(sampleMap, sampleMap$mouse), function(pair) {
      pair$ind[match(c("uns", "stim"), pair$condition)]
    })
    for (population in c("CD4", "CD8")) {
      markers <- names(.omip111Markers(population))
      data <- stats::setNames(lapply(sampleMap$sample, function(sid) {
        input <- prepared$input[[sid]][[population]]
        input$expression[input$nonNaive, markers, drop = FALSE]
      }), sampleMap$sample)
      pathStim <- file.path(pathTmp, strain, population, "non-naive", "stimgate")
      message("Running StimGate on non-naive cells: ", strain, " ", population)
      control <- do.call(
        stimgate::stimControl,
        settings$control[setdiff(names(settings$control), "calcCytPosGates")]
      )
      .analysis_with_seed(settings$seed, stimgate::gateStim(
        pathProject = pathStim, .data = data, batchList = batchList,
        marker = markers, biasUns = NULL, calcCytPosGates = isTRUE(settings$control$calcCytPosGates),
        control = control
      ))
      gates <- stimgate::getStimGates(pathStim)
      gates <- gates[gates$gateName == "loc_minClust", ]
      if (nrow(gates) != 5L * length(markers)) stop("Incomplete non-naive StimGate gate table.")
      for (batch in batchList) {
        iUns <- batch[[1]]
        iStim <- batch[[2]]
        for (marker in markers) {
          xUns <- data[[iUns]][, marker]
          xStim <- data[[iStim]][, marker]
          gate <- gates[as.integer(as.character(gates$ind)) == iStim & gates$marker == marker, ]
          if (nrow(gate) != 1L) stop("Missing/duplicate non-naive StimGate gate.")
          fits <- list(StimGate = data.frame(
            threshold = gate$gate, gate = gate$gate, gateCyt = gate$gateCyt,
            thresholdOrigin = as.character(gate$locSource),
            thresholdFallbackUsed = !isTRUE(gate$locGeneratedDirect), error = NA_character_
          ))
          for (method in c("F-beta", "Tailgate")) {
            comparator <- .acsCytofComparatorSettings(
              if (method == "F-beta") "fbeta_default" else "tailgate_simulation_allcells"
            )
            threshold <- .acsCytofThresholdOne(
              method = method, xUns = xUns, xStim = xStim,
              settings = comparator, pathFbeta = pathFbeta, fbetaEnv = fbetaEnv
            )
            threshold$error <- if (threshold$thresholdOrigin == "runtime_error") {
              "Comparator estimation failed; see render log"
            } else {
              NA_character_
            }
            fits[[method]] <- threshold
          }
          for (method in names(fits)) {
            g <- fits[[method]]$threshold
            rows[[length(rows) + 1L]] <- dplyr::bind_cols(
              data.frame(
                strain = strain, mouse = sampleMap$mouse[[iStim]], population = population,
                marker = marker, sampleStim = sampleMap$sample[[iStim]],
                sampleUns = sampleMap$sample[[iUns]], method = method
              ),
              fits[[method]],
              data.frame(
                countStim = if (is.finite(g)) sum(xStim > g) else NA_real_, nCellStim = length(xStim),
                countUns = if (is.finite(g)) sum(xUns > g) else NA_real_, nCellUns = length(xUns)
              )
            )
          }
        }
      }
      # With cytokine-positive refinement off, the package's counts are the
      # plain x > gate counts used above.
      stats <- stimgate::getStimStats(pathStim)
      stats <- stats[stats$gateName == "loc_minClust", ]
      sg <- dplyr::bind_rows(rows)
      sg <- sg[sg$method == "StimGate" & sg$strain == strain & sg$population == population, ]
      for (i in seq_len(nrow(sg))) {
        iStim <- match(sg$sampleStim[[i]], sampleMap$sample)
        marginal <- stats[as.integer(as.character(stats$ind)) == iStim &
          grepl(paste0(sg$marker[[i]], "~+~"), stats$cytCombn, fixed = TRUE), ]
        if (sum(marginal$countStim) != sg$countStim[[i]] || sum(marginal$countUns) != sg$countUns[[i]]) {
          stop("Non-naive StimGate counts differ from the package statistics.")
        }
      }
    }
  }
  out <- dplyr::bind_rows(rows)
  if (nrow(out) != 240L || anyDuplicated(out[c("strain", "mouse", "population", "marker", "method")])) {
    stop("Incomplete non-naive comparison cohort.")
  }
  out
}

# Rescore both fits on either evaluation population, against either manual
# reference. The fit ("cd44": each cytokine's CD44 parent; "nonnaive": every
# non-naive cell) only sets the threshold. Positives are always the CD44+
# parent cells above it (the authors' CD44 restriction, never inferred from
# the cytokine gate); evaluate sets the denominator: the CD44+ parent ("cd44")
# or every non-naive cell ("nonnaive"). reference: "bivariate" (the authors'
# CD44-by-cytokine gate) or "univariate" (their cytokine cutoff within the
# CD44+ parent). Counts on each fit's own population must reproduce the
# saved counts.
.omip111FitEval <- function(prepared, results, resultsNonNaive, fit, evaluate, reference) {
  fitted <- if (fit == "cd44") results else resultsNonNaive
  fitted <- fitted[c("strain", "mouse", "population", "marker", "method",
    "sampleStim", "sampleUns", "threshold", "countStim", "countUns")]
  accepted <- results[c("strain", "mouse", "population", "marker", "method", "referenceAccepted")]
  tube <- function(sid, population, marker, threshold) {
    input <- prepared$input[[sid]][[population]]
    parent <- input$parent[, marker]
    x <- input$expression[, marker]
    keep <- if (evaluate == "cd44") parent else input$nonNaive
    pos <- keep & parent & x > threshold
    own <- if (fit == "cd44") parent & x > threshold else input$nonNaive & x > threshold
    manual <- if (reference == "bivariate") input$manual[, marker] else parent & input$cytokineCut[, marker]
    list(prop = sum(pos) / sum(keep), own = sum(own), manual = sum(manual & keep) / sum(keep))
  }
  rows <- lapply(seq_len(nrow(fitted)), function(i) {
    r <- fitted[i, ]
    s <- tube(r$sampleStim, r$population, r$marker, if (is.finite(r$threshold)) r$threshold else Inf)
    u <- tube(r$sampleUns, r$population, r$marker, if (is.finite(r$threshold)) r$threshold else Inf)
    if (is.finite(r$threshold) && (s$own != r$countStim || u$own != r$countUns)) {
      stop("Recounted gates differ from saved counts: ", r$sampleStim, " ", r$marker, " ", r$method)
    }
    data.frame(r[c("strain", "mouse", "population", "marker", "method", "sampleStim", "sampleUns", "threshold")],
      propStim = s$prop, propUns = u$prop,
      propBs = if (is.finite(r$threshold)) s$prop - u$prop else NA_real_,
      manualStim = s$manual, manualUns = u$manual, manualBs = s$manual - u$manual
    )
  })
  out <- dplyr::left_join(dplyr::bind_rows(rows), accepted,
    by = c("strain", "mouse", "population", "marker", "method"))
  # The cytokine cutoff is applied to raw values exactly; only the 2-D gate
  # reconstruction carries the import tolerance.
  if (reference == "univariate") out$referenceAccepted <- TRUE
  out$fit <- fit
  out$evaluate <- evaluate
  out$reference <- reference
  out
}

.omip111Summary <- function(results) {
  results |>
    dplyr::group_by(.data$strain, .data$population, .data$marker, .data$method) |>
    dplyr::summarise(
      nMice = dplyr::n(), nFinite = sum(is.finite(.data$propBs)),
      nAgreementFinite = sum(is.finite(.data$errorPp)),
      nErrors = sum(!is.na(.data$error)),
      nFallback = sum(.data$thresholdFallbackUsed, na.rm = TRUE),
      meanErrorPp = if (any(is.finite(.data$errorPp))) mean(.data$errorPp, na.rm = TRUE) else NA_real_,
      meanAbsoluteErrorPp = if (any(is.finite(.data$absoluteErrorPp))) {
        mean(.data$absoluteErrorPp, na.rm = TRUE)
      } else {
        NA_real_
      }, .groups = "drop"
    )
}

# Reconstruct calls from saved gates and prepared float32 expression and score
# them against the author 2-D gate masks. Source omip016-methods.R first.
.omip111Classification <- function(prepared, results) {
  required <- c("strain", "mouse", "population", "marker", "method",
                "sampleStim", "threshold", "gate", "gateCyt", "countStim", "nCellStim")
  if (!all(required %in% names(results)) || is.null(prepared$input)) {
    stop("Saved OMIP-111 inputs/results lack essential classification fields.")
  }
  key <- c("strain", "mouse", "population", "marker", "method")
  if (anyDuplicated(results[key])) stop("Duplicate OMIP-111 classification outcomes.")
  groups <- split(seq_len(nrow(results)), interaction(
    results$strain, results$population, results$mouse, results$method, drop = TRUE
  ))
  rows <- lapply(groups, function(indices) {
    g <- results[indices, , drop = FALSE]
    markers <- names(.omip111Markers(g$population[[1]]))
    if (!setequal(g$marker, markers) || length(unique(g$sampleStim)) != 1L) {
      stop("Incomplete marker gates for an OMIP-111 stimulated tube.")
    }
    input <- prepared$input[[g$sampleStim[[1]]]][[g$population[[1]]]]
    if (is.null(input$expression) || is.null(input$manual) || is.null(input$parent) ||
        !identical(dim(input$expression), dim(input$manual)) ||
        !identical(dim(input$parent), dim(input$manual)) ||
        !all(markers %in% colnames(input$manual)) || !all(markers %in% colnames(input$parent)) ||
        anyNA(input$manual) || !is.logical(input$manual) || any(!is.finite(input$expression))) {
      stop("Missing or inconsistent prepared OMIP-111 expression/manual masks.")
    }
    # Each cytokine is scored on its own parent (its author gate's CD44
    # interval), as it was gated.
    dplyr::bind_rows(lapply(seq_len(nrow(g)), function(i) {
      row <- g[i, , drop = FALSE]
      cells <- .omip111MarkerInput(input, row$marker)
      x <- cells$expression
      manual <- cells$manual
      if (is.na(row$nCellStim) || row$nCellStim != nrow(x)) {
        stop("Missing or inconsistent prepared OMIP-111 expression/manual masks.")
      }
      if (startsWith(row$method, "StimGate")) {
        if (!is.finite(row$gate) || is.na(row$gateCyt)) {
          stop("Essential saved StimGate base/conditional gates are missing.")
        }
        pos <- .omip016Classify(x, row$gate, row$gateCyt)[, 1]
        if (is.na(row$countStim) || sum(pos) != row$countStim) {
          stop("StimGate reconstructed positive counts differ from saved marginal counts.")
        }
      } else {
        pos <- .omip016Classify(x, row$threshold, NA_real_)[, 1]
      }
      counts <- data.frame(
        tp = sum(pos & manual), fp = sum(pos & !manual),
        fn = sum(!pos & manual), tn = sum(!pos & !manual)
      )
      if (!startsWith(row$method, "StimGate") && !is.finite(row$threshold)) {
        counts[] <- NA_real_
      }
      nDefined <- if (is.finite(rowSums(counts))) nrow(x) else 0L
      dplyr::bind_cols(row[key], data.frame(n = nrow(x), n_defined = nDefined), counts,
                       .omipClassificationMetrics(counts))
    }))
  })
  dplyr::bind_rows(rows)
}

.omip111ClassificationSummary <- function(classification) {
  counts <- c("tp", "fp", "fn", "tn")
  metrics <- c("fdp", "precision", "sensitivity", "specificity", "f1")
  summary <- classification |>
    dplyr::group_by(.data$strain, .data$population, .data$marker, .data$method) |>
    dplyr::summarise(
      n_mice = dplyr::n(),
      n_defined = sum(is.finite(.data$tp + .data$fp + .data$fn + .data$tn)),
      n_cells_defined = sum(.data$n[is.finite(.data$tp + .data$fp + .data$fn + .data$tn)]),
      n = sum(.data$n),
      dplyr::across(dplyr::all_of(counts),
        ~ if (any(is.finite(.x))) sum(.x[is.finite(.x)]) else NA_real_),
      dplyr::across(dplyr::all_of(metrics), list(
        median = ~ if (any(is.finite(.x))) stats::median(.x[is.finite(.x)]) else NA_real_,
        n_defined = ~ sum(is.finite(.x))
      )), .groups = "drop"
    )
  # Pooled ratios use summed counts; medians give each defined mouse equal weight.
  dplyr::bind_cols(summary, .omipClassificationMetrics(summary))
}

# Supplementary Figures S11 (C57BL/6) and S12 (BALB/c): every non-naive
# stimulated CD4 T cell on CD44 against each cytokine, with the author
# rectangles and their frequencies (count / non-naive count, as printed).
# Axes use the analysis scale, asinh(raw / cofactor); points are coloured by
# the cell count in their bin, as in FlowJo's pseudocolour plots. The same
# panels carry each method's final gate (.omip111MethodFigureData()).
.omip111AuthorFigureRows <- c(
  IFNg = "IFNγ AF488", TNF = "TNF PE-Cy7", IL2 = "IL-2 PE",
  IL4_5 = "IL-4/5 BV421", IL17A = "IL-17 R718"
)

.omip111AuthorFigureLabel <- function(x) {
  vapply(x, function(v) {
    formatC(v, digits = if (v < 1) 2L else 3L, format = "fg", flag = "#")
  }, character(1))
}

# Non-naive cells of each stimulated tube with pseudocolour densities, and
# the author gate of every panel (rectangle, counts and denominators).
.omip111FigurePanels <- function(prepared, strain, cofactor, population = "CD4", bins = 128L) {
  markers <- intersect(names(.omip111AuthorFigureRows), names(.omip111Markers(population)))
  stim <- prepared$samples[prepared$samples$strain == strain &
    prepared$samples$condition == "stim", ]
  if (nrow(stim) != 5L) stop("Expected five stimulated tubes for ", strain, ".")
  stim$column <- paste("Mouse", sub(".*_M", "", stim$mouse))
  refs <- prepared$references[prepared$references$population == population, ]
  cells <- list()
  gates <- list()
  for (i in seq_len(nrow(stim))) {
    input <- prepared$input[[stim$sample[[i]]]][[population]]
    if (is.null(input$cd44Raw) || is.null(input$nonNaive)) {
      stop("Prepared OMIP-111 inputs lack CD44 or non-naive membership; rerun preprocessing.")
    }
    keep <- input$nonNaive
    for (marker in markers) {
      ref <- refs[refs$sample == stim$sample[[i]] & refs$marker == marker, ]
      if (nrow(ref) != 1L) stop("Missing author gate for ", stim$sample[[i]], " ", marker)
      x <- asinh(input$cd44Raw[keep] / cofactor)
      y <- input$expression[keep, marker]
      # Bin counts over the panel's range give FlowJo-style pseudocolour.
      bx <- cut(x, bins, labels = FALSE)
      by <- cut(y, bins, labels = FALSE)
      key <- (bx - 1L) * bins + by
      density <- tabulate(key, bins * bins)[key]
      ord <- order(density)
      cells[[length(cells) + 1L]] <- data.frame(
        column = stim$column[[i]], marker = marker,
        x = x[ord], y = y[ord], density = density[ord]
      )
      gates[[length(gates) + 1L]] <- data.frame(
        sample = stim$sample[[i]], column = stim$column[[i]], marker = marker,
        xmin = asinh(ref$cd44LowerRaw / cofactor), xmax = asinh(ref$cd44UpperRaw / cofactor),
        ymin = asinh(ref$cytokineLowerRaw / cofactor), ymax = asinh(ref$cytokineUpperRaw / cofactor),
        frequency = 100 * ref$importedCount / ref$importedNonNaive,
        flowjoFrequency = 100 * ref$flowjoCount / ref$flowjoNonNaive,
        parentFrequency = 100 * ref$importedCount / ref$importedParent,
        nonNaive = sum(keep), parent = ref$importedParent
      )
    }
  }
  rows <- function(x) {
    x$row <- factor(.omip111AuthorFigureRows[x$marker], levels = .omip111AuthorFigureRows[markers])
    x
  }
  list(cells = rows(dplyr::bind_rows(cells)), gates = rows(dplyr::bind_rows(gates)))
}

.omip111AuthorFigureData <- function(prepared, strain, cofactor, bins = 128L) {
  figureData <- .omip111FigurePanels(prepared, strain, cofactor, "CD4", bins)
  figureData$gates$label <- .omip111AuthorFigureLabel(figureData$gates$frequency)
  figureData
}

# Each method's stimulated-tube gate and how it counts, for one fit population
# ("cd44": each cytokine's CD44 parent; "nonnaive": every non-naive cell).
# Counts on the fit's own population must reproduce the saved counts.
.omip111FigureGate <- function(prepared, results, resultsNonNaive, sid, population,
                               marker, method, fit) {
  fitted <- if (fit == "cd44") results else resultsNonNaive
  row <- fitted[fitted$sampleStim == sid & fitted$population == population &
    fitted$marker == marker & fitted$method == method, , drop = FALSE]
  if (nrow(row) != 1L) stop("Missing saved ", method, " result for ", sid, " ", marker)
  input <- prepared$input[[sid]][[population]]
  parent <- input$parent[, marker]
  threshold <- row$threshold
  above <- is.finite(threshold) & input$expression[, marker] > threshold
  own <- above & (if (fit == "cd44") parent else input$nonNaive)
  if (is.finite(threshold) && sum(own) != row$countStim) {
    stop(method, " gate does not reproduce the saved stimulated count: ", sid, " ", marker)
  }
  # Positives are the CD44+ parent cells above the gate, whatever the fit.
  above <- above & parent
  list(
    threshold = threshold,
    nonNaiveFrequency = 100 * sum(above & input$nonNaive) / sum(input$nonNaive),
    parentFrequency = 100 * sum(above & parent) / sum(parent)
  )
}

# Bivariate panels (every non-naive cell, CD44 against the cytokine) with each
# method's threshold and the authors' manual gate, labelled with stimulated-
# tube frequencies over every non-naive cell. Every gate is drawn across the
# author rectangle's CD44 interval only.
.omip111MethodFigureData <- function(prepared, results, resultsNonNaive, strain, methods,
                                     cofactor, fit, population = "CD4", bins = 128L) {
  figureData <- .omip111FigurePanels(prepared, strain, cofactor, population, bins)
  author <- figureData$gates
  gates <- dplyr::bind_rows(lapply(seq_len(nrow(author)), function(i) {
    a <- author[i, , drop = FALSE]
    dplyr::bind_rows(lapply(methods, function(method) {
      g <- .omip111FigureGate(prepared, results, resultsNonNaive, a$sample, population,
        a$marker, method, fit)
      data.frame(a[c("column", "marker", "row", "xmin", "xmax", "ymax")],
        method = method, ymin = g$threshold, frequency = g$nonNaiveFrequency
      )
    }))
  }))
  gates$method <- factor(gates$method, levels = methods)
  gates$label <- ifelse(is.finite(gates$ymin),
    paste0(gates$method, ": ", .omip111AuthorFigureLabel(gates$frequency)),
    paste0(gates$method, ": no gate")
  )
  author$label <- paste0("Manual: ", .omip111AuthorFigureLabel(author$frequency))
  figureData$reference <- author
  figureData$gates <- gates
  figureData
}

.omip111FigureGateScales <- function(methods, methodColours) {
  if (!all(methods %in% names(methodColours))) stop("Missing method colours.")
  list(
    ggplot2::scale_colour_manual(values = c(methodColours[methods], Manual = "black"), name = NULL),
    ggplot2::scale_linetype_manual(
      values = c(stats::setNames(rep("solid", length(methods)), methods), Manual = "22"), name = NULL
    )
  )
}

# Per-panel label blocks: the manual gate first, then each method, each in
# its gate's colour.
.omip111FigureLabels <- function(author, gates, methods, methodColours, keys) {
  dplyr::bind_rows(
    data.frame(author[c(keys, "label")], line = 0L, colour = "black"),
    data.frame(gates[c(keys, "label")],
      line = match(as.character(gates$method), methods),
      colour = unname(methodColours[as.character(gates$method)])
    )
  )
}

# All methods' gates and the authors' manual gate on one bivariate figure,
# over grey density shading. The author rectangle is drawn as it is (sides and
# top in grey, its cytokine cutoff dashed black); each method's gate is its
# threshold line across the rectangle's CD44 interval, coloured by method.
.omip111MethodsFigure <- function(figureData, cofactor, methodColours) {
  gates <- figureData$gates
  author <- figureData$reference
  methods <- levels(gates$method)
  cells <- figureData$cells
  # Light-to-dark grey by log bin count (the colour scale is kept for gates).
  shade <- scales::rescale(log10(cells$density))
  cells$grey <- grDevices::grey(0.85 - 0.6 * shade)
  thresholds <- dplyr::bind_rows(
    data.frame(gates[is.finite(gates$ymin), c("column", "row", "xmin", "xmax", "ymin")],
      gate = as.character(gates$method[is.finite(gates$ymin)])),
    data.frame(author[c("column", "row", "xmin", "xmax", "ymin")], gate = "Manual")
  )
  thresholds$gate <- factor(thresholds$gate, levels = c(methods, "Manual"))
  # Draw the manual line last so it stays visible where a method coincides.
  thresholds <- thresholds[order(thresholds$gate == "Manual"), ]
  labels <- .omip111FigureLabels(author, gates, methods, methodColours, c("column", "marker", "row"))
  ggplot2::ggplot(cells, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point(colour = cells$grey, size = 0.12, shape = 16) +
    ggplot2::geom_segment(
      data = dplyr::bind_rows(
        dplyr::transmute(author, .data$column, .data$row, x = .data$xmin, xend = .data$xmin,
          y = .data$ymin, yend = .data$ymax),
        dplyr::transmute(author, .data$column, .data$row, x = .data$xmax, xend = .data$xmax,
          y = .data$ymin, yend = .data$ymax),
        dplyr::transmute(author, .data$column, .data$row, x = .data$xmin, xend = .data$xmax,
          y = .data$ymax, yend = .data$ymax)
      ),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      inherit.aes = FALSE, colour = "grey30", linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      data = thresholds,
      ggplot2::aes(x = .data$xmin, xend = .data$xmax, y = .data$ymin, yend = .data$ymin,
        colour = .data$gate, linetype = .data$gate),
      inherit.aes = FALSE, linewidth = 0.7
    ) +
    .omip111FigureGateScales(methods, methodColours) +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = -Inf, y = Inf, label = .data$label, vjust = 1.3 + 1.25 * .data$line),
      colour = labels$colour, inherit.aes = FALSE, hjust = -0.04, size = 2.5
    ) +
    .omip111FigureScales(cofactor) +
    ggplot2::theme(legend.position = "bottom")
}

# Univariate panels: each cytokine's CD44 parent in the stimulated tube, with
# each method's threshold and the authors' cytokine cutoff, labelled with
# stimulated-tube frequencies over that CD44 parent. The manual frequency is
# the authors' 2-D gate count over the same parent.
.omip111UnivariateFigureData <- function(prepared, results, resultsNonNaive, strain, methods,
                                         cofactor, fit, population = "CD4") {
  markers <- intersect(names(.omip111AuthorFigureRows), names(.omip111Markers(population)))
  stim <- prepared$samples[prepared$samples$strain == strain &
    prepared$samples$condition == "stim", ]
  stim$column <- paste("Mouse", sub(".*_M", "", stim$mouse))
  refs <- prepared$references[prepared$references$population == population, ]
  cells <- list()
  author <- list()
  gates <- list()
  for (i in seq_len(nrow(stim))) {
    input <- prepared$input[[stim$sample[[i]]]][[population]]
    for (marker in markers) {
      ref <- refs[refs$sample == stim$sample[[i]] & refs$marker == marker, ]
      parent <- input$parent[, marker]
      cells[[length(cells) + 1L]] <- data.frame(
        column = stim$column[[i]], marker = marker, x = input$expression[parent, marker]
      )
      author[[length(author) + 1L]] <- data.frame(
        column = stim$column[[i]], marker = marker, x = asinh(ref$cytokineLowerRaw / cofactor),
        frequency = 100 * ref$importedCount / ref$importedParent
      )
      for (method in methods) {
        g <- .omip111FigureGate(prepared, results, resultsNonNaive, stim$sample[[i]], population,
          marker, method, fit)
        gates[[length(gates) + 1L]] <- data.frame(
          column = stim$column[[i]], marker = marker, method = method,
          x = g$threshold, frequency = g$parentFrequency
        )
      }
    }
  }
  panel <- function(x) {
    x$row <- factor(.omip111AuthorFigureRows[x$marker], levels = .omip111AuthorFigureRows[markers])
    x$panel <- factor(paste0(x$row, ", ", x$column),
      levels = as.vector(t(outer(.omip111AuthorFigureRows[markers], sort(unique(stim$column)), paste, sep = ", "))))
    x
  }
  gates <- panel(dplyr::bind_rows(gates))
  gates$method <- factor(gates$method, levels = methods)
  gates$label <- ifelse(is.finite(gates$x),
    paste0(gates$method, ": ", .omip111AuthorFigureLabel(gates$frequency)),
    paste0(gates$method, ": no gate")
  )
  author <- panel(dplyr::bind_rows(author))
  author$label <- paste0("Manual: ", .omip111AuthorFigureLabel(author$frequency))
  list(cells = panel(dplyr::bind_rows(cells)), gates = gates, reference = author)
}

.omip111UnivariateFigure <- function(figureData, cofactor, methodColours, bins = 80L) {
  gates <- figureData$gates
  author <- figureData$reference
  methods <- levels(gates$method)
  raw <- c(-1e3, 0, 1e3, 1e4, 1e5, 1e6)
  thresholds <- dplyr::bind_rows(
    data.frame(gates[is.finite(gates$x), c("panel", "x")], gate = as.character(gates$method[is.finite(gates$x)])),
    data.frame(author[c("panel", "x")], gate = "Manual")
  )
  thresholds$gate <- factor(thresholds$gate, levels = c(methods, "Manual"))
  thresholds <- thresholds[order(thresholds$gate == "Manual"), ]
  labels <- .omip111FigureLabels(author, gates, methods, methodColours, "panel")
  ggplot2::ggplot(figureData$cells, ggplot2::aes(.data$x)) +
    ggplot2::geom_histogram(bins = bins, fill = "grey60", colour = NA) +
    # Keep every gate, however high, inside its panel.
    ggplot2::geom_blank(data = thresholds, ggplot2::aes(x = .data$x + 0.3)) +
    ggplot2::geom_vline(
      data = thresholds, ggplot2::aes(xintercept = .data$x, colour = .data$gate, linetype = .data$gate),
      linewidth = 0.6
    ) +
    .omip111FigureGateScales(methods, methodColours) +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = -Inf, y = Inf, label = .data$label, vjust = 1.3 + 1.25 * .data$line),
      colour = labels$colour, inherit.aes = FALSE, hjust = -0.04, size = 2.3
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), ncol = 5, scales = "free") +
    ggplot2::scale_x_continuous(breaks = asinh(raw / cofactor),
      labels = c("-10³", "0", "10³", "10⁴", "10⁵", "10⁶")) +
    ggplot2::scale_y_continuous(trans = "sqrt") +
    ggplot2::labs(x = "Cytokine", y = "Cells (square-root scale)") +
    .analysis_theme(grid = "none") +
    ggplot2::theme(
      legend.position = "bottom", strip.text = ggplot2::element_text(size = 7),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey40"),
      axis.text = ggplot2::element_text(size = 6)
    )
}

.omip111FigureScales <- function(cofactor) {
  raw <- c(-1e3, 0, 1e3, 1e4, 1e5, 1e6)
  breaks <- asinh(raw / cofactor)
  labels <- c("-10³", "0", "10³", "10⁴", "10⁵", "10⁶")
  list(
    ggplot2::facet_grid(ggplot2::vars(.data$row), ggplot2::vars(.data$column),
      scales = "free_y", switch = "y"
    ),
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels),
    ggplot2::scale_y_continuous(breaks = breaks, labels = labels),
    ggplot2::labs(x = "CD44 RB744", y = NULL),
    .analysis_theme(grid = "none"),
    ggplot2::theme(
      strip.placement = "outside", panel.border = ggplot2::element_rect(fill = NA, colour = "grey40"),
      axis.text = ggplot2::element_text(size = 7)
    )
  )
}

.omip111AuthorFigure <- function(figureData, cofactor) {
  gates <- figureData$gates
  ggplot2::ggplot(figureData$cells, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$density), size = 0.12, shape = 16) +
    ggplot2::geom_rect(
      data = gates, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      inherit.aes = FALSE, fill = NA, colour = "red", linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = gates, ggplot2::aes(x = .data$xmin, y = .data$ymax, label = .data$label),
      inherit.aes = FALSE, colour = "red", hjust = 1.1, vjust = 1.2, size = 3
    ) +
    ggplot2::scale_colour_gradientn(
      colours = c("#0000FF", "#00BFFF", "#00FF00", "#FFFF00", "#FF0000"),
      trans = "log10", guide = "none"
    ) +
    .omip111FigureScales(cofactor)
}

