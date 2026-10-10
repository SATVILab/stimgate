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
    if (any(abs(references$importedParent / references$flowjoParent - 1) >
      preprocessing$parentRelativeTolerance)) {
      stop("Parent count import exceeds recorded tolerance.")
    }
    references$frequencyDifferencePp <- 100 * (references$importedCount / references$importedParent -
      references$flowjoCount / references$flowjoParent)
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
        expression <- asinh(rawExpression[indices, unname(markers), drop = FALSE] / preprocessing$cofactor)
        colnames(expression) <- names(markers)
        expression[] <- .omip111Float32(as.numeric(expression))
        manual <- as.matrix(membership[names(markers)]) == 1L
        if (!all(is.finite(expression))) stop("Non-finite processed expression.")
        refs <- references[references$sample == sid & references$population == population, ]
        if (any(colSums(manual) != refs$projectedCount[match(names(markers), refs$marker)])) {
          stop("Manual membership/reference counts disagree.")
        }
        directReference <- sweep(rawExpression[indices, unname(markers), drop = FALSE], 2L,
          refs$cytokineLowerRaw[match(names(markers), refs$marker)],
          FUN = ">"
        )
        if (!identical(unname(directReference), unname(manual))) {
          stop("Python reference masks do not match the raw FCS values read by R.")
        }
        author <- as.matrix(membership[paste0("author_", names(markers))]) == 1L
        colnames(author) <- names(markers)
        sampleInput[[population]] <- list(
          expression = expression, manual = manual,
          author = author, eventIndex = indices
        )
      }
      input[[sid]] <- sampleInput
    }
    manifest <- list(
      semantics = "omip111-v4", sourceSha256 = fingerprint,
      importer = importer, importerSourceSha256 = digest::digest(file = pathImporter, algo = "sha256"),
      preprocessing = preprocessing,
      scale = paste0("asinh(raw unmixed fluorescence / ", preprocessing$cofactor, "); float32"),
      compensation = "none; source FCS already spectrally unmixed",
      populations = "reconstructed author-gated non-naive conventional CD4/CD8 T cells"
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
        markers <- names(.omip111Markers(population))
        data <- stats::setNames(lapply(sampleMap$sample, function(sid) {
          prepared$input[[sid]][[population]]$expression
        }), sampleMap$sample)
        pathStim <- file.path(pathTmp, strain, population, "stimgate")
        message("Running StimGate: ", strain, " ", population)
        control <- do.call(stimgate::stimControl, settings$control)
        .analysis_with_seed(settings$seed, stimgate::gateStim(
          pathProject = pathStim, .data = data, batchList = batchList,
          marker = markers, biasUns = NULL, control = control
        ))
        gates <- stimgate::getStimGates(pathStim)
        gates <- gates[gates$gateName == "loc_minClust", ]
        stats <- stimgate::getStimStats(pathStim)
        stats <- stats[stats$gateName == "loc_minClust", ]
        if (nrow(gates) != 5L * length(markers)) stop("Incomplete StimGate gate table.")
        for (batch in batchList) {
          iUns <- batch[[1]]
          iStim <- batch[[2]]
          sidUns <- sampleMap$sample[[iUns]]
          sidStim <- sampleMap$sample[[iStim]]
          xUns <- data[[iUns]]
          xStim <- data[[iStim]]
          refUns <- prepared$input[[sidUns]][[population]]$manual
          refStim <- prepared$input[[sidStim]][[population]]$manual
          for (marker in markers) {
            info <- data.frame(
              strain = strain, mouse = sampleMap$mouse[[iStim]],
              population = population, marker = marker, sampleStim = sidStim, sampleUns = sidUns
            )
            gate <- gates[as.integer(as.character(gates$ind)) == iStim & gates$marker == marker, ]
            if (nrow(gate) != 1L) stop("Missing/duplicate final StimGate gate.")
            marginal <- stats[as.integer(as.character(stats$ind)) == iStim &
              grepl(paste0(marker, "~+~"), stats$cytCombn, fixed = TRUE), ]
            if (nrow(marginal) != 2^(length(markers) - 1L) ||
              any(marginal$nCellStim != nrow(xStim)) ||
              any(marginal$nCellUns != nrow(xUns))) {
              stop("Incomplete or inconsistent StimGate combination counts.")
            }
            # Marginalise the actual package combination counts: cytokine-positive
            # refinement is conditional, so a single base threshold is insufficient.
            counts <- .omip111Counts(
              sum(marginal$countStim), nrow(xStim),
              sum(marginal$countUns), nrow(xUns), sum(refStim[, marker]), sum(refUns[, marker])
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
                if (method == "F-beta") "fbeta_default" else "tailgate_default"
              )
              threshold <- .acsCytofThresholdOne(
                method = method, xUns = xUns[, marker], xStim = xStim[, marker],
                settings = comparator, pathFbeta = pathFbeta, fbetaEnv = fbetaEnv
              )
              counts <- .omip111Counts(
                if (is.finite(threshold$threshold)) sum(xStim[, marker] > threshold$threshold) else NA_real_,
                nrow(xStim),
                if (is.finite(threshold$threshold)) sum(xUns[, marker] > threshold$threshold) else NA_real_,
                nrow(xUns), sum(refStim[, marker]), sum(refUns[, marker])
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
    results <- dplyr::bind_rows(resultRows)
    writeLines(capture.output(utils::sessionInfo()), file.path(pathTmp, "session-info.txt"))
    originals <- prepared$references |>
      dplyr::mutate(flowjoProp = .data$flowjoCount / .data$flowjoParent) |>
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
        author2dDifferencePp = 100 * (.data$propBs - .data$flowjoBs),
        projectedVsAuthorDifferencePp = 100 * (.data$manualBs - .data$flowjoBs)
      )
    results$errorPp <- 100 * (results$propBs - results$manualBs)
    results$authorReferenceAccepted <- results$referenceAcceptedStim & results$referenceAcceptedUns
    # The primary reference is a raw cytokine lower cutoff applied directly to
    # the common parent cells; the 2-D CD44 reconstruction audit is contextual.
    results$referenceAccepted <- TRUE
    results$unscreenedErrorPp <- results$errorPp
    results$errorPp[!results$referenceAccepted] <- NA_real_
    results$absoluteErrorPp <- abs(results$errorPp)
    key <- results[c("strain", "mouse", "population", "marker", "method")]
    if (nrow(results) != 240L || anyDuplicated(key)) stop("Incomplete comparison cohort.")
    saveRDS(
      list(results = results, settings = settings, preprocessing = prepared$manifest),
      file.path(pathTmp, "result.rds")
    )
    utils::write.csv(results, file.path(pathTmp, "comparison.csv"), row.names = FALSE)
    writeLines("Complete: 240 method/sample/population/marker outcomes.", file.path(pathTmp, "COMPLETE"))
  })
  invisible(pathResults)
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
