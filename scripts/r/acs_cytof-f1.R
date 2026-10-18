# A frequency-derived reference on the same stimulated cells as every method.
# Strict x > gate: tied expression values cannot be split. Among equally close
# attainable counts choose the lowest threshold (the larger positive set).
.acsCytofF1One <- function(x, freqManual, thresholds, failed) {
  if (any(!is.finite(x))) stop("Non-finite ACS stimulated expression.")
  if (length(freqManual) != 1L) stop("ACS F1 requires one manual net frequency.")
  referenceReason <- if (!is.finite(freqManual)) {
    "missing_manual_frequency"
  } else if (freqManual < 0) {
    "negative_manual_net_frequency"
  } else if (freqManual > 100) {
    "manual_net_frequency_above_100"
  } else if (!length(x)) {
    "empty_stimulated_tube"
  } else {
    "valid"
  }
  referenceValid <- referenceReason == "valid"
  if (length(thresholds) != length(failed) || anyNA(failed)) {
    stop("Missing ACS gate failure provenance.")
  }
  n <- length(x)
  referenceGate <- NA_real_
  reference <- rep(FALSE, n)
  if (referenceValid) {
    runs <- rle(sort(x))
    gates <- c(-Inf, runs$values)
    counts <- c(n, n - cumsum(runs$lengths))
    distance <- abs(counts - n * freqManual / 100)
    closest <- which(distance <= min(distance) + 1e-10)[[1]]
    referenceGate <- gates[[closest]]
    reference <- x > referenceGate
  }
  nReference <- if (referenceValid) sum(reference) else NA_integer_
  dplyr::bind_rows(lapply(seq_along(thresholds), function(i) {
    valid <- !failed[[i]] && is.finite(thresholds[[i]]) && n > 0L
    predicted <- if (valid) x > thresholds[[i]] else rep(NA, n)
    tp <- if (valid && referenceValid) sum(predicted & reference) else NA_integer_
    np <- if (valid) sum(predicted) else NA_integer_
    denominator <- np + nReference
    tibble::tibble(
      referenceGate = referenceGate, referenceReason = referenceReason,
      freqManualNet = freqManual, nCell = n, nReference = nReference,
      freqReference = if (n) 100 * nReference / n else NA_real_,
      referenceFrequencyError = if (n) 100 * nReference / n - freqManual else NA_real_,
      nPredicted = np, tp = tp, fp = np - tp, fn = nReference - tp,
      f1 = if (valid && referenceValid && denominator > 0) 2 * tp / denominator else NA_real_
    )
  }))
}

# Counts and frequencies must describe the same saved applied gate.
.acsCytofF1ValidateCounts <- function(d, scores, n) {
  if (anyNA(d$nCellStim) || any(d$nCellStim != n)) {
    stop("ACS F1 cell count differs from saved statistics.")
  }
  # Some valid per-marker Tailgate gates have no combination-derived count.
  # Retain their cell-based score and explicitly report the unavailable audit.
  valid <- is.finite(scores$nPredicted) & is.finite(d$countStim)
  if (any(scores$nPredicted[valid] != d$countStim[valid])) {
    stop("ACS F1 gates do not reproduce saved positive counts.")
  }
  error <- if (n) 100 * scores$nPredicted / n - d$freq_stim_auto else rep(NA_real_, nrow(d))
  if (any(!is.finite(error[valid])) || any(abs(error[valid]) > 1e-8)) {
    stop("ACS F1 gates do not reproduce saved stimulated frequencies.")
  }
  error[!valid] <- NA_real_
  error
}

# Read the prepared population GatingSets; never re-gate or sample cells.
.acsCytofF1Table <- function(comparisonTbl, pathGsBase) {
  .acsCytofValidationValidateComparisonTable(comparisonTbl)
  required <- c("popCode", "ind", "threshold", "thresholdFailed", "nCellStim", "countStim",
                "freq_stim_auto", "freq_stim_man", "freq_uns_man")
  if (!all(required %in% names(comparisonTbl))) {
    stop("ACS F1 requires saved gates and failure provenance; re-run Analysis 9b.")
  }
  manifests <- attr(comparisonTbl, "manifest")$methods
  dplyr::bind_rows(lapply(unique(as.character(comparisonTbl$popCode)), function(popCode) {
    path <- file.path(pathGsBase, popCode)
    gs <- flowWorkspace::load_gs(path, backend_readonly = TRUE)
    preprocessing <- .acsCytofReadPreprocessing(path, gs)
    saved <- manifests[[popCode]]
    .acsCytofValidateManifests(saved)
    if (!identical(preprocessing, saved[[1]]$context$preprocessing)) {
      stop("ACS F1 GatingSet differs from the saved method inputs: ", popCode)
    }
    rows <- comparisonTbl[as.character(comparisonTbl$popCode) == popCode, ]
    rows |>
      dplyr::group_by(.data$ind, .data$cyt) |>
      dplyr::group_modify(function(d, key) {
        map <- preprocessing$sampleMap
        tube <- map[as.character(map$ind) == as.character(key$ind), ]
        channel <- names(.acsCytofChannelMap())[.acsCytofChannelMap() == as.character(key$cyt)]
        if (nrow(tube) != 1L || length(channel) != 1L ||
            any(d$SampleID != tube$SampleID) || any(d$stim != tube$stim) ||
            length(unique(d$freq_stim_man)) != 1L ||
            length(unique(d$freq_uns_man)) != 1L) {
          stop("Inconsistent ACS F1 tube/manual reference keys.")
        }
        x <- .acs_get_expr(gs, key$ind, channel)
        scores <- .acsCytofF1One(
          x, d$freq_stim_man[[1]] - d$freq_uns_man[[1]], d$threshold, d$thresholdFailed
        )
        scores$savedCountAvailable <- is.finite(d$countStim)
        scores$predictedFrequencyError <- .acsCytofF1ValidateCounts(d, scores, length(x))
        dplyr::bind_cols(d, scores)
      }) |>
      dplyr::ungroup()
  }))
}

# Equal weight per matched tube/donor, with the existing manual signal rules.
.acsCytofF1Summary <- function(table) {
  out <- table |>
    dplyr::group_by(.data$method, .data$pop, .data$cyt, .data$stim) |>
    dplyr::summarise(
      eligible = .acsCytofValidationEligible(
        .data$freq_stim_man, .data$freq_uns_man, .data$freq_bs_man
      ),
      n_total = dplyr::n(),
      n_saved_count_missing = sum(!.data$savedCountAvailable), n_failed = sum(.data$thresholdFailed),
      n_scored = sum(is.finite(.data$f1)),
      n_invalid_reference = sum(.data$referenceReason != "valid"),
      n_undefined = sum(!.data$thresholdFailed & !is.finite(.data$f1)),
      f1_median = if (any(is.finite(.data$f1))) stats::median(.data$f1[is.finite(.data$f1)]) else NA_real_,
      f1_q10 = if (any(is.finite(.data$f1))) {
        unname(stats::quantile(.data$f1[is.finite(.data$f1)], 0.1, type = 7))
      } else NA_real_,
      f1 = if (any(is.finite(.data$f1))) mean(.data$f1[is.finite(.data$f1)]) else NA_real_,
      .groups = "drop"
    )
  excluded <- out[!out$eligible, ]
  excluded[c("f1", "f1_median", "f1_q10")] <- NA_real_
  out <- out[out$eligible, ]
  attr(out, "excludedStrata") <- excluded
  out
}
