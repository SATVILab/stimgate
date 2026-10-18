# Counts use unrounded scores and eligibility already computed by the ACS
# correlation or F1 summary helpers, exactly as in the thesis heatmaps.
.acsCytofCccSummary <- function(correlationTbl, tolerance = 1e-12, metric = c("ccc", "f1")) {
  metric <- match.arg(metric)
  methods <- c("stimgate", "fbeta_default", "tailgate")
  keys <- c("pop", "cyt", "stim")
  required <- c("method", keys, metric)
  if (!is.data.frame(correlationTbl) || !all(required %in% names(correlationTbl)) ||
      !is.numeric(correlationTbl[[metric]])) {
    stop("ACS summary requires method, pop, cyt, stim and numeric ", metric, ".")
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("ACS margin tolerance must be one finite nonnegative number.")
  }
  displayed <- function(tbl) {
    tbl[tbl$method %in% methods &
      tbl$pop %in% c(.acsCytofValidationRealPopulations(), "B cells", "NK cells") &
      tbl$cyt %in% c("IFNg", "IL2", "TNF", "IL17", "IL22", "IL6") &
      tbl$stim %in% c("p1", "mtb", "ebv"), , drop = FALSE]
  }
  eligible <- displayed(as.data.frame(correlationTbl))
  if (anyNA(eligible[c("method", keys)]) || anyDuplicated(eligible[c("method", keys)])) {
    stop("ACS summary requires unique nonmissing method/population/cytokine/stimulus keys.")
  }
  strata <- unique(eligible[keys])
  excluded <- attr(correlationTbl, "excludedStrata")
  observed <- strata
  if (!is.null(excluded)) {
    observed <- unique(rbind(observed, displayed(as.data.frame(excluded))[keys]))
  }
  for (method in methods) {
    part <- eligible[eligible$method == method, c(keys, metric), drop = FALSE]
    names(part)[names(part) == metric] <- method
    strata <- merge(strata, part, by = keys, all.x = TRUE, sort = FALSE)
  }
  values <- as.matrix(strata[methods])
  complete <- rowSums(is.finite(values)) == length(methods)
  comparison <- values[complete, , drop = FALSE]
  metrics <- c(paste0("best_", metric), "win_one", "win_both", "lose_0_1", "lose_0_25", "fail")
  counts <- data.frame(method = methods, stringsAsFactors = FALSE)
  for (column in metrics) counts[[column]] <- rep(NA_integer_, length(methods))
  ties <- 0L
  if (nrow(comparison)) {
    maxima <- apply(comparison, 1L, max)
    ties <- sum(rowSums(comparison == maxima) > 1L)
    for (j in seq_along(methods)) {
      own <- comparison[, j]
      others <- comparison[, -j, drop = FALSE]
      margins <- own - others
      best_other <- apply(others, 1L, max)
      counts[j, metrics] <- list(
        sum(own == maxima),
        sum(own == maxima & rowSums(margins >= 0.1 - tolerance) >= 1L),
        sum(own == maxima & rowSums(margins >= 0.1 - tolerance) == 2L),
        sum(best_other - own >= 0.1 - tolerance),
        sum(best_other - own >= 0.25 - tolerance),
        sum(own < 0.25 & best_other >= 0.75)
      )
    }
  }
  strata$complete_three_finite <- complete
  coverage <- data.frame(
    displayed_observed = nrow(observed), eligible = nrow(strata),
    complete_three_finite = sum(complete), incomplete = sum(!complete),
    exact_top_ties = ties
  )
  methodCoverage <- data.frame(
    method = methods,
    eligible = vapply(methods, function(m) sum(eligible$method == m), integer(1)),
    finite_ccc = colSums(is.finite(values)),
    missing_or_nonfinite_ccc = nrow(strata) - colSums(is.finite(values))
  )
  names(methodCoverage)[3:4] <- c(paste0("finite_", metric), paste0("missing_or_nonfinite_", metric))
  list(counts = counts, coverage = coverage, method_coverage = methodCoverage,
       strata = strata, margin_tolerance = tolerance)
}

.acsCytofCccSummaryMarkdown <- function(summary) {
  counts <- summary$counts
  coverage <- summary$coverage
  metric <- if ("best_f1" %in% names(counts)) "f1" else "ccc"
  score <- if (metric == "f1") "F1" else "CCC"
  methods <- c("stimgate", "fbeta_default", "tailgate")
  counts <- counts[match(methods, counts$method), , drop = FALSE]
  cells <- cbind(
    c("StimGate", "$F_{\\beta}$ (all cells)", "Tailgate"),
    vapply(counts[-1L], .analysis_md_num, character(nrow(counts)), digits = 0, big.mark = ",")
  )
  caption <- c(
    sprintf(if (metric == "ccc") "Comparison of concordance correlation coefficients (CCC) across %s population, cytokine and stimulation combinations." else
      "Comparison of mean finite per-tube F1 scores across %s population, cytokine and stimulation combinations.",
      .analysis_md_num(coverage$complete_three_finite, 0, big.mark = ",")),
    sprintf("Best %s counts combinations in which each method achieves the highest %s, including exact ties.", score, score),
    paste0("Win columns count combinations in which a method has the best ", score, " that is higher than ",
      "at least one (unspecified) or both (\"(both)\") competitors by at least 0.1."),
    "Loss columns count combinations in which a method is lower than the best competitor by at least 0.1 or 0.25, respectively.",
    sprintf("Failure counts combinations in which a method achieves %s $< 0.25$ when at least one competing method achieves %s $\\geq 0.75$.", score, score),
    paste0("Only eligible displayed strata (p4 excluded) with three finite",
      if (metric == "ccc") ", signed CCCs" else " F1 scores", " are compared; ",
      if (summary$margin_tolerance > 0) paste0("margin tolerance is $",
        .analysis_md_pow10(summary$margin_tolerance), "$.") else "no margin tolerance is applied.")
  )
  if (metric == "f1") {
    caption <- c(caption, paste0("Per-tube F1 uses the strict upper-tail reference closest to the manual ",
      "background-subtracted frequency. It compares the reference and method cell masks ",
      "on the same stimulated cells, with no further background subtraction; manual event labels are unavailable."))
  }
  if (coverage$incomplete || !coverage$complete_three_finite) {
    caption <- c(caption, sprintf("Of %s eligible strata, %s incomplete strata are omitted; counts are NA if none are complete.",
      .analysis_md_num(coverage$eligible, 0, big.mark = ","),
      .analysis_md_num(coverage$incomplete, 0, big.mark = ",")))
  }
  .analysis_md_table(
    cells,
    c("Method", "Best", "Win $\\geq 0.1$", "Win $\\geq 0.1$ (both)",
      "Loss $\\geq 0.1$", "Loss $\\geq 0.25$", "Fail"),
    caption, paste0("tbl-stimgate-results-acs-", metric, "-wins"), align = c("left", rep("right", 6))
  )
}
