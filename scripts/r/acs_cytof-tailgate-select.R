# Analysis 9a (`9a-real-tune-acs-cytof-tailgate.qmd`): score the eight
# Tailgate configurations of `.acsCytofTailgateConfigs()` against the ACS
# manual frequencies and select one by mean Pearson correlation (PCC).
#
# Scoring cohort: every population/cytokine/stimulus stratum that passes the
# manuscript signal eligibility rule (`.acsCytofValidationEligible()`), without
# the non-secreted-protein stimulus `p4`, which the heatmaps also omit. Within
# a stratum every configuration is scored on the same donors: those with a
# finite manual frequency and a successful gate under every configuration. A
# stratum enters the mean only if every configuration's PCC is defined there,
# so dropping difficult strata cannot raise one configuration's mean. Failed
# gates stay missing; they are never scored as zero frequencies.

.acsCytofTailgateSelectionSemantics <- function() {
  "acs-tailgate-selection-v2"
}

# One row per stratum and configuration, including ineligible strata
# (`eligible = FALSE`, no statistics) so heatmaps can mark them.
.acsCytofTailgateScoreStrata <- function(
    comparisonTbl, methods = .acsCytofTailgateSelectionMethods()) {
  .acsCytofValidationValidateComparisonTable(comparisonTbl, requiredMethods = methods)
  if (!"thresholdFailed" %in% names(comparisonTbl)) {
    stop("ACS comparison rows lack thresholdFailed provenance.")
  }
  tbl <- comparisonTbl |>
    dplyr::mutate(
      method = as.character(.data$method), pop = as.character(.data$pop),
      cyt = as.character(.data$cyt), stim = as.character(.data$stim)
    ) |>
    dplyr::filter(.data$method %in% methods, .data$stim != "p4") |>
    dplyr::mutate(freq_bs_auto = dplyr::if_else(
      .data$thresholdFailed %in% FALSE, .data$freq_bs_auto, NA_real_
    ))

  tbl |>
    dplyr::group_by(.data$pop, .data$cyt, .data$stim) |>
    dplyr::group_modify(function(rows, key) {
      manual <- rows[rows$method == methods[[1]], ]
      eligible <- .acsCytofValidationEligible(
        manual$freq_stim_man, manual$freq_uns_man, manual$freq_bs_man
      )
      finiteAll <- rows |>
        dplyr::group_by(.data$SampleID) |>
        dplyr::summarise(
          common = all(is.finite(.data$freq_bs_auto)) &&
            all(is.finite(.data$freq_bs_man)) &&
            dplyr::n() == length(methods),
          .groups = "drop"
        )
      commonIds <- finiteAll$SampleID[finiteAll$common]
      purrr::map_dfr(methods, function(m) {
        r <- rows[rows$method == m, ]
        own <- is.finite(r$freq_bs_auto) & is.finite(r$freq_bs_man)
        com <- r$SampleID %in% commonIds
        tibble::tibble(
          method = m, eligible = eligible,
          n_rows = nrow(r), n_failed = sum(!(r$thresholdFailed %in% FALSE)),
          n_finite = sum(own), n_common = sum(com),
          pcc = if (eligible) .acsCytofTailgatePcc(r$freq_bs_auto[com], r$freq_bs_man[com]) else NA_real_,
          ccc = if (eligible) .acsCytofValidationCcc(r$freq_bs_auto[com], r$freq_bs_man[com]) else NA_real_,
          mean_abs_diff = if (eligible && any(com)) mean(abs(r$freq_bs_auto[com] - r$freq_bs_man[com])) else NA_real_,
          pcc_own = if (eligible) .acsCytofTailgatePcc(r$freq_bs_auto[own], r$freq_bs_man[own]) else NA_real_
        )
      })
    }) |>
    dplyr::ungroup()
}

# NA unless at least two pairs and both variables vary.
.acsCytofTailgatePcc <- function(x, y) {
  if (length(x) < 2L) return(NA_real_)
  out <- suppressWarnings(stats::cor(x, y))
  if (is.finite(out)) out else NA_real_
}

# Mean PCC per configuration over strata defined for every configuration,
# with coverage, failure and diagnostic columns, ranked best first.
.acsCytofTailgateRanking <- function(strataTbl) {
  eligible <- strataTbl[strataTbl$eligible, ]
  scored <- eligible |>
    dplyr::group_by(.data$pop, .data$cyt, .data$stim) |>
    dplyr::summarise(scored = all(is.finite(.data$pcc)), .groups = "drop")
  eligible <- dplyr::left_join(eligible, scored, by = c("pop", "cyt", "stim"))
  configs <- .acsCytofTailgateConfigs()
  ranking <- eligible |>
    dplyr::group_by(.data$method) |>
    dplyr::summarise(
      n_eligible = dplyr::n(),
      n_defined = sum(is.finite(.data$pcc)),
      n_undefined = sum(!is.finite(.data$pcc)),
      n_scored = sum(.data$scored),
      mean_pcc = if (any(.data$scored)) mean(.data$pcc[.data$scored]) else NA_real_,
      mean_ccc = if (any(.data$scored)) mean(.data$ccc[.data$scored]) else NA_real_,
      mean_abs_diff = if (any(.data$scored)) mean(.data$mean_abs_diff[.data$scored]) else NA_real_,
      mean_pcc_own = if (any(is.finite(.data$pcc_own))) mean(.data$pcc_own, na.rm = TRUE) else NA_real_,
      n_pcc_own = sum(is.finite(.data$pcc_own)),
      n_rows = sum(.data$n_rows),
      n_failed = sum(.data$n_failed),
      n_common = sum(.data$n_common),
      .groups = "drop"
    )
  ranking <- configs |>
    dplyr::inner_join(ranking, by = "method") |>
    dplyr::arrange(dplyr::desc(.data$mean_pcc), match(.data$method, configs$method))
  top <- max(ranking$mean_pcc, -Inf, na.rm = TRUE)
  ranking$rank <- rank(-ranking$mean_pcc, ties.method = "min", na.last = "keep")
  ranking$best <- is.finite(ranking$mean_pcc) & ranking$mean_pcc == top
  ranking[c("rank", setdiff(names(ranking), "rank"))]
}

# "selected" only when exactly one configuration has the highest finite mean.
.acsCytofTailgateSelectionStatus <- function(ranking) {
  if (!any(ranking$best)) return("unavailable")
  if (sum(ranking$best) > 1L) return("tie")
  "selected"
}

# Correlation table in the shape .acsCytofValidationPlotCorrelation() takes,
# using the common-donor PCC and CCC that the ranking scores.
.acsCytofTailgateCorrelationTable <- function(strataTbl) {
  keys <- c("method", "pop", "cyt", "stim")
  out <- strataTbl[strataTbl$eligible, c(keys, "n_common", "pcc", "ccc")]
  attr(out, "excludedStrata") <- strataTbl[!strataTbl$eligible, keys]
  out
}

.acsCytofTailgateSelectionRecord <- function(
    ranking, comparisonTbl, runId = NA_character_) {
  status <- .acsCytofTailgateSelectionStatus(ranking)
  best <- ranking$method[ranking$best]
  method <- if (identical(status, "selected")) best else NA_character_
  configs <- .acsCytofTailgateConfigs()
  manifest <- attr(comparisonTbl, "manifest")
  list(
    status = status,
    method = method,
    tied = if (identical(status, "tie")) best else character(),
    config = if (!is.na(method)) as.list(configs[configs$method == method, ]) else NULL,
    params = if (!is.na(method)) .acsCytofComparatorSettings(method)$params else NULL,
    criterion = paste(
      "Highest arithmetic mean PCC over eligible population x cytokine x",
      "stimulus strata (p4 omitted), each stratum weighted equally, scored",
      "on donors with successful gates under every configuration."
    ),
    ranking = ranking,
    configs = configs,
    cytoUtils = list(
      version = as.character(utils::packageVersion("cytoUtils")),
      sha = utils::packageDescription("cytoUtils")$RemoteSha
    ),
    # One preprocessing context per population, without the git SHA; 9b
    # requires its own comparison to have identical contexts.
    preprocessing = .acsCytofTailgateContexts(manifest$methods),
    manualInputHash = manifest$manualInputHash,
    analysis_semantics_version = .acsCytofTailgateSelectionSemantics(),
    run_id = runId,
    gitSha = .git_sha(),
    createdAt = Sys.time()
  )
}

.acsCytofTailgateContexts <- function(methodManifests) {
  lapply(methodManifests, function(population) {
    context <- population[[1]]$context
    context$gitSha <- NULL
    context
  })
}

# Analysis 9b's comparison must use the GatingSets Analysis 9a selected on.
.acsCytofCheckSelectionPreprocessing <- function(selection, comparisonTbl) {
  contexts <- .acsCytofTailgateContexts(attr(comparisonTbl, "manifest")$methods)
  shared <- intersect(names(contexts), names(selection$preprocessing))
  if (!length(shared) || !setequal(names(contexts), shared) ||
      !identical(contexts[shared], selection$preprocessing[shared])) {
    stop(
      "Analysis 9b's populations or preprocessing differ from those Analysis 9a ",
      "selected Tailgate's configuration on; re-render Analysis 9a, then 9b."
    )
  }
  invisible(TRUE)
}

# Display formatting for the thesis table.
.acsCytofTailgateFormat <- function(ranking) {
  tibble::tibble(
    rank = as.character(ranking$rank),
    configuration = ranking$configLabel,
    zeros = ifelse(ranking$removeZero, "Yes", "No"),
    auto_tol = ifelse(ranking$autoTol, "Yes", "No"),
    tol = vapply(ranking$tol, .acsCytofTailgateNumber, character(1)),
    bias = format(ranking$bias, trim = TRUE, drop0trailing = TRUE),
    mean_pcc = ifelse(is.finite(ranking$mean_pcc),
      sub("^-", "−", formatC(ranking$mean_pcc, format = "f", digits = 3)), "NA"),
    strata = paste0(ranking$n_scored, "/", ranking$n_eligible),
    failed = paste0(ranking$n_failed, "/", ranking$n_rows),
    selected = ranking$best & sum(ranking$best) == 1L
  )
}

# 0.01 stays decimal; small values use scientific notation.
.acsCytofTailgateNumber <- function(x) {
  if (x >= 1e-3) return(format(x, scientific = FALSE, drop0trailing = TRUE))
  exponent <- floor(log10(x))
  paste0(format(signif(x / 10^exponent, 3)), "e", exponent)
}

# The ranking as a Quarto pipe table with a cross-referenceable caption
# (needs analysis-thesis-table.R); the HTML renders the same source.
.acsCytofTailgateMdTable <- function(ranking, label = "tbl-stimgate-appendix-acs-tailgate") {
  f <- .acsCytofTailgateFormat(ranking)
  tol <- ifelse(
    grepl("e", f$tol),
    sub("^([0-9.]+)e(-?[0-9]+)$", "$\\1\\\\times10^{\\2}$", f$tol),
    f$tol
  )
  configuration <- .analysis_md_escape(f$configuration)
  configuration[f$selected] <- paste0(configuration[f$selected], " (selected)")
  selected <- f$configuration[f$selected]
  selectedText <- if (length(selected)) {
    paste0("The selected configuration is in bold. ")
  } else {
    "No single configuration has the highest mean PCC, so none is selected. "
  }
  .analysis_md_table(
    cbind(f$rank, configuration, f$zeros, f$auto_tol, tol, f$bias, f$mean_pcc,
      f$strata, f$failed),
    header = c("Rank", "Configuration", "Zeros removed", "Automatic tolerance",
      "Tolerance", "Bias", "Mean PCC", "Strata scored", "Failed gates"),
    align = c("right", "left", "left", "left", rep("right", 5)),
    bold = f$selected,
    caption = paste0(
      "**Choosing Tailgate's configuration on the ACS CyTOF data**. Tailgate configurations ",
      "ranked by mean Pearson correlation (PCC) with the manual background-subtracted ",
      "frequencies over eligible population, cytokine and stimulus strata, each stratum ",
      "weighted equally. ", selectedText,
      "With automatic tolerance, the tolerance is a multiple of the largest absolute first ",
      "derivative of the density; otherwise it is an absolute bound on the derivative. ",
      "Zeros (cells with an expression of exactly zero) are removed only to estimate the ",
      "threshold. Strata scored: eligible strata where every configuration has a PCC on the ",
      "donors all gated successfully. Failed gates: failed donor-cytokine estimates in ",
      "eligible strata."
    ),
    label = label
  )
}
