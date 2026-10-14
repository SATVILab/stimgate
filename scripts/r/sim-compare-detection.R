# Sample-level evidence from final applied gates, independent of simulation and
# display settings. Source after analysis-plot-style.R and analysis-mcse.R.
# The source-compatibility reader also requires sim-misc.R.
.simCompareDetectionCutoffs <- sort(unique(c(
  0.001, 0.005, 0.01, 0.025,
  0.05, 0.1, 0.2, 10^seq(-12, 0, length.out = 121)
)))

# P(X >= nPosStim | fixed margins), the greater alternative of Fisher's test.
# Work on the log scale: exp(log_p) may underflow but decisions and evidence do not.
.simCompareFisherScores <- function(data) {
  counts <- c("nPosStim", "nPosUns", "nCellStim", "nCellUns")
  if (!all(c(counts, "propRespEst") %in% names(data))) {
    stop("Detection requires saved final applied-gate counts and propRespEst.")
  }
  m <- as.matrix(data[counts])
  valid <- rowSums(!is.finite(m) | m < 0 | m != floor(m)) == 0L &
    data$nCellStim > 0 & data$nCellUns > 0 &
    data$nPosStim <= data$nCellStim & data$nPosUns <= data$nCellUns &
    is.finite(data$propRespEst)
  if ("error" %in% names(data)) {
    valid <- valid & (is.na(data$error) | !nzchar(trimws(data$error)))
  }
  valid[is.na(valid)] <- FALSE
  data$detection_valid <- valid
  data$fisher_log_p <- NA_real_
  data$fisher_log_p[valid] <- stats::phyper(
    data$nPosStim[valid] - 1, data$nPosStim[valid] + data$nPosUns[valid],
    data$nCellStim[valid] + data$nCellUns[valid] -
      data$nPosStim[valid] - data$nPosUns[valid],
    data$nCellStim[valid],
    lower.tail = FALSE, log.p = TRUE
  )
  if (any(!is.finite(data$fisher_log_p[valid]))) {
    stop("A valid finite count table produced a non-finite Fisher log tail.")
  }
  data$fisher_p <- exp(data$fisher_log_p)
  data$fisher_log10 <- -data$fisher_log_p / log(10)
  data$fisher_log2 <- -data$fisher_log_p / log(2)
  # SBB lower bound on B01; not an estimated or exact Bayes factor. Log form
  # avoids underflow and cancellation when p is tiny.
  data$sbb_log_b01_min <- ifelse(data$fisher_log_p < -1,
    1 + data$fisher_log_p + log(-data$fisher_log_p), 0
  )
  data$sbb_b01_min <- exp(data$sbb_log_b01_min)
  data$sbb_log10 <- -data$sbb_log_b01_min / log(10)
  .analysis_method_data(data)
}

.simCompareDetectionSummary <- function(
  scores, group_cols,
  cutoffs = .simCompareDetectionCutoffs, mcse = FALSE, retain_draws = FALSE
) {
  if (any(!is.finite(cutoffs) | cutoffs <= 0 | cutoffs > 1)) {
    stop("Detection cutoffs must be in (0, 1].")
  }
  scores |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(rows, key) {
      family <- unique(paste0("sim_seed:", if ("sim_seed" %in% names(rows)) rows$sim_seed else key$sim_seed))
      if (length(family) != 1L) stop("Detection scenario must have one bootstrap family.")
      ids <- sort(unique(as.character(rows$iter)))
      valid <- vapply(ids, function(id) sum(rows$detection_valid[as.character(rows$iter) == id]), numeric(1))
      # Exact pooled-sample bootstrap via dataset multiplicities, computed once
      # for all cutoffs. This equals expanding each sampled dataset's rows.
      weights <- if (isTRUE(mcse)) {
        indices <- .analysis_mcse_bootstrap_indices(rows$iter, family)
        matrix(vapply(
          seq_len(ncol(indices)), function(i) tabulate(indices[, i], nbins = length(ids)),
          integer(length(ids))
        ), nrow = length(ids))
      } else {
        matrix(numeric(), length(ids), 0L)
      }
      denominators <- if (isTRUE(mcse)) as.numeric(valid %*% weights) else numeric()
      out <- dplyr::bind_rows(lapply(sort(unique(cutoffs)), function(alpha) {
        detected <- rows$detection_valid & rows$propRespEst > 0 & rows$fisher_log_p <= log(alpha)
        counts <- vapply(ids, function(id) sum(detected[as.character(rows$iter) == id], na.rm = TRUE), numeric(1))
        estimate <- if (sum(valid)) sum(counts) / sum(valid) else NA_real_
        draws <- if (isTRUE(mcse)) {
          ifelse(denominators > 0,
            as.numeric(counts %*% weights) / denominators, NA_real_
          )
        } else {
          numeric()
        }
        dplyr::bind_cols(
          tibble::tibble(
            alpha = alpha,
            n_sample = nrow(rows), n_dataset = length(ids),
            n_valid = sum(valid), n_missing = sum(!rows$detection_valid),
            n_missing_estimate = sum(!is.finite(rows$propRespEst)),
            n_detected = sum(counts)
          ),
          .analysis_mcse_bootstrap_cols(estimate, draws, "detection_rate",
            length(ids), sum(valid > 0),
            mcse = mcse
          )
        )
      }))
      if (!isTRUE(retain_draws)) out <- out[, !startsWith(names(out), ".boot_"), drop = FALSE]
      out
    }) |>
    dplyr::ungroup()
}

.simCompareEvidenceSummary <- function(scores, group_cols, mcse = FALSE) {
  scores |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(rows, key) {
      family <- unique(paste0("sim_seed:", if ("sim_seed" %in% names(rows)) rows$sim_seed else key$sim_seed))
      if (length(family) != 1L) stop("Evidence scenario must have one bootstrap family.")
      out <- tibble::tibble(
        n_sample = nrow(rows),
        n_valid = sum(rows$detection_valid), n_missing = sum(!rows$detection_valid)
      )
      for (score in c("fisher_log10", "fisher_log2", "sbb_log10")) {
        for (p in c(0.5, 0.9)) {
          nm <- paste0(score, if (p == 0.5) "_median" else "_p90")
          cols <- .analysis_mcse_pooled_cols(
            rows[[score]], as.character(rows$iter),
            function(x) .analysis_mcse_quantile_finite(x, p), nm, family, mcse
          )
          out <- dplyr::bind_cols(out, cols[!startsWith(names(cols), ".boot_")])
        }
        out[[paste0("n_above_cap_", score)]] <- sum(rows[[score]] >
          .simCompareEvidenceCap(score), na.rm = TRUE)
      }
      out
    }) |>
    dplyr::ungroup()
}

.simCompareEvidenceCap <- function(score) {
  switch(score,
    fisher_log10 = 12,
    fisher_log2 = 12 * log2(10),
    sbb_log10 = -(1 - 12 * log(10) + log(12 * log(10))) / log(10),
    stop("Unknown evidence scale.")
  )
}

.simCompareEvidenceReference <- function(score, alpha = 0.05) {
  switch(score,
    fisher_log10 = -log10(alpha),
    fisher_log2 = -log2(alpha),
    sbb_log10 = ifelse(alpha < 1 / exp(1), -log10(-exp(1) * alpha * log(alpha)), 0)
  )
}

# Match biological scenarios explicitly, including requested and realised
# background. Independent analyses keep separate bootstrap streams; methods
# within each analysis retain paired dataset draws. No interpolation at 5% FPR.
.simCompareDetectionRoc <- function(
  null_scores, response_scores, mcse = FALSE,
  cutoffs = .simCompareDetectionCutoffs
) {
  keys <- c(
    "transformation", "mean_pos_setting", "mean_pos", "n_cell",
    "sample_perturbation_sd", "condition_perturbation_sd", "cluster_perturbation_sd",
    "n_cell_uns_relative_to_stim"
  )
  required <- c(keys, "propUnsTruth", "propRespTruth", "prob_response")
  if (!all(c(required, "prob_background") %in% names(null_scores)) ||
    !all(c(required, "background_relative_to_response") %in% names(response_scores))) {
    stop("ROC matching requires biological settings and saved realised background truth.")
  }
  null_scores <- dplyr::mutate(null_scores,
    background_requested = .data$prob_background,
    background_realised = .data$propUnsTruth
  )
  response_scores <- dplyr::filter(response_scores, .data$propRespTruth > 0) |>
    dplyr::mutate(
      background_requested = .data$prob_response * .data$background_relative_to_response,
      background_realised = .data$propUnsTruth
    )
  join_keys <- c(keys, "background_requested", "background_realised")
  null_keys <- dplyr::distinct(null_scores, dplyr::across(dplyr::all_of(join_keys)))
  response_keys <- dplyr::distinct(
    response_scores,
    dplyr::across(dplyr::all_of(c(join_keys, "prob_response")))
  )
  audit <- dplyr::left_join(response_keys,
    dplyr::mutate(null_keys, matched = TRUE),
    by = join_keys
  ) |>
    dplyr::mutate(matched = .data$matched %in% TRUE)
  matched <- dplyr::filter(audit, .data$matched)
  if (!nrow(matched)) {
    return(list(audit = audit, curves = tibble::tibble(), at_5 = tibble::tibble()))
  }
  cols <- c(join_keys, "method")
  null <- .simCompareDetectionSummary(
    dplyr::semi_join(null_scores, matched, by = join_keys),
    cols, cutoffs, mcse
  )
  response <- .simCompareDetectionSummary(dplyr::semi_join(response_scores, matched,
    by = c(join_keys, "prob_response")
  ), c(cols, "prob_response"), cutoffs, mcse)
  curves <- dplyr::inner_join(response, null, by = c(cols, "alpha"), suffix = c("_response", "_null"))
  # A supported decision threshold, not an interpolated sensitivity. Boundary
  # cases (all cutoffs below/above 5%) cannot locate a 5% operating region.
  at_5 <- curves |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(cols, "prob_response")))) |>
    dplyr::group_modify(function(rows, key) {
      fpr <- rows$detection_rate_null
      bracket <- any(fpr <= 0.05, na.rm = TRUE) && any(fpr > 0.05, na.rm = TRUE)
      eligible <- which(is.finite(fpr) & fpr <= 0.05)
      chosen <- if (bracket && length(eligible)) eligible[which.max(rows$alpha[eligible])] else integer()
      tibble::tibble(
        region_bracketed = bracket,
        alpha_selected = if (length(chosen)) rows$alpha[chosen] else NA_real_,
        fpr_bracket_low = if (length(eligible)) max(fpr[eligible]) else NA_real_,
        fpr_bracket_high = if (any(fpr > 0.05, na.rm = TRUE)) min(fpr[fpr > 0.05 & is.finite(fpr)]) else NA_real_,
        fpr_observed = if (length(chosen)) fpr[chosen] else NA_real_,
        sensitivity_at_supported_cutoff = if (length(chosen)) rows$detection_rate_response[chosen] else NA_real_,
        fpr_lower = if (length(chosen)) rows$detection_rate_lower_null[chosen] else NA_real_,
        fpr_upper = if (length(chosen)) rows$detection_rate_upper_null[chosen] else NA_real_,
        sensitivity_lower = if (length(chosen)) rows$detection_rate_lower_response[chosen] else NA_real_,
        sensitivity_upper = if (length(chosen)) rows$detection_rate_upper_response[chosen] else NA_real_
      )
    }) |>
    dplyr::ungroup()
  list(audit = audit, curves = curves, at_5 = at_5)
}

# Source-derived compatibility checks for Analysis 7, like the Analysis 13
# settings reader. Evaluate settings/grid only; never run simulations or collation.
.simCompareDetectionResponseParams <- function(
  qmd_path, root_dir, simulation_seed,
  analysis_quick = FALSE, analysis_dev = FALSE, sim_size = "final"
) {
  lines <- readLines(qmd_path, warn = FALSE)
  chunk <- function(label) {
    start <- grep(paste0("^#\\| label: ", label, "$"), lines)
    if (length(start) != 1L) stop("Missing/ambiguous Analysis 7 chunk: ", label)
    end <- start + which(lines[-seq_len(start)] == "```")[1L]
    if (is.na(end)) stop("Unclosed Analysis 7 chunk: ", label)
    parse(text = lines[(start + 1L):(end - 1L)])
  }
  env <- new.env(parent = environment())
  env$root_dir <- root_dir
  env$analysis_quick <- analysis_quick
  env$analysis_dev <- analysis_dev
  env$sim_size <- sim_size
  env$simulation_seed <- simulation_seed
  env$run_simulations <- FALSE
  env$sim_grid_n_chunks <- 1L
  source(file.path(root_dir, "scripts", "r", "sim-compare-qmd7-presentation.R"), local = env)
  setup <- chunk("setup")
  assigned <- function(expr, name) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name(name))
  }
  semantics <- Filter(function(expr) assigned(expr, "comparison_semantics_version"), as.list(setup))
  if (length(semantics) != 1L) stop("Missing Analysis 7 comparison semantics.")
  eval(semantics[[1]], env)
  eval(chunk("actual-settings"), env)
  invisible(utils::capture.output(eval(chunk("compare-grid"), env)))
  params <- Filter(
    function(expr) assigned(expr, "analysis_result_params"),
    as.list(chunk("compare-parallel"))
  )
  if (length(params) != 1L) stop("Missing Analysis 7 result parameters.")
  eval(params[[1]], env)
  env$analysis_result_params
}
