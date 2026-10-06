# Monte Carlo uncertainty of simulation summaries, from the replicates that
# were already simulated (dataset bootstrap resamples saved results only).
# Source after analysis-plot-style.R. The functions return tibbles or plain
# vectors; the error-bar layer returns a ggplot2 layer and writes no files.
#
# Conventions. Wide summaries use the columns `<stat>` (the plotted value),
# `<stat>_lower`, `<stat>_upper` (a 95% interval) and `<stat>_mcse` (the Monte
# Carlo standard error, or for percentiles the half-width-based equivalent).
#
# - Means and proportions: MCSE = sd(x) / sqrt(n), with the usual n - 1 sd;
#   for a proportion of binary outcomes, sqrt(p (1 - p) / n). Interval:
#   estimate +/- 1.96 MCSE.
# - Percentiles (median, 90th, 95th, 10th, ...): a distribution-free interval
#   from order statistics (see `.analysis_mcse_quantile_index()`), reported
#   with an MCSE-equivalent of (upper - lower) / (2 x 1.96).
# - Maxima and minima: no interval (NA); the samples say little about how much
#   larger the maximum could be.
# - Equal-weight averages of k independent scenarios: MCSE =
#   sqrt(sum(MCSE_i^2)) / k (`.analysis_mcse_average()`).
# - Between-unit MCSE (`.analysis_mcse_between_units()`): for a mean of
#   dataset statistics (including paired method differences), MCSE is
#   sd(dataset statistics) / sqrt(D); NA when D < 5. This is distinct from
#   the uncertainty of a pooled sample percentile.
# - Pooled sample distributions in analyses 7/8: 999 deterministic whole-
#   dataset bootstrap draws recompute the plotted statistic, preserving CRN
#   pairing across methods/mismatches and joint scenario averages. Intervals
#   are bootstrap percentiles; MCSE is bootstrap SD. Require five contributing
#   datasets and at least 95% finite draws. Conditional dataset maxima additionally
#   require five affected eligible datasets; occurrence requires five eligible.
# Everything assumes the units (samples or datasets) are independent draws.

.analysis_mcse_z <- stats::qnorm(0.975)

# Fewest values for which an interval is reported.
.analysis_mcse_min_n <- 5L

.analysis_mcse_finite <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.finite(x)]
}

# Mean of the finite values with MCSE sd(x) / sqrt(n); NA MCSE when n < 2.
.analysis_mcse_mean <- function(x) {
  x <- .analysis_mcse_finite(x)
  n <- length(x)
  est <- if (n > 0L) mean(x) else NA_real_
  mcse <- if (n >= 2L) stats::sd(x) / sqrt(n) else NA_real_
  tibble::tibble(
    estimate = est,
    mcse = mcse,
    lower = est - .analysis_mcse_z * mcse,
    upper = est + .analysis_mcse_z * mcse,
    n = n
  )
}

# Proportion of TRUE among the non-missing values of a logical vector, with
# MCSE sqrt(p (1 - p) / n); NA MCSE when n < 2.
.analysis_mcse_proportion <- function(x) {
  x <- as.logical(x)
  x <- x[!is.na(x)]
  n <- length(x)
  p <- if (n > 0L) mean(x) else NA_real_
  mcse <- if (n >= 2L) sqrt(p * (1 - p) / n) else NA_real_
  tibble::tibble(
    estimate = p,
    mcse = mcse,
    lower = max(0, p - .analysis_mcse_z * mcse),
    upper = min(1, p + .analysis_mcse_z * mcse),
    n = n
  )
}

# Order-statistic indices of a distribution-free interval for the p-th
# quantile of n independent values: (x_(l), x_(u)) with
#   l = qbinom(alpha / 2, n, p) and u = qbinom(1 - alpha / 2, n, p) + 1.
# The number B of values at or below the true quantile is Binomial(n, p);
# x_(l) <= q exactly when B >= l, and x_(u) >= q when B <= u - 1, so the
# coverage is F(u - 1) - F(l - 1) >= (1 - alpha / 2) - alpha / 2 = 1 - alpha
# (F is the binomial distribution function). No randomisation, so every render
# gives the same interval. Returns c(lower = NA, upper = NA) when n < 5 or the
# interval would need a value beyond the sample (l < 1 or u > n), rather than
# clipping, which would lose the coverage guarantee.
.analysis_mcse_quantile_index <- function(n, p, level = 0.95) {
  alpha <- 1 - level
  na <- c(lower = NA_integer_, upper = NA_integer_)
  if (is.na(n) || n < .analysis_mcse_min_n) {
    return(na)
  }
  l <- as.integer(stats::qbinom(alpha / 2, n, p))
  u <- as.integer(stats::qbinom(1 - alpha / 2, n, p)) + 1L
  if (l < 1L || u > n) {
    return(na)
  }
  c(lower = l, upper = u)
}

# p-th quantile of the finite values (type 7, as `stats::quantile()`), with
# its order-statistic interval and MCSE-equivalent half-width.
.analysis_mcse_quantile <- function(x, p, level = 0.95) {
  x <- sort(.analysis_mcse_finite(x))
  n <- length(x)
  est <- if (n > 0L) {
    unname(stats::quantile(x, probs = p, names = FALSE))
  } else {
    NA_real_
  }
  idx <- .analysis_mcse_quantile_index(n, p, level)
  lower <- if (is.na(idx[["lower"]])) NA_real_ else x[[idx[["lower"]]]]
  upper <- if (is.na(idx[["upper"]])) NA_real_ else x[[idx[["upper"]]]]
  tibble::tibble(
    estimate = est,
    mcse = (upper - lower) / (2 * stats::qnorm(1 - (1 - level) / 2)),
    lower = lower,
    upper = upper,
    n = n
  )
}

# `.analysis_mcse_quantile()` bounds as one row of columns `<name>_lower`,
# `<name>_upper` and `<name>_mcse`, for splicing into `dplyr::summarise()`.
.analysis_mcse_quantile_cols <- function(x, p, name) {
  q <- .analysis_mcse_quantile(x, p)
  out <- tibble::tibble(q$lower, q$upper, q$mcse)
  names(out) <- paste0(name, c("_lower", "_upper", "_mcse"))
  out
}

# Maximum of the finite values, with no interval.
.analysis_mcse_max <- function(x) {
  x <- .analysis_mcse_finite(x)
  tibble::tibble(
    estimate = if (length(x)) max(x) else NA_real_,
    mcse = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    n = length(x)
  )
}

# Equal-weight average of k independent scenario estimates, each with its own
# MCSE: estimate = mean, MCSE = sqrt(sum(mcse^2)) / k. Scenarios with a
# non-finite estimate are left out (as the plotted averages leave them out); the
# combined MCSE is NA if any remaining scenario has no MCSE.
.analysis_mcse_average <- function(estimate, mcse) {
  keep <- is.finite(estimate)
  estimate <- estimate[keep]
  mcse <- mcse[keep]
  k <- length(estimate)
  est <- if (k > 0L) mean(estimate) else NA_real_
  se <- if (k > 0L && !anyNA(mcse)) sqrt(sum(mcse^2)) / k else NA_real_
  tibble::tibble(
    estimate = est,
    mcse = se,
    lower = est - .analysis_mcse_z * se,
    upper = est + .analysis_mcse_z * se,
    k = k
  )
}

# MCSE of a statistic from its spread between independent units (simulated
# datasets): `stat` is applied to the values of each unit, and the MCSE is
# sd(unit statistics) / sqrt(D) over the D units with a finite statistic.
# NA when D < `min_units`.
.analysis_mcse_between_units <- function(
    x, unit, stat, min_units = .analysis_mcse_min_n) {
  ok <- !is.na(unit)
  if (!any(ok)) {
    return(NA_real_)
  }
  vals <- vapply(
    split(x[ok], unit[ok]),
    function(v) {
      out <- suppressWarnings(stat(v))
      if (length(out) != 1L) NA_real_ else as.numeric(out)
    },
    numeric(1)
  )
  vals <- vals[is.finite(vals)]
  if (length(vals) < min_units) {
    return(NA_real_)
  }
  stats::sd(vals) / sqrt(length(vals))
}

# Dataset bootstrap streams are shared by true common-random-number families.
# Numeric sim_seed is the production family; textual keys support fixtures.
.analysis_mcse_family_seed <- function(family) {
  key <- paste(as.character(family), collapse = "|")
  value <- 20261005
  for (byte in utf8ToInt(key)) value <- (value * 131 + byte) %% 2147483646
  as.integer(value + 1)
}

.analysis_mcse_bootstrap_indices <- function(unit, family, reps = 999L) {
  ids <- sort(unique(as.character(unit[!is.na(unit)])))
  if (!length(ids)) return(matrix(integer(), 0L, reps))
  if (length(reps) != 1L || !is.finite(reps) || reps < 2L || reps != as.integer(reps)) {
    stop("Dataset bootstrap needs at least two integer replicates.")
  }
  withr::with_seed(.analysis_mcse_family_seed(family),
    matrix(sample.int(length(ids), length(ids) * reps, replace = TRUE), nrow = length(ids)),
    .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
}

.analysis_mcse_block_draws <- function(x, unit, stat, family, reps = 999L) {
  if (length(x) != length(unit)) stop("Dataset bootstrap needs one unit ID per value.")
  blocks <- split(x[!is.na(unit)], as.character(unit[!is.na(unit)]))
  indices <- .analysis_mcse_bootstrap_indices(unit, family, reps)
  if (!length(blocks)) return(rep(NA_real_, reps))
  vapply(seq_len(reps), function(i) {
    value <- suppressWarnings(stat(unlist(blocks[indices[, i]], use.names = FALSE)))
    if (length(value) != 1L) NA_real_ else as.numeric(value)
  }, numeric(1))
}

# Bootstrap coverage is explicit. Missing conditional statistics are never zero.
.analysis_mcse_bootstrap_cols <- function(
    estimate, draws, name, n_dataset, n_contributing,
    mcse = TRUE, min_contributing = 5L, coverage_min = 0.95) {
  finite <- draws[is.finite(draws)]
  coverage <- if (length(draws)) length(finite) / length(draws) else NA_real_
  available <- isTRUE(mcse) && is.finite(estimate) && n_dataset >= 5L &&
    n_contributing >= min_contributing && is.finite(coverage) && coverage >= coverage_min
  bounds <- if (available) stats::quantile(finite, c(0.025, 0.975), names = FALSE) else c(NA_real_, NA_real_)
  out <- tibble::tibble(value = estimate)
  names(out) <- name
  out[[paste0(name, "_mcse")]] <- if (available) stats::sd(finite) else NA_real_
  out[[paste0(name, "_lower")]] <- bounds[1]
  out[[paste0(name, "_upper")]] <- bounds[2]
  out[[paste0("n_bootstrap_", name)]] <- length(draws)
  out[[paste0("n_bootstrap_finite_", name)]] <- if (isTRUE(mcse)) length(finite) else NA_integer_
  out[[paste0("bootstrap_coverage_", name)]] <- if (isTRUE(mcse)) coverage else NA_real_
  out[[paste0("interval_available_", name)]] <- available
  out[[paste0("n_dataset_", name)]] <- n_contributing
  out[[paste0(".boot_", name)]] <- list(draws)
  out
}

.analysis_mcse_pooled_cols <- function(x, unit, stat, name, family, mcse = TRUE) {
  blocks <- split(x[!is.na(unit)], as.character(unit[!is.na(unit)]))
  contributing <- sum(vapply(blocks, function(v) is.finite(suppressWarnings(stat(v))), logical(1)))
  estimate <- suppressWarnings(stat(x[!is.na(unit)]))
  draws <- if (isTRUE(mcse)) .analysis_mcse_block_draws(x, unit, stat, family) else numeric()
  .analysis_mcse_bootstrap_cols(estimate, draws, name, length(blocks), contributing, mcse)
}

# Average a fixed set of originally finite scenario statistics in every draw.
# Vectors are aligned across CRN families; independent families have independent
# streams. An undefined constituent makes the whole draw undefined.
.analysis_mcse_bootstrap_average <- function(tbl, group_cols, stats, mcse = TRUE) {
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(rows, keys) {
      result <- tibble::tibble(n_scenario = nrow(rows))
      for (name in stats) {
        keep <- is.finite(rows[[name]])
        estimate <- if (any(keep)) mean(rows[[name]][keep]) else NA_real_
        vectors <- rows[[paste0(".boot_", name)]][keep]
        draws <- if (isTRUE(mcse) && length(vectors) && all(lengths(vectors) == 999L)) {
          matrix <- do.call(cbind, vectors)
          apply(matrix, 1L, function(v) if (all(is.finite(v))) mean(v) else NA_real_)
        } else numeric()
        contributing <- if (any(keep)) min(rows[[paste0("n_dataset_", name)]][keep]) else 0L
        cols <- .analysis_mcse_bootstrap_cols(estimate, draws, name, contributing, contributing, mcse)
        result <- dplyr::bind_cols(result, cols)
        result[[paste0("n_scenario_", name)]] <- sum(keep)
      }
      result
    }) |>
    dplyr::ungroup()
}

# Quantile of the finite values, NA when there are none.
.analysis_mcse_quantile_finite <- function(x, p) {
  x <- .analysis_mcse_finite(x)
  if (length(x) == 0L) {
    return(NA_real_)
  }
  unname(stats::quantile(x, probs = p, names = FALSE))
}

# Add `<stat>_lower` and `<stat>_upper` as `<stat> +/- 1.96 <stat>_mcse`,
# clipped to `range` (e.g. c(0, 1) for proportions, c(0, Inf) for sizes).
.analysis_mcse_add_bounds <- function(tbl, stats, range = c(-Inf, Inf)) {
  for (s in stats) {
    se <- tbl[[paste0(s, "_mcse")]]
    tbl[[paste0(s, "_lower")]] <- pmax(range[[1]], tbl[[s]] - .analysis_mcse_z * se)
    tbl[[paste0(s, "_upper")]] <- pmin(range[[2]], tbl[[s]] + .analysis_mcse_z * se)
  }
  tbl
}

# Average the `stats` columns equally over the rows of each group, combining
# their `<stat>_mcse` columns with `.analysis_mcse_average()`, and add
# `<stat>_lower` / `<stat>_upper`. Statistics without an `_mcse` column are
# averaged with NA bounds.
.analysis_mcse_average_cols <- function(tbl, group_cols, stats) {
  for (s in stats) {
    if (!paste0(s, "_mcse") %in% names(tbl)) {
      tbl[[paste0(s, "_mcse")]] <- NA_real_
    }
  }
  grouped <- dplyr::group_by(tbl, dplyr::across(dplyr::all_of(group_cols)))
  out <- dplyr::summarise(grouped, n_scenario = dplyr::n(), .groups = "drop")
  for (s in stats) {
    avg <- dplyr::summarise(
      grouped,
      .avg = list(.analysis_mcse_average(
        .data[[s]], .data[[paste0(s, "_mcse")]]
      )),
      .groups = "drop"
    )$.avg
    avg <- dplyr::bind_rows(
      tibble::tibble(
        estimate = numeric(), mcse = numeric(), lower = numeric(),
        upper = numeric(), k = integer()
      ),
      avg
    )
    out[[paste0("n_scenario_", s)]] <- avg$k
    out[[s]] <- avg$estimate
    out[[paste0(s, "_mcse")]] <- avg$mcse
    out[[paste0(s, "_lower")]] <- avg$lower
    out[[paste0(s, "_upper")]] <- avg$upper
  }
  out
}

# Error bars for the finite `ymin`/`ymax` columns of `data`, in the colour of
# the line they belong to (the colour aesthetic is inherited). Small caps
# and opaque strokes keep intervals visible even with a single x value.
.analysis_mcse_errorbar <- function(
    data, ymin = "lower", ymax = "upper", width = 0.1, alpha = 1,
    linewidth = 0.4) {
  data <- data[
    is.finite(data[[ymin]]) & is.finite(data[[ymax]]), ,
    drop = FALSE
  ]
  .analysis_mcse_layer(ggplot2::geom_errorbar(
    data = data,
    mapping = ggplot2::aes(ymin = .data[[ymin]], ymax = .data[[ymax]]),
    width = width,
    alpha = alpha,
    linewidth = linewidth
  ))
}
