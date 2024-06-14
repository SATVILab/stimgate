# extract and save plots
# ---------------------------------------

# save to temp directory
if (plot && FALSE) {
  # add bias label
  cp_uns_plot_list <- cp_uns_gate_combn_obj[["p_list"]]
  for (i in seq_along(cp_uns_plot_list)) {
    for (j in seq_along(cp_uns_plot_list[[i]])) {
      names(cp_uns_plot_list[[i]][[j]]) <- stringr::str_replace(
        names(cp_uns_plot_list[[i]][[j]]),
        "loc",
        paste0("locb", bias)
      )
    }
  }

  dir_save <- file.path(
    tempdir(),
    params$data_name,
    paste0("cp_locb", bias, "_plots")
  )
  if (!dir.exists(dir_save)) {
    dir.create(dir_save, recursive = TRUE)
  }
  saveRDS(
    cp_uns_plot_list,
    file.path(dir_save, paste0(
      names(cp_uns_plot_list),
      ".rds"
    ))
  )
}


# okay, so this all ties into the `prob_smooth`
# values that we actually end up modelling.
# presumably there were some issues here.
# okay, I think we're just smoothing values
# just the mean of all values within a given range around it, I think.
# not clear how important this is.
# could add a parameter `smooth` for this.
# can also add a parameter `monotonic`
# for the later smoothing.
# I don't know how much this helps, but it's
# included for now!
# I'm sure we could get a faster smoother, this seems rather
# slow as it's a loop.
# I think this probably helps prevent odd things happening.

# this really should get kicked out, can't we just average within
# narrow bins?
# surely it doesn't change that much anyway?
# don't smooth

if (FALSE) {
  break_vec <- seq(min(data_mod$expr), max(data_mod_expr$expr), by = margin)
  break_vec <- c(break_vec, max(data_mod$expr)) |> unique()
  prob_smooth_vec <- data_mod$expr
  for (i in seq_len(length(break_vec) - 1)) {
    if (i == length(break_vec) - 1) {
      ind_vec <- data_mod$expr >= break_vec[i] &
        data_mod$expr <= break_vec[i + 1]
    } else {
      ind_vec <- data_mod$expr >= break_vec[i] &
        data_mod$expr < break_vec[i + 1]
    }
    if (sum(ind_vec) == 0L) {
      next
    }
    val_vec <- data_mod$expr[ind_vec]
    prob_smooth_vec[ind_vec] <- mean(val_vec)
  }
  data_mod[, "prob_smooth"] <- prob_smooth_vec
}
if (FALSE) {
  data_mod <- data_mod |>
    dplyr::mutate(
      prob_smooth = purrr::map_dbl(
        seq_along(expr),
        function(i) {
          x <- data_mod$expr[i]

          # swm: smaller within margin
          # lwm: larger within margin
          prob_tbl_pos_margin |>
            dplyr::mutate(
              swm = x_stim > x - margin & x_stim <= x,
              lwm = x_stim < x + margin & x_stim >= x
            ) |>
            dplyr::filter(swm | lwm) |>
            dplyr::pull("prob_stim_norm") |>
            mean()
        }
      )
    )
}

if (FALSE) {
  data_plot <- data_mod |>
    dplyr::select(expr, prob_smooth, pred) |>
    dplyr::rename(x_stim = expr) |>
    tidyr::pivot_longer(
      names_to = "type",
      values_to = "prob",
      -x_stim
    ) |>
    dplyr::bind_rows(prob_tbl_pos |>
      dplyr::select(x_stim, prob_stim_norm) |>
      dplyr::mutate(type = "prob_stim_norm") |>
      dplyr::rename(prob = prob_stim_norm))

  ggplot(data_plot, aes(x = x_stim, col = type)) +
    geom_line(aes(y = prob)) +
    geom_rug(data = data_plot |>
      dplyr::filter(type == "pred")) +
    scale_colour_brewer(palette = "Accent")
}

if (nrow(data_mod) >= 10) {
  .debug(debug, "Smoothing I")

  # monotonic increasing
  fit <- try(
    scam::scam(
      prob_smooth ~ s(expr, bs = "mpi"),
      family = "binomial",
      data = data_mod |>
        dplyr::mutate(
          prob_smooth = pmin(prob_smooth, 0.999),
          prob_smooth = pmax(prob_smooth, 0.001)
        ),
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        devtol.fit = 0.5,
        steptol.fit = 1e-1,
        maxHalf = 5,
        bfgs = list(steptol.bfgs = 1e-1),
        maxit = 1e1
      )
    ),
    silent = TRUE
  )
  mean_abs_error <- 0
  pred_vec <- 0

  if (!inherits(fit, "try-error")) {

  }

  # monotonic increasing and convex
  if (
    inherits(fit, "try-error") ||
      all(pred_vec > 0.99) ||
      mean_abs_error > 0.3
  ) {
    .debug(debug, "Smoothing II")
    fit <- try(
      scam::scam(
        prob_smooth ~ s(expr, bs = "micv"),
        family = "binomial",
        data = data_mod,
        control = scam::scam.control(
          print.warn = FALSE,
          trace = FALSE,
          devtol.fit = 0.01 # ,
          # steptol.fit = 1e-3,
          # bfgs = list(steptol.bfgs = 1e-4),
          # maxit = 5
        )
      ),
      silent = TRUE
    )

    if (!inherits(fit, "try-error")) {
      pred_vec <- predict(fit, type = "response")
      mean_abs_error <- mean(abs(pred_vec - data_mod$prob_smooth))
    }
  }

  # smoothes
  if (
    inherits(fit, "try-error") ||
      all(pred_vec > 0.99) ||
      mean_abs_error > 0.3) {
    # .debug(debug, "Smoothing III")
    .debug(debug, "Skipping mgcv smoothing")
    # It's very slow, and returned
    # a WAAAY larger error than scam
    # when examined.
    # can try isotonic regression here instead.
    # (though perhaps could do a moving average first...)
    # fit <- try(
    #   mgcv::gam(prob_smooth ~ s(expr),
    #     family = "binomial",
    #     data = data_mod
    #   ),
    #   silent = TRUE
    # )
    if (!inherits(fit, "try-error")) {
      pred_vec <- predict(fit, type = "response")
      mean_abs_error <- mean(abs(pred_vec - data_mod$prob_smooth))
    }
  }
  if (
    inherits(fit, "try-error") ||
      all(pred_vec > 0.99) ||
      mean_abs_error > 0.3) {
    .debug(debug, "Failed to smooth")
    data_mod <- data_mod |>
      dplyr::mutate(pred = prob_smooth - 0.0001)
  } else {
    .debug(debug, "Smoothed")
    data_mod <- data_mod |>
      dplyr::mutate(pred = pred_vec)
  }
} else {
  .debug(debug, "Failed to smooth")
  data_mod <- data_mod |>
    dplyr::mutate(pred = prob_smooth - 1e-4)
}

if (FALSE) {
  ggplot(
    data_threshold |>
      tidyr::pivot_longer(
        names_to = "type",
        values_to = "prop",
        c(prop_stim, prop_uns, prop_bs)
      ),
    aes(x = expr, y = prop * 1e2, col = type)
  ) +
    geom_line()
}

if (plot && FALSE) {
  p_list <- .plot_cp_loc_fdr(
    dens_tbl = dens_tbl_raw,
    params = params,
    prob_tbl = prob_tbl,
    pred_tbl = all_cell_pred_tbl,
    cut_stim = orig_list$stim$expr,
    sample = orig_list$stim$sample[1],
    min_x_pos_prob = min(all_cell_pred_tbl$cut_stim),
    path_project = path_project
  )
} else {
  p_list <- .get_cp_uns_loc_p_list_empty()
}

# print('done getting loc fdr gate for this single sample')

# min_diff_tbl
# min_diff_tbl$prop_bs * 1e2

# check that marginal gain in signal-to-noise ratio is positive
# if not, move back up until a maximum of margin units
# data_threshold |>
#  dplyr::mutate(prop_bs_sd_std = sqrt(prop_bs_sd * pmin(prop_bs_sd, min_diff_bs_sd))) |>
#  dplyr::mutate(prop_bs_diff_sd_std = abs(prop_bs_diff)/prop_bs_sd_std) |>
#  dplyr::select(expr, prop_stim, prop_uns, prop_bs, prop_bs_diff, prop_bs_diff_sd_std) |>
#  dplyr::filter(prop_bs_diff_sd_std <= min_diff_sd + 0.75 * min_diff_bs_sd) |>
#  dplyr::arrange(expr) |>
#  dplyr::slice(1)

# explore s2n
# data_threshold |>
#  dplyr::mutate(s2n = prop_bs/prop_bs_sd) |>
#  dplyr::mutate(s2n_marginal = diff(c(0, s2n))) |>
#  dplyr::select(expr, prop_bs, prop_bs_sd, s2n, s2n_marginal)

.get_cell_specific_response_probs <- function(prob_tbl, nrow_level, params,
                                              prob_min, cut_stim) {
  # =====================
  # Fit model
  # =====================

  # fit a model to get smoothed probabilities
  prob_mod <- .get_prob_mod(prob_tbl, nrow_level = nrow_level)

  # if model fails, then return it
  if ("ultimate error" %in% class(prob_mod)) {
    out <- list(
      pred_tbl = NA,
      prob_mod = NA
    )

    # out <- FALSE
    # class(out) <- 'ultimate error'
    return(out)
  }

  # =====================
  # Get predictions
  # =====================

  # get fitted values
  pred_vec <- fitted.values(prob_mod)

  # add them to table
  pred_tbl <- prob_tbl |>
    dplyr::mutate(pred = pred_vec)

  out <- list(
    pred_tbl = pred_tbl,
    prob_mod = prob_mod
  )

  out
}

.get_prob_mod <- function(prob_tbl, nrow_level){

  prob_mod <- try(scam::scam(prob_stim_norm ~ s(x_stim, bs = "mpi"),
                         family = "binomial",
                         data = prob_tbl |>
                           dplyr::filter(prob_stim_norm > 0.25)))

  if(inherits(prob_mod, "try-error")){
    prob_mod <- try(scam::scam(prob_stim_norm ~ s(x_stim, bs = "micx"),
                               family = "binomial",
                               data = prob_tbl))
  } else return(prob_mod)

  if(inherits(prob_mod, "try-error")){
    prob_mod <- try(mgcv::gam(prob_stim_norm ~ s(x_stim),
                               family = "binomial",
                               data = prob_tbl))
  }

  prob_mod
}