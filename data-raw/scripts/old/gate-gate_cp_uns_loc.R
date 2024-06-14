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

  prob_mod <- try(scam::scam(prob_stim_norm ~ s(x_stim, bs = "mpi"), #
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

#' @title Plot gating plots for local fdr method
#'
#' @description
#' Plot the following three plots:
#' \describe{
#'   \item{p_loc_dens}{Density of the stim and unstim (pos and neg) distributions.}
#'   \item{p_loc_prob}{Raw, normalised and smoothed probability of positivty by marker value.}
#'   \item{p_loc_ctb}{Contribution to total count by binned marker values}.
#' }
#'
#' @param dens_tbl dataframe. Dataframe equal to \code{dens_tbl_raw}
#' in the \code{.get_cp_uns_loc_ind} function.
#' @inheritParams plot_cp
#' @param prob_mod object of class \code{glm}. Models probability of positivity
#' as a function of \code{x_stim}. Created within \code{.get_cp_uns_loc_ind}
#' function.
#' @param prob_tbl dataframe. Dataframe contains columns x_stim,
#' prob_stim, prob_stim_norm and pred, where x_stim are marker values
#' and prob_stim, prob_stim_norm and pred are the raw, normalised and
#' smoothed probabilities of posititivity, respectively. Dataframe is
#' produced appropriately within \code{.get_cp_uns_loc_ind} function.
#' @param cut_stim numeric vector. Each element is an observed value of the marker
#' in the stim condition.
#' @param min_x_pos_prob numeric. Minimum value for \code{cut_stim} such that
#' the probability of positivity is greater than the value specified by the
#' \code{prob_min} parameter in the \code{.get_cp_uns_loc_ind} function.
#'
#' @return A list with elements named p_loc_dens, p_loc_prob and
#' p_loc_ctb, where each element is a ggplot2 plot.
#'
#' @examples
#' For example of use, set \code{debugonce(:::.get_cp_uns_loc_ind)}
#' before running \code{::gate}. Step through the debuggec function
#' until near the end.
.plot_cp_loc_fdr <- function(dens_tbl,
                             params,
                             prob_tbl,
                             pred_tbl,
                             cut_stim,
                             sample,
                             min_x_pos_prob,
                             path_project){

  dir_base <- stimgate_dir_base_create( # nolint
    params = params, dir_base_init = path_project
  )
  dir_base <- file.path(dir_base, "gating_plots", "gate")
  dir_len <- stringr::str_length(dir_base) + 10
  underscore_loc <- stringr::str_locate(sample, "_")[1,"start"][[1]]
  subject_visittype <- stringr::str_sub(sample, end = underscore_loc - 1)
  stim <- stringr::str_sub(sample, underscore_loc + 1)
  subjectid_visittype_len <- 255 - dir_len - 4 - 4
  subjectid_visittype_len <- min(
    stringr::str_length(subject_visittype), subjectid_visittype_len
  )
  subject_visittype <- stringr::str_sub(
    subject_visittype,
    end = subjectid_visittype_len
  )
  fn <- paste0(subject_visittype, "-", stim, ".png")

  # Plot densities
  # -------------------------------
  p_loc_dens <- ggplot(dens_tbl |>
                         dplyr::filter(x_stim > min(min_x_pos_prob, 100)),
                       aes(x = x_stim, y = dens, linetype = stim)) +
    cowplot::theme_cowplot(font_size = 20) +
    geom_line() +
    scale_linetype_manual(values = c("yes" = "solid",
                                     "no" = "dotted"),
                          labels = c("yes" = "stim",
                                     "no" = "unstim")) +
    labs(x = params$chnl_lab[params$cut], y = "Density") +
    theme(legend.title = element_blank())

  dir_save <- file.path(dir_base, "p_loc_dens")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  cowplot::ggsave2(filename = file.path(file.path(dir_save, fn)),
                   height = 10, width = 10)

  # Plot probabilities of being stim-induced
  # -------------------------------

  plot_tbl_prob <- prob_tbl |>
    tidyr::pivot_longer(prob_stim:prob_stim_norm,
                 names_to = "prob_type",
                 values_to = "prob") |>
    dplyr::select(x_stim, prob_type, prob)

  plot_tbl_pred <- pred_tbl |>
    dplyr::rename(x_stim = cut_stim,
           prob = pred) |>
    dplyr::mutate(prob_type = "pred") |>
    dplyr::select(x_stim, prob_type, prob )

  plot_tbl <- plot_tbl_prob |>
    dplyr::bind_rows(plot_tbl_pred)
  p_loc_prob <- ggplot(plot_tbl,
                       aes(x = x_stim, y = prob, linetype = prob_type)) +
    cowplot::theme_cowplot(font_size = 20) +
    cowplot::background_grid(major = "y", minor = "y") +
    geom_hline(yintercept = 0) +
    geom_line(size = 2) +
    scale_linetype_manual(values = c("pred" = "solid",
                                     "prob_stim_norm" = "dotted",
                                     "prob_stim" = "dashed"),
                          labels = c("pred" = "Smoothed - Model",
                                     "prob_stim" = "Raw",
                                     "prob_stim_norm" = "Smoothed - Raw")) +
    lims(y = c(-1, 1)) +
    theme(legend.title = element_blank()) +
    labs(x = params$chnl_lab[params$cut],
         y = "Probability of being stimulation-induced")

  dir_save <- file.path(dir_base, "p_loc_prob")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  cowplot::ggsave2(filename = file.path(file.path(dir_save, fn)),
                   height = 10, width = 10)

  # Plot contribution to count by bin
  # -------------------------------

  # plot binned contributions to total count by range

  #  n_bin <- hist(new_pred_tbl$x_stim, plot = FALSE)$mids |> length
  n_bin <- hist(pred_tbl$cut_stim, plot = FALSE)$mids |> length()
  if(n_bin == 1){
    bin_tbl <- pred_tbl |>
      dplyr::mutate(bin = paste0("(", cut_stim, ",", cut_stim, "]"),
             ctb = pred)
  } else {
    bin_tbl <- pred_tbl |>
      dplyr::mutate(bin = cut(.data$cut_stim, breaks = n_bin)) |>
      dplyr::group_by(bin) |>
      dplyr::summarise(ctb = sum(pred))
  }

  #' @title Get the midpoint of a bin returned by base::cut
  .get_bin_mid <- function(bin) {
    purrr::map_dbl(bin, function(bin_curr){
      comma_loc <- stringr::str_locate(bin_curr, ",")[1, "start"][[1]]
      num_1 <- stringr::str_sub(bin_curr, 2, comma_loc - 1) |>
        as.numeric()
      num_2 <- stringr::str_sub(
        bin_curr, comma_loc + 1, stringr::str_length(bin_curr) - 1
        ) |>
        as.numeric()
      mean(c(num_1, num_2))
    })
  }

  bin_tbl <- bin_tbl |>
    dplyr::mutate(bin_mid = .get_bin_mid(bin))

  # plot of contribution per bin
  p_loc_ctb <- ggplot(bin_tbl,
                      aes(x = bin_mid, y = ctb)) +
    geom_bar(stat = "identity",
             col = "gray25",
             fill = "gray85")

  dir_save <- file.path(dir_base, "p_loc_ctb")
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  cowplot::ggsave2(filename = file.path(file.path(dir_save, fn)),
                   height = 10, width = 10)

  #list(p_loc_dens = p_loc_dens, p_loc_prob = p_loc_prob,
  #     p_loc_ctb = p_loc_ctb)
  list()

}


