#' @title Calculate the local fdr-based cut
.get_cp_uns_loc <- function(ex_list,
                            ind_gate,
                            ind_uns,
                            gate_combn,
                            pop_root = NULL,
                            data,
                            bias_uns = 0,
                            noise_sd = NULL,
                            min_bw = 80,
                            cp_min,
                            min_cell,
                            params,
                            plot,
                            path_project,
                            debug = FALSE){

  # get cutpoints for each level of bias
  .get_cp_uns_loc_bias(
    ex_list = ex_list,
    ind_gate = ind_gate,
    ind_uns = ind_uns,
    data = data,
    bias_uns = bias_uns,
    noise_sd = NULL,
    cp_min = cp_min,
    gate_combn = gate_combn,
    min_bw = min_bw,
    min_cell = min_cell,
    params = params,
    plot = plot,
    path_project = path_project,
    debug = debug
  )

}

#' @title Get the unstim-based local fdr-method cutpoint for each level of bias
.get_cp_uns_loc_bias <- function(ex_list,
                                 ind_gate, ind_uns,
                                 data,
                                 bias_uns,
                                 noise_sd,
                                 cp_min,
                                 gate_combn,
                                 min_bw,
                                 min_cell,
                                 params,
                                 plot,
                                 path_project,
                                 debug){

  # get ecdf of uns
  purrr::map(bias_uns, function(bias){
    .debug(debug, "bias_uns", bias)
    #print('getting loc fdr gate for a given bias')

    # get ecdf of uns sample
    # -------------------------------------
    cut_tbl_uns <- .get_cut_list(
      ex_list = ex_list,
      ind = ind_uns,
      exc_min = TRUE,
      bias = bias,
      debug = debug,
      noise_sd = NULL)[[1]]

    # get gates for given level of bias across gate combination methods
    # --------------------------------------
    cp_uns_gate_combn_obj <- .get_cp_uns_loc_gate_combn(
      ex_list =  ex_list,
      cut_uns = cut_tbl_uns,
      ind_uns = ind_uns,
      ind_gate = ind_gate,
      exc_min = TRUE,
      gate_combn = gate_combn,
      cp_min = cp_min,
      min_bw = min_bw,
      min_cell = min_cell,
      params = params,
      plot = plot,
      bias = bias,
      path_project = path_project,
      debug = debug
    )
    # extract and save plots
    # ---------------------------------------

    # save to temp directory
    if (plot && FALSE){
      # add bias label
      cp_uns_plot_list <- cp_uns_gate_combn_obj[['p_list']]
      for(i in seq_along(cp_uns_plot_list)){
        for(j in seq_along(cp_uns_plot_list[[i]])){
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
      if(!dir.exists(dir_save)) {
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

    # extract and add bias label to gates
    # ---------------------------------------
    cp_uns_gate_combn_list <- cp_uns_gate_combn_obj[['cp_uns']]
    names(cp_uns_gate_combn_list) <-
      paste0(names(cp_uns_gate_combn_list), "b", bias)
    cp_uns_gate_combn_list

    #print('done getting loc fd
    # gate across gate combinations for a given bias')

    cp_uns_gate_combn_list

  }) |>
    purrr::flatten()
}

#' @title Get the unstim-based local fdr-method cutpoint for a given bias across gate combination methods
.get_cp_uns_loc_gate_combn <- function(ex_list,
                                       cut_uns,
                                       ind_uns,
                                       ind_gate,
                                       exc_min = TRUE,
                                       gate_combn,
                                       cp_min,
                                       min_bw,
                                       min_cell,
                                       params,
                                       plot,
                                       bias,
                                       path_project,
                                       debug = FALSE){
  .debug(debug, "getting gate_combn")

  #print('getting loc fdr gate across gate combinations')

  # get cutpoints for prejoin gate combination method
  if ('prejoin' %in% gate_combn) {
    .debug(debug, "prejoin")

    # get marker expression for stim samples,
    # join and then sort into descending order
    cut_list_stim <- .get_cut_list(
      ex_list = ex_list,
      ind = setdiff(ind_gate, ind_uns),
      exc_min = TRUE,
      bias = 0,
      noise_sd = NULL,
      debug = debug
    ) |>
      unlist() |>
        sort() |>
        rev() |>
        list()

    # get cutpoints for gate combn method for a range of fdr's
    cp_uns_list_prejoin <- .get_cp_uns_loc_sample(
      cut_stim = cut_list_stim,
      cut_uns = cut_uns,
      cp_min = cp_min,
      ind_uns = ind_uns,
      ind_gate = ind_gate,
      min_bw = min_bw,
      min_cell = min_cell,
      params = params,
      plot = plot,
      bias = bias,
      path_project = path_project
    ) |>
      purrr::map(function(x) list('prejoin' = x))


  } else cp_uns_list_prejoin <- list()

  # get cutpoint if group method is not only prejoin
  non_prejoin_combn_vec <- setdiff(gate_combn, 'prejoin')

  if (length(non_prejoin_combn_vec) > 0){
    .debug(debug, "non-prejoin")

    # get marker expression for stim samples, and sort into descending order
    cut_list_stim <- .get_cut_list(
      ex_list = ex_list,
      ind = setdiff(ind_gate, ind_uns),
      exc_min = TRUE,
      bias = 0,
      noise_sd = NULL
      ) |>
      purrr::map(function(x) x |> dplyr::arrange(desc(expr)))

    cp_uns_list_nonjoin <- .get_cp_uns_loc_sample(
      cut_stim = cut_list_stim,
      cut_uns = cut_uns,
      cp_min = cp_min,
      ind_uns = ind_uns,
      ind_gate = ind_gate,
      min_bw = min_bw,
      min_cell = min_cell,
      params = params,
      plot = plot,
      bias = bias,
      path_project = path_project
    )

    # get list of plots organised
    # ---------------------------
    .debug(debug, "Organising plots")

    indices_name_vec <- paste0(
      names(cp_uns_list_nonjoin[["loc"]]),
      collapse = "_"
    )
    sample_name_vec <- purrr::map_chr(
      names(cp_uns_list_nonjoin[['p_list']]),
      function(ind){
        paste0(ex_list[[ind]]$batch_sh[1],
              "_", ex_list[[ind]]$stim[1])
    })

    p_list_sample_level <- stats::setNames(
      cp_uns_list_nonjoin[["p_list"]],
      sample_name_vec
    )

    p_list <- stats::setNames(
      list(p_list_sample_level), indices_name_vec
    )


    # get list of cutpoints combined in the appropriate way
    # ---------------------------

    .debug(debug, "Combining cutpoints")
    cp_uns_list_nonjoin <- .combine_cp(
      cp = cp_uns_list_nonjoin[["loc"]],
      gate_combn = non_prejoin_combn_vec
    )
  } else cp_uns_list_nonjoin <- list()

  .debug(debug, "Getting unstim cutpoints")

  combined_list <- cp_uns_list_prejoin |>
    append(cp_uns_list_nonjoin)
  cp_uns_list <- purrr::map(
    unique(names(combined_list)),
    function(x){
      .debug(debug, "cutpoint name", paste0(x, collapse = "-"))
      cp_uns_list_prejoin[[x]] |>
        append(cp_uns_list_nonjoin[[x]])
    }) |>
      stats::setNames(unique(names(combined_list)))

  #print('done getting loc fdr gates across gate combinations')
  .debug(debug, "done getting gate_combn")

  list(
    "cp_uns" = list("loc" = cp_uns_list),
    "p_list" = p_list
  )
}

#' @title Get cutpoint for a range of samples given the q-value and fdr
#'
#' @description Calculate the cutpoint for each sample in a batch at a given FDR.
#'
#' @param cut_stim list. List where each element are the marker expression readings
#' of the marker to be cut on for the cells in a sample. Note that the i-th element
#' in \code{cut_stim} must correspond to the i-th element in \code{q_list}, i.e.
#' must be related to the same marker in same cell population from the same blood sample and stimulation.
#' @param fdr numeric. A value between 0 and 1 specifying the false discovery rate
#' the sample should be cut at.
#'
#' @return Numeric vector. A cutpoint for each sample.
.get_cp_uns_loc_sample <- function(cut_stim,
                                   cut_uns,
                                   cp_min,
                                   min_bw,
                                   ind_uns,
                                   ind_gate,
                                   min_cell,
                                   params,
                                   plot = TRUE,
                                   bias,
                                   path_project,
                                   debug = FALSE){
  .debug(debug, "getting loc gate at sample level")
  #print('getting loc fdr gates across sample')

  # get cutpoints for each sample
  cp_uns_loc_obj_list <- purrr::map(seq_along(cut_stim), function(i) {
    .debug(debug, "sample", i)

    if (nrow(cut_stim[[i]]) < min_cell || nrow(cut_uns) < min_cell) {
      .debug(debug, "Too few cells")
      p_list <- lapply(1:3, function(x) ggplot()) |>
        stats::setNames(c("p_loc_dens", "p_loc_prob", "p_loc_ctb"))
      return(
        list(p = NA, p_list = p_list)
      )
    }

    # remove any cytokine-positive cells from unstim using gates from
    # sample for which single-positive gates are required
    if (!is.null(params$gate_tbl)) {
      .debug(debug, "Removing cytokine-positive cells from unstim")

      ex_uns <- params$ex_uns

      # first filter
      gate_tbl_gn_ind <- params$gate_tbl |>
        # filter to use gates from cut_stim
        dplyr::filter(
          ind == cut_stim[[i]]$ind[1],
          gate_name == params$gate_name_curr
        )

      pos_ind_vec_but_single_pos_curr <-
        .get_pos_ind_but_single_pos_for_one_cyt(
          ex = ex_uns,
          gate_tbl = gate_tbl_gn_ind,
          chnl_single_exc = params$cut,
          chnl = NULL,
          gate_type_cyt_pos = ifelse(
            params$calc_cyt_pos_gates,
            'cyt', 'base'
          ),
          gate_type_single_pos = 'base'
        )

      ex_uns <- ex_uns[
        !pos_ind_vec_but_single_pos_curr, , drop = FALSE
      ]

      # filter unstim based on this
      cut_uns <- .get_cut_list(
        ex_list = stats::setNames(list(ex_uns), ex_uns$ind[1]),
        ind = ex_uns$ind[1],
        exc_min = TRUE,
        bias = bias,
        noise_sd = NULL
      )[[1]]
    }

    cp_uns_loc_obj <- .get_cp_uns_loc_ind(
      cut_stim = cut_stim[[i]],
      cut_uns = cut_uns,
      min_bw = min_bw,
      cp_min = cp_min,
      min_cell = min_cell,
      params = params,
      plot = plot,
      bias = bias,
      path_project = path_project,
      debug = debug
    )
  })

  # name sample
  # ------------------

  .debug(debug, "Possibly re-using calculated cutpoints")

  # extract vector of cutpoints
  cp_vec <- purrr::map_dbl(
    cp_uns_loc_obj_list,
    function(cp_uns_loc_obj) {
      cp_uns_loc_obj$cp
    }
  )

  # stimulation indices
  ind_stim <- setdiff(ind_gate, ind_uns)

  cp_vec <- purrr::map_dbl(
    cp_uns_loc_obj_list,
    function(cp_uns_loc_obj) {
      cp_uns_loc_obj$cp
    }
  )

  # Repeat cutpoint if it was prejoined
  if (length(cp_vec) == 1 && ind_stim != 1) {
    cp_vec <- stats::setNames(
      rep(cp_vec, length(ind_stim)), ind_stim
    )
  } else {
    # name gate indices if not prejoined
    cp_vec <- stats::setNames(cp_vec, ind_stim) 
  } 

  # add unstim if unstim gate required
  if(ind_uns %in% ind_gate){
    if(!all(purrr::map_lgl(cp_vec, is.na))){
      cp_vec <- c(cp_vec, stats::setNames(mean(cp_vec, na.rm = TRUE), ind_uns))
    } else cp_vec <- c(cp_vec, NA)
  }

  p_list <- purrr::map(cp_uns_loc_obj_list, function(x) x$p_list) |>
    stats::setNames(ind_stim)

  #print('done getting loc fdr gates for this batch')
  .debug(debug, "done getting loc gate at sample level")

  # collate plots
  list("loc" = cp_vec,
       "p_list" = p_list)
}

.get_cp_uns_loc_ind <- function(cut_uns,
                                cut_stim,
                                min_bw,
                                cp_min,
                                min_cell,
                                params,
                                plot = TRUE,
                                prob_min = 0.1,
                                bias,
                                path_project,
                                debug = FALSE) {

  .debug(debug, "getting loc gate for single sample")

  p_list_empty <- lapply(1:3, function(x) ggplot()) |>
    stats::setNames(c("p_loc_dens", "p_loc_prob", "p_loc_ctb"))

  # estimate densities for stim and unstim over stim range
  if (nrow(cut_stim) < min_cell) {
    .debug(debug, "Too few cells")
    return(list(
      cp = max(
        cp_min,
        cut_stim$expr +
          (max(cut_stim$expr) - min(cut_stim$expr)) / 5
        ),
      p_list = p_list_empty
    ))
  }

  cut_stim_orig <- cut_stim; cut_uns_orig <- cut_uns
  max_dens_x <- max(cut_stim$expr) - 0.05 * (diff(range(cut_stim$expr)))

  #max_dens_x <- min(180, max_dens_x)
  #while(sum(cut_stim_orig > max_dens_x) < 3) max_dens_x <- max_dens_x - 5

  if (max_dens_x <= cp_min) {
    .debug(debug, "Max value less than minimum threshold")
    return(
      list(
        cp = max(
          cp_min,
          cut_stim$expr + (max(cut_uns$expr) - min(cut_uns$expr)) / 5
          ),
        p_list = p_list_empty
      )
    )
  }

  cut_stim <- cut_stim |>
    dplyr::mutate(
      expr = ifelse(
        .data$expr > max_dens_x, max_dens_x, .data$expr
      )
    )
  cut_uns <- cut_uns |>
    dplyr::mutate(
      expr = ifelse(
        .data$expr > max_dens_x, max_dens_x, .data$expr
      )
    )

  # min_bw_x <- max_dens_x
  # while(
  #   sum(cut_stim$expr > min_bw_x) <= 2 ||
  #     sum(cut_stim$expr[cut_stim$expr != max_dens_x] > min_bw_x) <= 2) {
  #   min_bw_x <- min_bw_x - 1
  # }

  #bw_stim <-  density(cut_stim$expr[cut_stim$expr > min_bw_x])$bw
  #bw_uns <-  density(cut_uns$expr[cut_uns$expr > min_bw_x])$bw
  .debug(debug, "Calculating densities")
  bw_stim <-  try(density(cut_stim$expr, bw  = "SJ")$bw)
  bw_stim <- ifelse(class(bw_stim) == 'try-error', min_bw, bw_stim)
  bw_uns <-  try(density(cut_uns$expr, bw = "SJ")$bw)
  bw_uns <- ifelse(class(bw_uns) == 'try-error', min_bw, bw_uns)
  bw <- max(bw_stim, min_bw, bw_uns)
  dens_stim <- density(cut_stim$expr, bw = bw)
  dens_uns <- density(
    cut_uns$expr, from = min(dens_stim$x),
    to = max(dens_stim$x), bw = bw
  )

  # put raw densities into table
  dens_tbl_raw <- tibble::tibble(
    x_stim = dens_stim$x, y_stim = dens_stim$y
  ) |>
    dplyr::mutate(
      y_uns = purrr::map_dbl(
        x_stim,
        function(marker) {
          .interp(
            val = marker,
            x = dens_uns$x,
            y = dens_uns$y
          )
        }
      )
     ) |>
    #dplyr::mutate(y_uns = .env$y_uns) |>
    tidyr::pivot_longer(
      y_stim:y_uns,
      names_to = "stim",
      values_to = "dens"
      ) |>
    dplyr::mutate(
      stim = ifelse(stim == "y_stim", "yes", "no")
    )

  .debug(debug, "Normalising probabilities")

  # calculate raw and
  # normed probability based on densities for density measurements
  prob_tbl <- dens_tbl_raw |>
    tidyr::pivot_wider(
      id_cols = x_stim,
      names_from = stim,
      values_from = dens
      ) |>
    dplyr::mutate(
      prob_stim = 1 - no/yes,
      prob_stim = ifelse(yes == 0 & no == 0, 0, prob_stim),
      prob_stim_norm = pmin(1, prob_stim),
      prob_stim_norm = pmax(0, prob_stim_norm)
    ) |>
    dplyr::filter(x_stim > cp_min)

  .debug(debug, "Filtering before smoothing")

  density_exc_min <- density(cut_stim_orig$expr)
  dens_tbl <- tibble::tibble(
    x = density_exc_min$x,
    y = density_exc_min$y
  )
  peak <- dens_tbl |>
    dplyr::filter(y == max(y)) |>
    dplyr::pull("x")

  prob_tbl <- prob_tbl |>
    dplyr::filter(
      x_stim > peak + 0.02 * diff(range(x_stim))
    )

  # get range of values for which we'd want to calculate
  # probability:
  # - those cells for which the probability
  # of responding from their position onwards
  # is 0.025 or more
  prob_tbl_pos <- prob_tbl |>
    dplyr::mutate(
      ge10 = prob_stim_norm >= 0.025,
      ge10 = cumsum(ge10) > 0
      ) |>
    dplyr::filter(ge10)

  # we then look at the minimum protein expression to achieve this
  min_prob_x <- min(prob_tbl_pos$x_stim)

  # we set the threshold larger than the
  # maximum observed value if the following conditions are met:
  # - there are no cells that mark the start of responding in
  # the way defined above
  # - the maximum value of the stimulated sample is less
  # than the minimum probability
  #   - not clear why we do this?...
  #     - Well, it just means there's no response, #
  #       as teh predicted probability is too high
  # why did we choose the imputed threshold we chose, then?
  # what is the threshold?
  # - well, it's a bit to the right of the mximum
  # - and we set it quite high - half the difference in the
  # original range!
  #   - I guess sometimes this could help, and we don't lose out
  # too much if we gate more aggressively later on
  #   - I suppose we set it high to avoid the gates being too low
  if (nrow(prob_tbl_pos) == 0 || max(cut_stim_orig$expr) < min_prob_x) {
    .debug(debug, "No responding cells")
    return(list(
      cp = max(
        cp_min,
        max(cut_stim_orig$expr) + 0.5 * diff(range(cut_stim_orig$expr))
        ),
      p_list = p_list_empty
    ))
  }

  # okay, so here get the probability
  # that was observed at the point
  # where the protein is lowest and included
  min_prob <- min(
    prob_tbl_pos$prob_stim[prob_tbl_pos$x_stim == min_prob_x]
  )
  # we then do something funky
  # - take a 40th of the 6 - what?
  # okay, so that's just the typical range of expression
  # for our case, it was 6.
  # we could probably do differenty.
  # we then take only values in the probability
  # table for which the protein measuremnts are
  # a bit greater than the minimum
  # probability
  # we then filter out the
  # model in

  range_len <- 6; margin <- 0.025 * range_len
  prob_tbl_pos_margin <- prob_tbl |>
    dplyr::filter(x_stim > min_prob_x - margin)
  data_mod <- cut_stim_orig |>
    dplyr::filter(expr > min(min_prob_x) - margin)

  if (nrow(prob_tbl_pos) == 0) {
    .debug(debug, "No responding cells")
    return(
      list(
        cp = max(cut_stim_orig$expr) + 0.5 * diff(range(cut_stim_orig$expr)),
        p_list = p_list_empty
      )
    )
  }

  min_prob_x <- min(data_mod$expr)
  max_prob_x <- max(data_mod$expr)

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
  data_mod <- data_mod |>
    dplyr:::mutate(prob_smooth = expr)
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
        }))
  }
  

  if (nrow(data_mod) >= 10){
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

    if (!inherits(fit, "try-error")){
      pred_vec <- predict(
        fit, type = 'response'
      )
      mean_abs_error <- mean(
        abs(pred_vec - data_mod$prob_smooth)
      )
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
            devtol.fit = 0.01#,
            #steptol.fit = 1e-3,
            #bfgs = list(steptol.bfgs = 1e-4),
            #maxit = 5
            )),
            silent = TRUE
        )

      if (!inherits(fit, "try-error")) {
        pred_vec <- predict(fit, type = 'response')
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
        pred_vec <- predict(fit, type = 'response')
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
      dplyr::mutate(pred = prob_smooth - 0.0001)
  } 

  if (FALSE) {
    data_plot <- data_mod |>
      dplyr::select(expr, prob_smooth, pred) |>
      dplyr::rename(x_stim = expr) |>
      tidyr::pivot_longer(names_to = "type",
                   values_to = "prob",
                   -x_stim) |>
      dplyr::bind_rows(prob_tbl_pos |>
                  dplyr::select(x_stim, prob_stim_norm) |>
                  dplyr::mutate(type = 'prob_stim_norm') |>
                  dplyr::rename(prob = prob_stim_norm))

    ggplot(data_plot, aes(x = x_stim, col = type)) +
      geom_line(aes(y = prob)) +
      geom_rug(data = data_plot |>
                 dplyr::filter(type == 'pred')) +
      scale_colour_brewer(palette = "Accent")
  }

  data_count <- data_mod |>
    dplyr::filter(expr >= min_prob_x) |>
    dplyr::arrange(expr) |>
    dplyr::mutate(n_row = seq_len(dplyr::n())) |>
    dplyr::filter(cumsum(pred > prob_smooth) != n_row)

  count <- sum(data_count$pred)

  prop_bs_est <- count/nrow(cut_stim_orig)

  data_threshold <- data_count |>
    dplyr::arrange(desc(expr)) |>
    dplyr::mutate(count_stim = seq_len(dplyr::n())) |>
    dplyr::mutate(
      prop_stim = count_stim/nrow(cut_stim_orig),
      prop_uns = purrr::map_dbl(expr, function(x){
        sum((cut_uns$expr - bias) >= x) / nrow(cut_uns_orig)
      }),
      prop_bs = prop_stim - prop_uns,
      prop_bs_diff = prop_bs - prop_bs_est,
      prop_stim_pos = pmax(prop_stim, 0.5 / nrow(cut_stim_orig)),
      prop_uns_pos = pmax(prop_uns, 0.5 / nrow(cut_uns_orig)),
      prop_stim_sd =
        sqrt(prop_stim_pos * (1 - prop_stim_pos)) /
          nrow(cut_stim_orig),
      prop_uns_sd =
        sqrt((prop_uns_pos * (1 - prop_uns_pos)) / 
          nrow(cut_uns_orig)),
      prop_bs_sd =
        sqrt(prop_stim_sd^2 + prop_uns_sd^2)
    )

  if (FALSE) {
    ggplot(data_threshold |>
              tidyr::pivot_longer(names_to = "type",
                          values_to = "prop",
                          c(prop_stim, prop_uns, prop_bs)),
            aes(x = expr, y = prop * 1e2, col = type)) +
      geom_line()
  }

  min_diff_tbl <- data_threshold |>
    dplyr::filter(
      abs(prop_bs_diff) ==  min(abs(prop_bs_diff))
    ) |>
    dplyr::slice(1)

  min_diff <- min_diff_tbl$prop_bs_diff
  min_diff_bs_sd <- min_diff_tbl$prop_bs_sd
  min_diff_sd <- min_diff/min_diff_bs_sd

  cp <- min_diff_tbl$expr
  # min_diff_tbl
  #min_diff_tbl$prop_bs * 1e2

  # check that marginal gain in signal-to-noise ratio is positive
  # if not, move back up until a maximum of margin units
  #data_threshold |>
  #  dplyr::mutate(prop_bs_sd_std = sqrt(prop_bs_sd * pmin(prop_bs_sd, min_diff_bs_sd))) |>
  #  dplyr::mutate(prop_bs_diff_sd_std = abs(prop_bs_diff)/prop_bs_sd_std) |>
  #  dplyr::select(expr, prop_stim, prop_uns, prop_bs, prop_bs_diff, prop_bs_diff_sd_std) |>
  #  dplyr::filter(prop_bs_diff_sd_std <= min_diff_sd + 0.75 * min_diff_bs_sd) |>
  #  dplyr::arrange(expr) |>
  #  dplyr::slice(1)

  # explore s2n
  #data_threshold |>
  #  dplyr::mutate(s2n = prop_bs/prop_bs_sd) |>
  #  dplyr::mutate(s2n_marginal = diff(c(0, s2n))) |>
  #  dplyr::select(expr, prop_bs, prop_bs_sd, s2n, s2n_marginal)



  if (plot && FALSE) {
    p_list <- .plot_cp_loc_fdr(
      dens_tbl = dens_tbl_raw,
      params = params,
      prob_tbl = prob_tbl,
      pred_tbl = all_cell_pred_tbl,
      cut_stim = cut_stim_orig$expr,
      sample = cut_stim_orig$sample[1],
      min_x_pos_prob = min(all_cell_pred_tbl$cut_stim),
      path_project = path_project
    )
  } else p_list <- p_list_empty

  #print('done getting loc fdr gate for this single sample')

  .debug(debug, "Completed loc gate for single sample")

  list(cp = cp, p_list = p_list_empty)

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

  dir_base <- stim_gate_dir_base_create(
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
    dplyr::mutate(prob_type = 'pred') |>
    dplyr::select(x_stim, prob_type, prob )

  plot_tbl <- plot_tbl_prob |>
    dplyr::bind_rows(plot_tbl_pred)
  p_loc_prob <- ggplot(plot_tbl,
                       aes(x = x_stim, y = prob, linetype = prob_type)) +
    cowplot::theme_cowplot(font_size = 20) +
    cowplot::background_grid(major = 'y', minor = 'y') +
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
  # new_pred_tbl <- tibble::tibble(x_stim = cut_stim[cut_stim >= min_x_pos_prob])
  # new_pred_vec <- predict(prob_mod, newdata = new_pred_tbl, type = 'response')
  # new_pred_tbl <- new_pred_tbl |> dplyr::mutate(pred = new_pred_vec)

  #  n_bin <- hist(new_pred_tbl$x_stim, plot = FALSE)$mids |> length
  n_bin <- hist(pred_tbl$cut_stim, plot = FALSE)$mids |> length()
  if(n_bin == 1){
    bin_tbl <- pred_tbl |>
      dplyr::mutate(bin = paste0("(", cut_stim, ",", cut_stim, "]"),
             ctb = pred)
  } else{
    bin_tbl <- pred_tbl |>
      dplyr::mutate(bin = cut(.data$cut_stim, breaks = n_bin)) |>
      dplyr::group_by(bin) |>
      dplyr::summarise(ctb = sum(pred))
  }

  #' @title Get the midpoint of a bin returned by base::cut
  .get_bin_mid <- function(bin){
    purrr::map_dbl(bin, function(bin_curr){
      comma_loc <- stringr::str_locate(bin_curr, ",")[1,"start"][[1]]
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
    geom_bar(stat = 'identity',
             col = 'gray25',
             fill = 'gray85')

  dir_save <- file.path(dir_base, "p_loc_ctb")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  cowplot::ggsave2(filename = file.path(file.path(dir_save, fn)),
                   height = 10, width = 10)

  #list(p_loc_dens = p_loc_dens, p_loc_prob = p_loc_prob,
  #     p_loc_ctb = p_loc_ctb)
  list()

}

.get_prob_mod <- function(prob_tbl, nrow_level){

  prob_mod <- try(scam::scam(prob_stim_norm ~ s(x_stim, bs = "mpi"),
                         family = "binomial",
                         data = prob_tbl |>
                           dplyr::filter(prob_stim_norm > 0.25)))

  if(class(prob_mod) == 'try-error'){
    prob_mod <- try(scam::scam(prob_stim_norm ~ s(x_stim, bs = "micx"),
                               family = "binomial",
                               data = prob_tbl))
  } else return(prob_mod)

  if(class(prob_mod) == 'try-error'){
    prob_mod <- try(mgcv::gam(prob_stim_norm ~ s(x_stim),
                               family = "binomial",
                               data = prob_tbl))
  }

  prob_mod
}

.get_cell_specific_response_probs <- function(prob_tbl, nrow_level, params,
                                              prob_min, cut_stim){

  # =====================
  # Fit model
  # =====================

  # fit a model to get smoothed probabilities
  prob_mod <- .get_prob_mod(prob_tbl, nrow_level = nrow_level)

  # if model fails, then return it
  if('ultimate error' %in% class(prob_mod)){

    out <- list(pred_tbl = NA,
                prob_mod = NA)

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

  out <- list(pred_tbl = pred_tbl,
              prob_mod = prob_mod)

  out
}


