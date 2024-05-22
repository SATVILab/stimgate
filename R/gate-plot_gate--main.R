#' @title Plot gating-related plots for a list of cp objects
#'
#' @description
#' Creates the plots needed to inspect the cutpoints selected for a
#' given cp object. If the cp_obj was gated as a batch, then
#' the plots are created for each sample that went into the batch as well.
#'
#' @param cp_obj object of class \code{cp}. Output from get_cp.
#' @param save logical. If \code{TRUE}, then the plots are saved. Default is \code{TRUE}.
#' @param exc_min logical. If \code{TRUE}, then all cells with values for at least one
#' of the plotted channels equal to the minimum across any cell in those channels are excluded
#' from the plots.
#' @param node character vector. If provided, then plots are made using the cutpoints provided by \code{cp_obj_list} for each of the populations specified in node, in addition to the population initially gated on. If \code{NULL}, then  only the population initially gated on has plots made for it.
#' @param cp named numeric vector. The names are the names of the cutpoint methods and the values are the actual cutpoints. Seems incorrect - below seems right.
#' @param cp character vector. Names of cutpoints to plot.
#' @param ind_uns_in_batch integer. Integer specifying which index of a given batch is unstim. For example, if ind_uns is 3, then the third element in cp$params$ind is taken to be the unstim distribution. If cp$params$ind were 26:30, then 28 would be the index of the unstim distribution for that batch in cp$data.
#' @param auto_y_lim logical. If \code{TRUE}, the limits of the y-axis for the count and frequency plots are set equal to double the count and frequency, respectively, of the lowest specified cutpoint.
#' @param save logical. If \code{TRUE}, then the plots are saved as PNGs. Default is \code{TRUE}.
#' @param save_r_obj logical. If \code{TRUE}, then the plots are saved as an R-readable object. Default is \code{FALSE}.
#' @param return logical. If \code{TRUE}, then the list of plots generated are returned. If \code{FALSE}, then the function returns \code{TRUE} invisibly. Default is \code{FALSE}.

#' @return A list, where each element is of class \code{cp_plot}.
#'
#' @param cp_obj_list list. List where each element is an object of class \code{cp}.
plot_gate <- function(gate_tbl, params, pop_sub,
                      gate_lab, col) {

  # get populations to get results for
  pop_res_vec_full <- c(
    params$pop_gate, paste0(params$pop_gate, pop_sub)
  ) |> unique()
  pop_res_name_vec <- c("gate", pop_sub)

  # update_and_save_locb_plots_with_lines(params = params, gate_tbl = gate_tbl,
  #                                        col = col, gate_lab = gate_lab)

  # get plots for each batch
  plot_list <- purrr::walk(seq_along(pop_res_vec_full), function(i) {
    pop_res_full_curr <- pop_res_vec_full[i]
    print("saving plots")
    purrr::walk(seq_along(params$ind_batch_list), function(j) {
      batch_curr <- params$ind_batch_list[[j]]

      plot_list_plots <- .plot_gate_batch(ind_batch = batch_curr,
                                          gate_tbl = gate_tbl |>
                                            dplyr::filter(gate_type != 'ctrl'),
                                          params = params,
                                          gate_lab = gate_lab, col = col,
                                          pop_res = pop_res_full_curr)
      plot_list_plots_batch <- stats::setNames(
        list(plot_list_plots), paste0(batch_curr, collapse = "_")
      )
      plot_list_plots_batch_pop <- stats::setNames(
        list(plot_list_plots_batch), pop_res_name_vec
      )
      #debugonce(.save_plot_gate)
      .save_plot_gate(plot_list = plot_list_plots_batch_pop,
                      params = params)
    }) #|>
      #stats::setNames(purrr::map_chr(params$ind_batch_list, function(x) paste0(x, collapse = "_"))
  }) #|>
  #stats::setNames(pop_res_name_vec)

  # save plots

  print("done saving gate plots")


  invisible(TRUE)
}




#' @title Plot gating-related plots for a single cp object
#'
#' @description
#'
#'
#' @inheritParams plot_cp
#' @param cp_obj object of class cp_obj.
#'
#' @return An object of class \code{cp_plot}. This object is a list,
#' which consists of a set of plots for each sample that went into creating
# a set of cutpoint(s) and the set of plots for the entire batch as well.
#' Within each such element are the following plots: plot_hist (histogram of channel),
#' plot_ll_scp (log-likelihood for first changepoint),
#' plot_prop (binned probability of positivity across expression range),
#' plot_ll_dcp (log-likelihood for second changepoint),
#' plot_2d (list of hex plots of cut channel against high channels),
#' plot_dcp_hist (histogram of marker expression above first changepoint) and
#' plot_count_stats_by_cp (list of plots of names 'p_prop', 'p_count', 'p_freq'
#' and 'p_std_prop' that show various stats for a range of cutpoints).
.plot_gate_batch <- function(gate_tbl, ind_batch, params,
                       gate_lab, col, pop_res) {

  # ==============================
  # Preparation
  # ==============================

  # spread nearby obs
  # gate_tbl <- .spread_nearby_obs_tbl(gate_tbl = gate_tbl,
  #                                   width = max(gate_tbl$gate)/400)

  # get plot indices
  if(params$data_name == "gs_proto" & identical(ind_batch, 26:30)) {
    ind_in_batch_lab_vec_curr <- params$ind_in_batch_lab_vec[c(3, 2, 1, 4, 5)]
    ind_plot_list <- stats::setNames(
      ind_batch,
      ind_in_batch_lab_vec_curr)[
        union(params$ind_in_batch_gate, params$ind_in_batch_uns)
        ] |>
      as.list()
  } else{
    ind_plot_list <- stats::setNames(
      ind_batch,
      params$ind_in_batch_lab_vec)[
        union(params$ind_in_batch_gate, params$ind_in_batch_uns)
      ] |>
      as.list()
  }

  ind_plot_list <- ind_plot_list |>
    append(list(grp = ind_plot_list |> stats::setNames(NULL) |> unlist()))

  # get response tibble as lists
  ex_list <- .get_ex_list(data = params$data, ind_batch = ind_batch,
                          ind_in_batch_gate = params$ind_in_batch_gate,
                          ind_in_batch_uns = params$ind_in_batch_uns,
                          ind_in_batch_lab_vec = params$ind_in_batch_lab_vec,
                          pop = pop_res, cut = params$cut, high = params$high,
                          data_name = params$data_name)

  # create combined ex
  ex <- dplyr::bind_rows(ex_list) |> tibble::as_tibble()

  # ==============================
  # Plot
  # ==============================

  plot_gate_batch_list <- purrr::map(ind_plot_list, function(ind_plot) {
    #print(ind_plot)
    #if(ind_plot == '19') debugonce(.plot_gate_ind)
    .plot_gate_ind(ex = ex |> dplyr::filter(is_uns | ind %in% ind_plot),
                  params = params,
                  gate_tbl = gate_tbl |> dplyr::filter(ind %in% ind_plot),
                  gate_lab = gate_lab,
                  col = col,
                  ind_plot = ind_plot)
  }) |>
    stats::setNames(paste0(ex[['batch_sh']][1], "_", names(ind_plot_list)))

  plot_gate_batch_list
}

#' @inheritParams plot_cp
#' @inheritParams .plot_cp
#' @param node character. Node for which to make the plots. If provided, then
#' the plots are made for this cell population rather than the population the cell
#' was gated for. Also if provided, then the log-likelihood plots are not longer applicable,
#' and so those plots are then not made.
#' @param cp named numeric vector. Names are names of cutpoints, and values are associated cutpoints.
#' @return A list, which consists of a set of plots for a given sample associated with a
#' given set of cutpoints.
#' Within this list are the following plots: plot_hist (histogram of channel),
#' plot_ll_scp (log-likelihood for first changepoint),
#' plot_prop (binned probability of positivity across expression range),
#' plot_ll_dcp (log-likelihood for second changepoint),
#' plot_2d (list of hex plots of cut channel against high channels),
#' plot_dcp_hist (histogram of marker expression above first changepoint) and
#' plot_count_stats_by_cp (list of plots of names 'p_prop', 'p_count', 'p_freq'
#' and 'p_std_prop' that show various stats for a range of cutpoints).
.plot_gate_ind <- function(ex, params, gate_tbl, gate_lab, col, ind_plot) {

  # ----------------------
  # Plots - base
  # ----------------------

  # get whether to exclude minimum values or note, if data cytof or not
  exc_min <- ifelse(str_detect_any(
    params$data_name, c(
      "gs_cytof", "gs_proto", "gs_cd8_base",
      "gs_cytof_acs"
    )
  ), TRUE, FALSE)


  # get the minimum value to plot for plot_hist_tail and plot_stats_by_cp
  lb <- quantile(
    ex[["cut"]],
    max(c(
      0, 1 - 2 * (1 - ecdf(ex[["cut"]])(min(gate_tbl$gate))),
      1 - (1 - ecdf(ex[["cut"]])(min(gate_tbl$gate))) - 0.05
    ))
  )[1]
  # ----------------------
  # Plots - base
  # ----------------------

  # plot histogram of cut channel
  #plot_hist <- .plot_hist(ex = ex |> dplyr::filter(ind %in% ind_plot),
  #                        plot_var = params$cut,
  #                        gate_tbl = gate_tbl,
  #                        exc_min = exc_min,
  #                        col = col, gate_lab = gate_lab,
  #                        axis_lab = params$chnl_lab)

  # plot first cutpoint log-likelihood
  #plot_ll_scp <- .plot_ll(ll = cp_obj$ll_scp,
  #                        cps = params$cps_scp,
  #                        cp = cp,
  #                        cut = params$cut,
  #                        axis_lab = params$chnl_lab,
  #                        col = col,
  #                        cp_lab = cp_lab)

  # plot proportion positive
  # plot_prop <- .plot_prop(high_ind_tbl,
  #                         cp = cp,
  #                        cut = params$cut,
  #                        axis_lab = params$chnl_lab,
  #                        col = col,
  #                        cp_lab = cp_lab,
  #                        rug = rug)

  # plot second cutpoint log-likelihood
  # plot_ll_dcp <- .plot_ll(ll = cp_obj$ll_dcp,
  #                        cps = as.numeric(names(cp_obj$ll_dcp)),
  #                        cp = cp, cut = params$cut,
  #                        axis_lab = params$chnl_lab, col = col,
  #                        cp_lab = cp_lab)

  # get 2d plots of cut channel against each of the high channels
  plot_2d <- .plot_2d(ex = ex |> dplyr::filter(ind %in% ind_plot),
                      plot_var_x = params$cut,
                      high = params$high,
                      gate_tbl = gate_tbl |> dplyr::filter(chnl == params$cut),
                      axis_lab = params$chnl_lab,
                      exc_min = exc_min,
                      col = col,
                      gate_lab = gate_lab)

  # plot histogram of cut channel filtered on cp_pw cut with cp_dcp cut overlaid
  # plot_hist_tail <- .plot_hist(ex = ex |> dplyr::filter(ind %in% ind_plot),
  #                             plot_var = params$cut,
  #                             gate_tbl = gate_tbl,
  #                             exc_min = exc_min,
  #                             col = col, gate_lab = gate_lab,
  #                             axis_lab = params$chnl_lab,
  #                             lb = lb)

  # plot
  # plot_logit <- .plot_logit(high_ind_tbl = high_ind_tbl, cp = cp,
  #                           plot_prop = plot_prop)

  # plot count stats by cutpoint
  # plot_count_stats_by_cp <- .plot_count_stats_by_cp(ex = ex,
  #                                                  plot_var = params$cut,
  #                                                  gate_tbl = gate_tbl,
  #                                                  col = col,
  #                                                  gate_lab = gate_lab,
  #                                                  axis_lab = params$chnl_lab,
  #                                                  lb = lb)

  # that function Tom wants of the stim vs unstim as bars for each person

  # ----------------------
  # Output
  # ----------------------

  out_list <- list(plot_hist = list(), #plot_hist,
                   plot_hist_tail = list()) |>#plot_hist_tail) |>
    append(plot_2d) #|>
    #append(plot_count_stats_by_cp)

  class(out_list) <- 'cp_plot'

  out_list
}

#' @title Save gating-related plots
#'
#' @description
#'
#' @param params named list. List containing elements with the
#' following names and values:
#' \begin{describe}{
#'   \item{cut}{character. Name of channel to gate on.}
#'   \item{chnl_lab}{named character vector. Names are channel names
#'   and values are channel descriptions (typically markers), e.g.
#'   \code{c("Nd146Di" = "TNFa", "Lu175Di" = "Perforin")}.
#'   \item{pop_gate}{character. Full name of population in \code{GatingSet} on
#'   which to set the gate.}
#'   \item{data_name}{character. Name of \code{GatingSet}.}
#'   \item{ind_in_batch_gate}{numeric vector. For a given batch, specifies the indices
#'   within the batch for which gates are required.}
#'   \item{ind_in_batch_lab_vec}{character vector. Specifies the names for the samples within
#'   a batch. Useful if they are, for instance, different stim conditions.}
#' }
#' @param plot_list named list.
#'
#' @return Returns \code{TRUE} invisibly.
#'
#' @examples
#' params <- list(cut = "Ho165Di", chnl_lab = c("Ho165Di" = "IFNg", "Nd146Di" = "TNFa"),
#'                pop_gate = "/CD33- CD3+/CD20-/CD14-/TCRgd-/NKT-/CD4 T Cells",
#'                data_name = "gs_cytof",
#'                ind_in_batch_gate = 1:3,
#'                ind_in_batch_lab_vec  = c("mtbaux", "ebv", "uns"))
#'  plot_list <- list("CD4 T cells" = list("p_hist" = list("1" = p_hist_1,
#'                                                         "2" = p_hist_2),
#'                                                         "3" = p_hist_3),
#'                                          "p_2d_tnfa" = list("))
#' .save_plot_gate(params = params, plot_list = )
.save_plot_gate <- function(params, plot_list, dir_base_init = NULL) {

  # --------------------------
  # Preparation
  # --------------------------

  # create base directory if need be
  dir_base <- stimgate_dir_base_create(params = params, dir_base_init = dir_base_init)
  dir_save <- file.path(dir_base, "gating_plots")
  #if(dir.exists(dir_save)) unlink(dir_save, recursive = TRUE, force = TRUE)
  if(!dir.exists(dir_save)) dir.create(dir_save)

  # --------------------------
  # Save plots as PNGs
  # --------------------------

  plot_type_vec <- names(plot_list[[1]][[1]][[1]])
  batch_stim_vec <- names(plot_list[[1]][[1]])
  ind_batch_vec <- names(plot_list[[1]])
  pop_res_vec <- names(plot_list)
  pop_res_vec <- stringr::str_remove_all(pop_res_vec,"/") |>
    stringr::str_remove_all(" ")

  # loop over the populations to plot for
  plot_list <- purrr::map(pop_res_vec, function(pop_res) {
    # loop over the populations to get types for
    purrr::map(plot_type_vec, function(plot_type) {
      # get names of batches
      batch_name_vec <- purrr::map_chr(
        names(plot_list[[pop_res]]),
        function(ind_batch) {
          purrr::map_chr(
            names(plot_list[[pop_res]][[ind_batch]]),
            function(x) stringr::str_split(x, "_")[[1]][1]
          )[1]
        }
      )

      ind_batch_name_vec <- purrr::map_chr(
        seq_along(batch_name_vec),
        function(i) {
        paste0(names(plot_list[[pop_res]])[i], "-", batch_name_vec[i])
      })

      # loop over batches
      purrr::map(ind_batch_vec, function(ind_batch) {
        # get names of individual samples
        plot_list_pop_ind <- plot_list[[pop_res]][[ind_batch]]
        batch_stim_vec <- names(plot_list_pop_ind)
        sample_name_vec <- purrr::map_chr(
          batch_stim_vec, function(x) stringr::str_split(x, "_")[[1]][2]
        )
        # get plots from given pop of given type
        # from given batch for all samples in batch
        purrr::map(batch_stim_vec, function(batch_stim) {
          plot_list_pop_ind[[batch_stim]][[plot_type]]
         }) |>
            stats::setNames(sample_name_vec)
      }) |>
        stats::setNames(ind_batch_name_vec)
    }) |>
      stats::setNames(plot_type_vec)
  }) |>
    stats::setNames(pop_res_vec)

  # save plots in groups above
  purrr::walk(pop_res_vec, function(pop_res) {
    print(pop_res)
    purrr::walk(plot_type_vec, function(plot_type) {
      print(plot_type)
      # get directory to save to
      # node_save <- ifelse(node_curr == nodes[1], "base", stringr::str_remove(node_curr, nodes[1]))
      # node_save <- stringr::str_replace_all(node_save, "[[//]]", "_")
      dir_save_curr <- file.path(dir_save, pop_res, plot_type)

      # create directory to save to if it doesn't yet exist
      #if(dir.exists(dir_save_curr)) unlink(dir_save_curr, recursive = TRUE)
      if(!dir.exists(dir_save_curr)) dir.create(dir_save_curr, recursive = TRUE)

      # create plots and save them

      purrr::walk(seq_along(plot_list[[pop_res]][[plot_type]]), function(j) {
        plot_list_grp <- plot_list[[pop_res]][[plot_type]][[j]]

        plot_list_grp <- suppressWarnings(
          .equalise_plot_lims(
            plot_list = plot_list_grp,
            plot_type = plot_type
          )
        )
        plot_grp <- cowplot::plot_grid(
          plotlist = plot_list_grp,
          labels = names(plot_list_grp),
          ncol = 3
        )

        plot_name <- stringr::str_split(
          names(plot_list[[pop_res]][[plot_type]])[j], "-"
        )[[1]][2]

        cowplot::ggsave2(
          filename = file.path(dir_save_curr, paste0(plot_name, ".png")),
          plot = plot_grp, height = 30, width = 55, units = "cm"
        )
      })
    })
  })

  invisible(TRUE)
}


#' @title Get axis limits for
#' x- and possibly y-axes for a range of plot
#' types for groups of plots
.get_plot_lims <- function(plot_list, plot_type) {

  # histogram plots
  if(plot_type %in% c('plot_hist', 'plot_hist_tail')) {
    lims_df <- purrr::map_df(plot_list, function(p) {

      # do nothing if plot blank
      if(length(p$layers) == 0) return(data.frame(min = NA, max = NA))

      # histogram limts
      data_x <- p$data$x
      lims_vec <-  range(data_x)

      # cutpoint limits
      cp_x <- p$layers[[2]]$data$cp
      lims_vec <- c(min(lims_vec, cp_x),
                    max(lims_vec, cp_x))

      # limits limits
      lims_lims_vec <- p$coordinates$limits$x
      if (!is.null(lims_lims_vec)) {
        lims_vec <- c(
          max(c(lims_vec[1], lims_lims_vec)),
          min(c(lims_vec[2], lims_lims_vec))
        )
      }

      data.frame(min = min(lims_vec), max = max(lims_vec))
    }) |>
      dplyr::filter(!is.na(min))

    if(!nrow(lims_df)) return(list(x = rep(NA, 2)))

    return(list(x = c(min(lims_df$min), max(lims_df$max))))
  }

  # prop plot
  if(plot_type == 'plot_prop') {

    lims_df <- purrr::map_df(plot_list, function(p) {

      # do nothing if plot blank
      if(length(p$layers) == 0) return(data.frame(min = NA, max = NA))

      # prop limits
      data_prop_x <- p$layers[[1]]$data$mid
      lims_vec <- range(data_prop_x)

      # rug limits
      rug_x <- p$layers[[3]]$data$cut
      lims_vec <- c(min(lims_vec, rug_x),
                    max(lims_vec, rug_x))

      # cutpoint limits
      cp_x <- p$layers[[4]]$data$cp
      lims_vec <- c(min(lims_vec, cp_x),
                    max(lims_vec, cp_x))

      # limits limits
      lims_x <- p$coordinates$limits$x
      if(!is.null(lims_x)) lims_vec <- c(max(c(lims_vec[1], lims_x)),
                                         min(c(lims_vec[2], lims_x)))

      # return
      data.frame(min = min(lims_vec), max = max(lims_vec))
    }) |>
      dplyr::filter(!is.na(min))

    if(!nrow(lims_df)) return(list(x = rep(NA, 2)))

    return(list(x = c(min(lims_df$min), max(lims_df$max))))
  }

  # 2d plots
  if(stringr::str_detect(plot_type, "plot2d")) {
    lims_df_x <- purrr::map_df(plot_list, function(p) {

      # do nothing if plot blank
      if(length(p$layers) == 0) return(data.frame(min = NA, max = NA))

      # hex limits
      hex_x <- p$layers[[1]]$data$cut
      lims_vec <- range(hex_x)

      # cutpoint limits
      cp_x <- p$layers[[2]]$data$cp
      lims_vec <- c(min(lims_vec, cp_x),
                    max(lims_vec, cp_x))

      # limits limits
      lims_lims_vec <- p$coordinates$limits$x
      if (!is.null(lims_lims_vec)) {
        lims_vec <- c(
          max(c(lims_vec[1], lims_lims_vec)),
          min(c(lims_vec[2], lims_lims_vec))
        )
      }

      data.frame(min = min(lims_vec), max = max(lims_vec))
    }) |>
      dplyr::filter(!is.na(min))

    lims_df_y <- purrr::map_df(plot_list, function(p) {

      # do nothing if plot blank
      if(length(p$layers) == 0) return(data.frame(min = NA, max = NA))

      # hex limits
      hex_y <- p$layers[[1]]$data$high
      lims_vec <- range(hex_y)

      # cutpoint limits
      lims_vec <- c(lims_vec[1], max(lims_vec, p$layers[[2]]$data$min))

      # limits limits
      lims_lims_vec <- p$coordinates$limits$y
      if (!is.null(lims_lims_vec)) {
        lims_vec <- c(
          max(c(lims_vec[1], lims_lims_vec)),
          min(c(lims_vec[2], lims_lims_vec))
        )
      }

      data.frame(min = min(lims_vec), max = max(lims_vec))
    }) |>
      dplyr::filter(!is.na(min))

    if(!nrow(lims_df_x)) return(list(x = rep(NA, 2),
                                     y = rep(NA, 2)))

    return(list(x = c(min(lims_df_x$min), max(lims_df_x$max)),
                y = c(min(lims_df_y$min), max(lims_df_y$max))))
  }

  # stat plots
  if(stringr::str_detect(plot_type, "plot_cp_stat")) {
    lims_df <- purrr::map_df(plot_list, function(p) {

      # do nothing if plot blank
      if(length(p$layers) == 0) return(data.frame(min = NA, max = NA))

      # histogram limts
      data_x <- p$data$cp
      lims_vec <-  range(data_x)

      # cutpoint limits
      cp_x <- p$layers[[3]]$data$cp
      lims_vec <- c(
        min(lims_vec, cp_x),
        max(lims_vec, cp_x)
      )

      # limits limits
      lims_lims_vec <- p$coordinates$limits$x
      if (!is.null(lims_lims_vec)) {
        lims_vec <- c(
          max(c(lims_vec[1], lims_lims_vec)),
          min(c(lims_vec[2], lims_lims_vec))
        )
      }

      data.frame(min = min(lims_vec), max = max(lims_vec))
    }) |>
      dplyr::filter(!is.na(min))

    if(!nrow(lims_df)) return(list(x = rep(NA, 2),
                                   y = rep(NA, 2)))

    return(list(x = c(min(lims_df$min), max(lims_df$max))))
  }

  # histogram plots
  if(plot_type == 'plot_logit') {
    lims_df <- purrr::map_df(plot_list, function(p) {

      # do nothing if plot blank
      if(length(p$layers) == 0) return(data.frame(min = NA, max = NA))

      # histogram limts
      data_x <- p$data$x
      lims_vec <-  range(data_x)

      # cutpoint limits
      cp_x <- p$layers[[2]]$data$cp
      lims_vec <- c(min(lims_vec, cp_x),
                    max(lims_vec, cp_x))

      # limits limits
      lims_lims_vec <- p$coordinates$limits$x
      if (!is.null(lims_lims_vec)) {
        lims_vec <- c(
          max(c(lims_vec[1], lims_lims_vec)),
          min(c(lims_vec[2], lims_lims_vec))
        )
      }

      data.frame(min = min(lims_vec), max = max(lims_vec))
    }) |>
      dplyr::filter(!is.na(min))

    if(!nrow(lims_df)) return(list(x = rep(NA, 2)))

    return(list(x = c(min(lims_df$min), max(lims_df$max))))
  }
}

#' @title Ensure that a group of plots have at least the same x-axis limits
.equalise_plot_lims <- function(plot_list, plot_type) {

  # do nothing if plot type is a log-likelihood plot (all limits same anyway)
  if(plot_type %in% c("plot_ll_scp", "plot_ll_dcp")) return(plot_list)

  # get limits for x and possibly y, depending on plot type
  plot_lims_list <- .get_plot_lims(plot_list = plot_list, plot_type = plot_type)

  # return plot list as is if the limits are all NA or list is NULL
  if(is.null(plot_lims_list)) return(plot_list)
  if(is.na(plot_lims_list$x[1])) return(plot_list)

  # readjust the limits if both x and y are specified
  if(length(setdiff(c("x", "y"), names(plot_lims_list))) == 0) {
    plot_list <- purrr::map(plot_list, function(p) p +
                              coord_cartesian(xlim = plot_lims_list$x,
                                              ylim = plot_lims_list$y))
    # readjust the limits if only x is specified
  } else{
    plot_list <- purrr::map(plot_list, function(p) p +
                              coord_cartesian(xlim = plot_lims_list$x))
  }
  plot_list
}








