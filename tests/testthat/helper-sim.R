# fs
get_fs <- function() {
  tryCatch(
    # don't want a package warning if file does not exist,
    # clearly handling such an error
    suppressWarnings(readRDS(
      testthat::test_path(
        "../../inst/extdata/bodenmiller_bcr_xl_fs.rds"
        ))),
    error = function(e) {
      path_hub <- file.path(Sys.getenv("HOME"), ".cache/R/ExperimentHub")
      if (!dir.exists(path_hub)) {
        dir.create(path_hub, recursive = TRUE)
      }
      # get an odd flowSet slot dropping warning,
      # which is unrelated to our functionality
      suppressWarnings(HDCytoData::Bodenmiller_BCR_XL_flowSet())
    }
  )
}
# chnl_list.R
get_chnl_list <- function(fs) {
  set.seed(123)
  batch_list <- lapply(1:8, function(i) seq((i - 1) * 2 + 1, i * 2))
  chnl_vec <- c("BC1(La139)Dd", "BC2(Pr141)Dd")
  args_list_bc1  <- list(
    batch_list = batch_list,
    chnl = "BC1(La139)Dd",
    prop_mean_pos = 0.05,
    prop_sd_pos = 0.01,
    prop_mean_neg = 0.005,
    prop_sd_neg = 0.0075,
    expr_mean_neg = 0,
    expr_mean_pos = 2,
    expr_sd_pos = 0.1
  )
  args_list_bc2 <- list(
    batch_list = batch_list,
    chnl = "BC2(Pr141)Dd",
    prop_mean_pos = 0.01,
    prop_sd_pos = 0.01, 
    prop_mean_neg = 0.005,
    prop_sd_neg = 0.0075,
    expr_mean_neg = 0,
    expr_mean_pos = 2,
    expr_sd_pos = 0.1
  )
  args_list <- list(
    "BC1(La139)Dd" = args_list_bc1, "BC2(Pr141)Dd" = args_list_bc2
  )

  sample_chnls(args_list = args_list, fs = fs)
}

# sample.R
sample_chnls <- function(args_list,
                         fs) {
  chnl_obj_list <- lapply(seq_along(args_list), function(x) NULL) |>
    stats::setNames(names(args_list))
  for (i in seq_along(args_list)) {
    args <- args_list[[i]]
    # use updated flowSet
    args$fs <- fs
    chnl_obj_list[[i]] <- sample_chnl(
      fs = args$fs,
      batch_list = args$batch_list,
      chnl = args$chnl,
      prop_mean_pos = args$prop_mean_pos,
      prop_sd_pos = args$prop_sd_pos,
      prop_mean_neg = args$prop_mean_neg,
      prop_sd_neg = args$prop_sd_neg,
      expr_mean_neg = args$expr_mean_neg,
      expr_mean_pos = args$expr_mean_pos,
      expr_sd_pos = args$expr_sd_pos,
      expr_sd_neg = args$expr_sd_neg
    )
    # update flowSet
    fs <- chnl_obj_list[[i]]$fs
  }
  chnl_obj_list
}

sample_chnl <- function(fs,
                        batch_list,
                        chnl,
                        prop_mean_pos,
                        prop_sd_pos,
                        prop_mean_neg,
                        prop_sd_neg = NULL,
                        expr_mean_neg,
                        expr_mean_pos,
                        expr_sd_pos,
                        expr_sd_neg = NULL) {
  ind_list <- lapply(seq_along(fs), function(x) NULL)
  resp_tbl <- tibble::tibble(
    chnl = chnl,
    sample_ind = seq_along(fs),
    batch = rep(NA_integer_, length(fs)),
    n_cell = rep(NA_integer_, length(fs)),
    prop_pos = rep(NA_real_, length(fs)),
  )
  for (i in seq_along(batch_list)) {

    ind_uns <- batch_list[[i]][length(batch_list[[i]])]
    ind_stim_vec <- setdiff(batch_list[[i]], ind_uns)
    for (ind_stim in ind_stim_vec) {
      fr_stim <- fs[[ind_stim]]
      ex_mat_stim <- flowCore::exprs(fr_stim)
      n_cell_stim <- nrow(ex_mat_stim)
      response_stim <- sample_response(
        n = n_cell_stim,
        prop_mean = prop_mean_pos,
        prop_sd = prop_sd_pos,
        expr_mean_neg = expr_mean_neg,
        expr_mean_pos = expr_mean_pos,
        expr_sd_pos = expr_sd_pos,
        expr_sd_neg = expr_sd_neg
      )
      ex_mat_stim[, chnl] <- response_stim$expr_vec
      flowCore::exprs(fr_stim) <- ex_mat_stim
      fs[[ind_stim]] <- fr_stim
      resp_tbl$batch[ind_stim] <- i
      resp_tbl$n_cell[ind_stim] <- nrow(ex_mat_stim)
      resp_tbl$prop_pos[ind_stim] <- response_stim$n_pos / n_cell_stim
      ind_list[[ind_stim]] <- response_stim$ind_pos
    }
    fr_uns <- fs[[ind_uns]]
    ex_mat_uns <- flowCore::exprs(fr_uns)
    n_cell_uns <- nrow(ex_mat_uns)
    prop_sd_neg <- if (is.null(prop_sd_neg)) prop_sd_pos else prop_sd_neg
    response_uns <- sample_response(
      n = n_cell_uns,
      prop_mean = prop_mean_neg,
      prop_sd = prop_sd_neg,
      expr_mean_neg = expr_mean_neg,
      expr_mean_pos = expr_mean_pos,
      expr_sd_pos = expr_sd_pos,
      expr_sd_neg = expr_sd_neg   
    )
    ex_mat_uns[, chnl] <- response_uns$expr_vec
    flowCore::exprs(fr_uns) <- ex_mat_uns
    fs[[ind_uns]] <- fr_uns
    resp_tbl$batch[ind_uns] <- i
    resp_tbl$n_cell[ind_uns] <- nrow(ex_mat_uns)
    resp_tbl$prop_pos[ind_uns] <- response_uns$n_pos / n_cell_uns
    ind_list[[ind_uns]] <- response_uns$ind_pos
  }
  list(
    fs = fs,
    ind_list = ind_list,
    resp_tbl = resp_tbl,
    batch_list = batch_list
  )
}



sample_response <- function(n,
                            prop_mean,
                            prop_sd,
                            expr_mean_neg,
                            expr_mean_pos,
                            expr_sd_pos,
                            expr_sd_neg = NULL) {
  n_pos <- sample_n_pos(
    n = n,
    prop_mean = prop_mean,
    prop_sd = prop_sd
  )
  ind_pos <- sample_ind_pos(n = n, n_cell_pos = n_pos)
  expr_vec <- sample_expr(
    n = n,
    n_pos = n_pos,
    ind_pos = ind_pos,
    mean_neg = expr_mean_neg,
    mean_pos = expr_mean_pos,
    sd_pos = expr_sd_pos,
    sd_neg = expr_sd_neg
  )
  list(
    n_pos = n_pos,
    ind_pos = ind_pos,
    expr_vec = expr_vec
  )
}
sample_n_pos <- function(n, prop_mean, prop_sd, eps = 1e-8) {
  # Compute the maximum SD and force prop_sd < max_sd
  max_sd <- sqrt(prop_mean * (1 - prop_mean))
  prop_sd <- min(prop_sd, max_sd * (1 - eps))
  
  nu    <- prop_mean * (1 - prop_mean) / prop_sd^2 - 1
  alpha <- prop_mean * nu
  beta  <- (1 - prop_mean) * nu
  
  round(n * rbeta(n = 1, shape1 = alpha, shape2 = beta))
}

sample_ind_pos <- function(n, n_cell_pos) {
  sample.int(n, n_cell_pos)
}
sample_expr <- function(n,
                        n_pos,
                        ind_pos,
                        mean_neg,
                        mean_pos,
                        sd_pos,
                        sd_neg = NULL) {
  n_neg <- n - n_pos
  expr_vec <- rep(NA_real_, n)
  # only simulate when group has members
  if (n_pos > 0) {
    expr_vec[ind_pos] <- rnorm(n_pos, mean = mean_pos, sd = sd_pos)
  }
  if (n_neg > 0) {
    if (is.null(sd_neg)) {
      sd_neg <- sd_pos
    }
    expr_vec[setdiff(seq_len(n), ind_pos)] <-
      rnorm(n_neg, mean = mean_neg, sd = sd_neg)
  }
  expr_vec
}

calc_resp_combn <- function(chnl_list) {
  combn_pos_list <- calc_resp_combn_get_combn_pos(chnl_list)
  resp_tbl <- lapply(combn_pos_list, function(chnl_pos) {
    calc_resp_combn_ind(chnl_pos, chnl_list)
  }) |>
    dplyr::bind_rows()
  resp_tbl |>
    dplyr::group_by(chnl, batch) |>
    dplyr::mutate(
      prop_bs = prop_pos - prop_pos[length(prop_pos)]
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      n_pos_chnl = stringr::str_count(chnl, stringr::fixed("~+~")),
    ) |>
    dplyr::select(
      n_pos_chnl, chnl, batch, sample_ind, n_cell, ind, prop_pos, prop_bs,
      everything()
    )
 }
calc_resp_combn_get_combn_pos <- function(chnl_list) {
  chnl_vec <- names(chnl_list)
  # get combinations
  combn_mat_list_pos_init <- lapply(seq_along(chnl_vec), function(n_pos) {
    utils::combn(chnl_vec, m = n_pos, simplify = FALSE)
  })
  # flatten it
  combn_mat_list_pos <- list()
  for (i in seq_along(combn_mat_list_pos_init)) {
    combn_mat_list_pos_init_curr <- combn_mat_list_pos_init[[i]]
    for (j in seq_along(combn_mat_list_pos_init_curr)) {
      combn_mat_list_pos <- combn_mat_list_pos |>
        append(combn_mat_list_pos_init_curr[j])
    }
  }
  combn_mat_list_pos
}
calc_resp_combn_ind <- function(chnl_pos, chnl_list) {
  chnl_neg <- setdiff(names(chnl_list), chnl_pos)
  chnl_pos_ind <- which(names(chnl_list) %in% chnl_pos)
  sign_ind <- rep("~-~", length(chnl_list))
  sign_ind[chnl_pos_ind] <- "~+~"
  combn_nm <- paste0(names(chnl_list), sign_ind, collapse = "")
  n_sample <- length(chnl_list[[chnl_pos[1]]]$ind_list)
  # get resp_tbl as a scaffold
  resp_tbl <- chnl_list[[chnl_pos[1]]]$resp_tbl
  resp_tbl$chnl <- combn_nm
  ind_list_pos <- lapply(seq_len(n_sample), function(x) NULL)
  for (i in seq_len(n_sample)) {
    # get indices for each channel for this sample
    ind_list_sample <- lapply(seq_along(chnl_list), function(chnl) {
      chnl_list[[chnl]]$ind_list[[i]]
    }) |>
      stats::setNames(names(chnl_list))
    # get indices for positive channels
    pos_vec_ind <- Reduce(intersect, lapply(chnl_pos, function(chnl) {
      ind_list_sample[[chnl]]
    }))
    # no need to remove negative indices if nothing was positive
    if (length(chnl_neg) > 0L && length(pos_vec_ind) > 0L) {
      # union as we remove if any negative indices
      # were present in the positive indices
      neg_vec_ind <- Reduce(union, lapply(chnl_neg, function(chnl) {
        ind_list_sample[[chnl]]
      }))
      # remove negative indices from positive indices
      pos_vec_ind <- setdiff(pos_vec_ind, neg_vec_ind)
    }
    ind_list_pos[[i]] <- pos_vec_ind
    resp_tbl$prop_pos[i] <- length(pos_vec_ind) / resp_tbl$n_cell[i]
  }
  resp_tbl$ind <- ind_list_pos
  resp_tbl
}

# gatingset.R
get_gatingset <- function(fs,
                          dir_cache) {
  frames_list <- lapply(seq_along(fs), function(i) fs[[i]])
  fs2 <- flowCore::flowSet(frames = frames_list)
  gs <- flowWorkspace::GatingSet(fs2)
  path_save <- file.path(dir_cache, "gs")
  flowWorkspace::save_gs(
    gs,
    path = path_save
  )
  path_save
}
