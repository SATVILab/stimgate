#' Get example GatingSet
#'
#' Create and save a complete example GatingSet for testing and examples.
#' This function internally creates a flowSet, samples channels, and saves the GatingSet.
#'
#' @param dir_cache Directory to save the GatingSet. If NULL, uses a temporary directory.
#' @return A list containing the path to the saved GatingSet, batch_list, and marker names
#' @export
get_test_data <- function(
  scenario = "default",
  dir_cache = NULL,
  clear = FALSE,
  n_ind = 8
) {
  # Set seed for reproducibility
  set.seed(123)

  if (is.null(dir_cache)) {
    dir_cache <- testthat::test_path("cache", "test_data", scenario)
  }
  if (dir.exists(dir_cache)) {
    if (clear) {
      unlink(dir_cache, recursive = TRUE)
    } else {
      cache_list <- try(.get_test_data_cache(dir_cache), silent = TRUE)
      if (!inherits(cache_list, "try-error")) {
        return(cache_list)
      } else {
        message("Cache not found or incomplete, regenerating test data...")
        unlink(dir_cache, recursive = TRUE)
      }
    }
    dir.create(dir_cache, recursive = TRUE)
  } else {
    dir.create(dir_cache, recursive = TRUE)
  }

  .get_test_data_fresh(
    scenario = scenario,
    dir_cache = dir_cache,
    n_ind = n_ind
  )
}

.get_test_data_cache <- function(dir_cache) {
  path_gs <- file.path(dir_cache, "gs")
  if (!dir.exists(path_gs)) {
    stop("No cached test data found at: ", dir_cache)
  }
  path_chnl <- file.path(dir_cache, "chnl.rds")
  path_marker <- file.path(dir_cache, "marker.rds")
  path_batch_list <- file.path(dir_cache, "batch_list.rds")
  if (
    !file.exists(path_chnl) ||
      !file.exists(path_marker) ||
      !file.exists(path_batch_list)
  ) {
    stop("Incomplete cached test data at: ", dir_cache)
  }
  list(
    path_gs = path_gs,
    chnl = readRDS(path_chnl),
    marker = readRDS(path_marker),
    batch_list = readRDS(path_batch_list)
  )
}

.get_test_data_fresh <- function(scenario, dir_cache, n_ind) {
  fs <- .get_fs_test(dir_cache)

  # Get flowSet
  fs_list_sim <- switch(
    scenario,
    "default" = sim_fs_default(fs, n_ind),
    "easy" = sim_fs_easy(fs, n_ind),
    "poor_separation" = sim_fs_poor_separation(fs, n_ind),
    "cyt_pos" = sim_fs_cyt_pos(fs, n_ind),
    stop("Unknown scenario: ", scenario)
  )

  # Extract the final flowSet from the last channel in the list
  fs_sim <- fs_list_sim[[length(fs_list_sim)]]$fs

  # Create and save GatingSet
  path_gs <- get_gatingset(
    fs = fs_sim,
    dir_cache = dir_cache
  )

  desc_df <- flowCore::parameters(fs_list_sim[[1]][[1]][[1]])@data
  chnl_vec <- names(fs_list_sim)
  saveRDS(chnl_vec, file = file.path(dir_cache, "chnl.rds"))
  marker_vec <- NULL
  for (chnl_curr in chnl_vec) {
    marker_vec <- c(
      marker_vec,
      desc_df$desc[which(desc_df$name == chnl_curr)]
    )
  }
  saveRDS(marker_vec, file = file.path(dir_cache, "marker.rds"))
  batch_list <- fs_list_sim[[1]]$batch_list
  saveRDS(batch_list, file = file.path(dir_cache, "batch_list.rds"))

  list(
    path_gs = path_gs,
    batch_list = batch_list,
    chnl = chnl_vec,
    marker = marker_vec
  )
}

# Internal helper functions (not exported)

#' @keywords internal
.get_fs_test <- function(dir_cache) {
  dir_fs <- file.path(dir_cache, "fs")
  if (!dir.exists(dir_fs)) {
    dir.create(dir_fs, recursive = TRUE)
  }
  # Try to load from installed package location first
  rds_path <- file.path(dir_fs, "bodenmiller_bcr_xl_fs.rds")
  if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  }

  # If not found, try to download from GitHub release
  message(
    "Test data not found locally. Attempting to download from GitHub release..."
  )

  tryCatch(
    {
      .download_fs_from_github(dir_fs)
    },
    error = function(e_github) {
      message("Failed to download from GitHub: ", conditionMessage(e_github))
      message("Falling back to HDCytoData...")

      # Fall back to HDCytoData
      tryCatch(
        {
          .download_fs_from_hdcytodata()
        },
        error = function(e_hdc) {
          stop(
            "Failed to obtain test data from all sources.\n",
            "GitHub error: ",
            conditionMessage(e_github),
            "\n",
            "HDCytoData error: ",
            conditionMessage(e_hdc)
          )
        }
      )
    }
  )
}

#' @keywords internal
.download_fs_from_github <- function(dir_fs) {
  repo <- "SATVILab/stimgate"
  tag <- "test_data"
  filename <- "bodenmiller_bcr_xl_fs.rds"

  # Build download URL (no authentication needed for public releases)
  download_url <- sprintf(
    "https://github.com/%s/releases/download/%s/%s",
    repo,
    tag,
    filename
  )

  temp_file <- file.path(dir_fs, filename)
  message("Downloading ", filename, " from GitHub release...")

  result <- tryCatch(
    {
      utils::download.file(
        url = download_url,
        destfile = temp_file,
        mode = "wb",
        quiet = FALSE
      )
      TRUE
    },
    error = function(e) {
      message("Download failed: ", conditionMessage(e))
      FALSE
    }
  )

  if (!result || !file.exists(temp_file)) {
    stop("Failed to download test data from GitHub")
  }

  # Read and return
  fs <- readRDS(temp_file)
  unlink(temp_file)

  message("Successfully downloaded test data from GitHub")
  fs
}

#' @keywords internal
.download_fs_from_hdcytodata <- function() {
  path_hub <- file.path(Sys.getenv("HOME"), ".cache/R/ExperimentHub")
  if (!dir.exists(path_hub)) {
    dir.create(path_hub, recursive = TRUE)
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    if (interactive()) {
      prompt_answer <- readline(
        prompt = paste0(
          "The 'BiocManager' package is required to download the example data. ",
          "Do you want to install it now? [y/n]: "
        )
      )
      if (tolower(prompt_answer) != "y") {
        stop("Cannot proceed without installing 'BiocManager' package.")
      }
    } else {
      stop("Cannot proceed without installing 'BiocManager' package.")
    }
    utils::install.packages("BiocManager")
  }

  if (!requireNamespace("HDCytoData", quietly = TRUE)) {
    if (interactive()) {
      prompt_answer <- readline(
        prompt = paste0(
          "The 'HDCytoData' package is required to download the example data. ",
          "Do you want to install it now? [y/n]: "
        )
      )
      if (tolower(prompt_answer) != "y") {
        stop("Cannot proceed without installing 'HDCytoData' package.")
      }
      BiocManager::install("HDCytoData")
    } else {
      stop("Cannot proceed without installing 'HDCytoData' package.")
    }
  }

  # Get an odd flowSet slot dropping warning,
  # which is unrelated to our functionality
  suppressWarnings(HDCytoData::Bodenmiller_BCR_XL_flowSet())
}

#' @keywords internal
.get_chnl_list <- function(fs) {
  batch_list <- lapply(1:8, function(i) seq((i - 1) * 2 + 1, i * 2))
  chnl_vec <- c("BC1(La139)Dd", "BC2(Pr141)Dd")
  args_list_bc1 <- list(
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
    "BC1(La139)Dd" = args_list_bc1,
    "BC2(Pr141)Dd" = args_list_bc2
  )

  .sample_chnls(args_list = args_list, fs = fs)
}

# Internal helper functions (not exported)

#' @keywords internal
.sample_chnls <- function(args_list, fs) {
  chnl_obj_list <- lapply(seq_along(args_list), function(x) NULL) |>
    stats::setNames(names(args_list))
  for (i in seq_along(args_list)) {
    args <- args_list[[i]]
    # use updated flowSet
    args$fs <- fs
    chnl_obj_list[[i]] <- .sample_chnl(
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

#' @keywords internal
.sample_chnl <- function(
  fs,
  batch_list,
  chnl,
  prop_mean_pos,
  prop_sd_pos,
  prop_mean_neg,
  prop_sd_neg = NULL,
  expr_mean_neg,
  expr_mean_pos,
  expr_sd_pos,
  expr_sd_neg = NULL
) {
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
      response_stim <- .sample_response(
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
    response_uns <- .sample_response(
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

#' @keywords internal
.sample_response <- function(
  n,
  prop_mean,
  prop_sd,
  expr_mean_neg,
  expr_mean_pos,
  expr_sd_pos,
  expr_sd_neg = NULL
) {
  n_pos <- .sample_n_pos(
    n = n,
    prop_mean = prop_mean,
    prop_sd = prop_sd
  )
  ind_pos <- .sample_ind_pos(n = n, n_cell_pos = n_pos)
  expr_vec <- .sample_expr(
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

#' @keywords internal
.sample_n_pos <- function(n, prop_mean, prop_sd, eps = 1e-8) {
  # Compute the maximum SD and force prop_sd < max_sd
  max_sd <- sqrt(prop_mean * (1 - prop_mean))
  prop_sd <- min(prop_sd, max_sd * (1 - eps))

  nu <- prop_mean * (1 - prop_mean) / prop_sd^2 - 1
  alpha <- prop_mean * nu
  beta <- (1 - prop_mean) * nu

  round(n * stats::rbeta(n = 1, shape1 = alpha, shape2 = beta))
}

#' @keywords internal
.sample_ind_pos <- function(n, n_cell_pos) {
  sample.int(n, n_cell_pos)
}

#' @keywords internal
.sample_expr <- function(
  n,
  n_pos,
  ind_pos,
  mean_neg,
  mean_pos,
  sd_pos,
  sd_neg = NULL
) {
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
