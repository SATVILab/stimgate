# data-raw/create_test_fixture.R
#
# Generate the tiny deterministic cytometry fixture shipped with the package.
# Run this script deliberately when the fixture needs to be regenerated:
#
#   source("data-raw/create_test_fixture.R")
#
# The resulting FCS files and metadata are stored under
#   inst/extdata/stimgate_test_fixture/
# and are loaded at runtime by .getTestFixture().

set.seed(42L)

# ---- Parameters -------------------------------------------------------
N_SAMPLE <- 2L # biological samples
N_CONDITION <- 2L # conditions per sample (1 = unstim, 2 = stim)
N_CELL <- 200L # cells per flow frame
N_MARKER <- 2L # markers per frame

CHNL_VEC <- c("BC1(La139)Dd", "BC2(Pr141)Dd")
MARKER_VEC <- paste0("MarkerF", seq_len(N_MARKER))

# Mean expression per cluster [negNeg, negPos, posNeg, posPos]
MEAN_MAT <- matrix(
  c(0, 0, 0, 8, 8, 0, 8, 8),
  nrow = 4L,
  byrow = TRUE
)
PROB_UNS <- c(0.90, 0.04, 0.04, 0.02)
# Probability shifts under stimulation
PROB_STIM_SHIFT <- c(-0.25, 0.2, 0.04, 0.01)
SD_CLUSTER <- 1.0

# ---- Helper: simulate one flow frame -----------------------------------
.sim_frame <- function(n_cell, prob_vec, sd_val = SD_CLUSTER) {
  n_cluster <- length(prob_vec)
  cluster_assign <- sample.int(n_cluster, n_cell, replace = TRUE, prob = prob_vec)
  expr_mat <- do.call(rbind, lapply(seq_len(n_cell), function(i) {
    cl <- cluster_assign[i]
    MEAN_MAT[cl, ] + stats::rnorm(N_MARKER, sd = sd_val)
  }))
  colnames(expr_mat) <- CHNL_VEC
  expr_mat
}

# ---- Build flow frames -------------------------------------------------
prob_stim <- pmax(PROB_UNS + PROB_STIM_SHIFT, 0)
prob_stim <- prob_stim / sum(prob_stim)

ff_list <- vector("list", N_SAMPLE * N_CONDITION)
idx <- 1L
for (s in seq_len(N_SAMPLE)) {
  for (cond in seq_len(N_CONDITION)) {
    prob_use <- if (cond == 1L) PROB_UNS else prob_stim
    expr_mat <- .sim_frame(N_CELL, prob_use)

    # Build flowFrame parameters
    params_df <- data.frame(
      name = CHNL_VEC,
      desc = MARKER_VEC,
      range = rep(2^18, N_MARKER),
      minRange = rep(min(expr_mat) - 1, N_MARKER),
      maxRange = rep(max(expr_mat) + 1, N_MARKER),
      stringsAsFactors = FALSE
    )
    rownames(params_df) <- paste0("$P", seq_len(N_MARKER))
    metadata <- new(
      "AnnotatedDataFrame",
      data = params_df,
      varMetadata = data.frame(
        labelDescription = rep(NA_character_, 5L),
        row.names = colnames(params_df)
      )
    )
    ff <- flowCore::flowFrame(
      exprs = expr_mat,
      parameters = metadata
    )
    ff_list[[idx]] <- ff
    idx <- idx + 1L
  }
}

# ---- Write FCS files and metadata -------------------------------------
out_dir <- file.path("inst", "extdata", "stimgate_test_fixture")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

fcs_names <- character(length(ff_list))
for (i in seq_along(ff_list)) {
  fname <- sprintf("sample%02d_cond%02d.fcs", ((i - 1L) %/% N_CONDITION) + 1L, ((i - 1L) %% N_CONDITION) + 1L)
  fcs_path <- file.path(out_dir, fname)
  flowCore::write.FCS(ff_list[[i]], filename = fcs_path)
  fcs_names[i] <- fname
}

# batchList: each element is c(unstim_idx, stim_idx)
batch_list <- lapply(seq_len(N_SAMPLE), function(s) {
  unstim_idx <- (s - 1L) * N_CONDITION + 1L
  stim_idx <- (s - 1L) * N_CONDITION + 2L
  c(unstim_idx, stim_idx)
})

metadata <- list(
  chnl = CHNL_VEC,
  marker = MARKER_VEC,
  batchList = batch_list,
  fcsNames = fcs_names,
  sampleNames = vapply(seq_along(fcs_names), function(i) {
    s <- ((i - 1L) %/% N_CONDITION) + 1L
    cond_suffix <- if ((i - 1L) %% N_CONDITION == 0L) "unstim" else "stim1"
    sprintf("sample%03d_%s", s, cond_suffix)
  }, character(1L))
)
saveRDS(metadata, file = file.path(out_dir, "metadata.rds"))

message("Fixture written to: ", out_dir)
message("Files: ", paste(c(fcs_names, "metadata.rds"), collapse = ", "))
