# data-raw/create_test_fixture.R
#
# Generate the canonical saved example dataset shipped with the package.
# Run this script deliberately when the fixture needs to be regenerated:
#
#   source("data-raw/create_test_fixture.R")
#
# The resulting FCS files and metadata are stored under
#   inst/extdata/stimgate_example_data/
# and are loaded at runtime by getExampleData().

source(file.path("scripts", "r", "functionsForBenchmarking-Cyt.R"))
set.seed(42L)

# ---- Parameters -------------------------------------------------------
N_SAMPLE <- 2L
N_CONDITION <- 2L
N_CELL <- 1e4L
N_MARKER <- 2L

CHNL_VEC <- c("BC1(La139)Dd", "BC2(Pr141)Dd")
MARKER_VEC <- paste0("MarkerF", seq_len(N_MARKER))
CLUSTER_LABEL_VEC <- c("negNeg", "negPos", "posNeg", "posPos")
MEAN_MAT <- matrix(
  c(0, 0, 0, 4, 4, 0, 4, 4),
  nrow = 4L,
  byrow = TRUE
)

# Use the canonical response model required by the package example dataset.
# Marker 1 responds on stimulated cells by +0.02; marker 2 by +0.05.
PROB_UNS <- c(0.988, 0.008, 0.002, 0.002)
PROB_STIM <- c(0.926, 0.050, 0.014, 0.010)
PROB_RESPONSE_VEC <- list(PROB_STIM - PROB_UNS)

# ---- Build canonical dataset -------------------------------------------
sim_res <- simCytExperiment(
  nSample = N_SAMPLE,
  nMarker = N_MARKER,
  nCondition = N_CONDITION,
  nCluster = length(CLUSTER_LABEL_VEC),
  nCellByCondition = N_CELL,
  transformationFunc = function(x) x,
  mixtureType = "gaussianOnly",
  meanExprMat = MEAN_MAT,
  clusterLabelVec = CLUSTER_LABEL_VEC,
  probVecUns = PROB_UNS,
  probExact = TRUE,
  probResponseVecByStimCondition = PROB_RESPONSE_VEC,
  clusterPerturbationSd = 0,
  conditionPerturbationSd = 0,
  samplePerturbationSd = 0
)

ff_list <- lapply(sim_res$flowFrameList, function(ff) {
  flowCore::colnames(ff) <- CHNL_VEC
  p <- flowCore::parameters(ff)
  p_data <- flowCore::pData(p)
  p_data$name <- CHNL_VEC
  p_data$desc <- MARKER_VEC
  flowCore::pData(p) <- p_data
  ff
})

# ---- Write FCS files and metadata -------------------------------------
out_dir <- file.path("inst", "extdata", "stimgate_example_data")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

fcs_names <- character(length(ff_list))
for (i in seq_along(ff_list)) {
  fname <- sprintf(
    "sample%02d_cond%02d.fcs",
    ((i - 1L) %/% N_CONDITION) + 1L,
    ((i - 1L) %% N_CONDITION) + 1L
  )
  fcs_path <- file.path(out_dir, fname)
  flowCore::write.FCS(ff_list[[i]], filename = fcs_path)
  fcs_names[i] <- fname
}

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

message("Canonical example data written to: ", out_dir)
message("Files: ", paste(c(fcs_names, "metadata.rds"), collapse = ", "))
