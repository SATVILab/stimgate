source("R/cytoUtils-cytokine_cutpoint.R")
source("R/openCyto-find_peaks_and_valleys.R")
x_small <- c(1, 2, 3)
res <- .cytokine_cutpoint(x_small, plot = FALSE)
