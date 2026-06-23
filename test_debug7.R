source("R/cytoUtils-cytokine_cutpoint.R")
source("R/openCyto-find_peaks_and_valleys.R")
set.seed(123)
x <- c(rnorm(100, mean = 2, sd = 0.5))
peaks <- sort(.find_peaks(x, num_peaks = 1, adjust = 1.5)[, "x"])

deriv_out <- .deriv_density(x = x, adjust = 1.5, deriv = 1)
# What cytoUtils actually did in cytokine_cutpoint:
# deriv_valleys <- .find_valleys(deriv_out$y) ? No... wait, .find_valleys needs the original data x?
# If we pass deriv_out$y to .find_valleys... wait, .find_valleys calculates density() of its input.
# The original openCyto flowClust/openCyto .find_valleys might have taken a density object or vector?

# In openCyto, find_valleys takes a density object or a numeric vector?
# Looking at original openCyto:
# https://github.com/RGLab/openCyto/blob/trunk/R/gating-functions.R
