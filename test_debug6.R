source("R/cytoUtils-cytokine_cutpoint.R")
source("R/openCyto-find_peaks_and_valleys.R")
set.seed(123)
x <- c(rnorm(100, mean = 2, sd = 0.5))
deriv_out <- .deriv_density(x = x, adjust = 1.5, deriv = 1)
# Original openCyto find_valleys might have worked on the y vector directly
# .find_valleys uses stats::density on the x values. Wait...
# if with(deriv_out, .find_valleys(x = x, adjust=adjust)), it passes `deriv_out$x` as x! Which is uniformly spaced points.
