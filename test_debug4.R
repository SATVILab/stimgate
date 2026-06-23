source("R/cytoUtils-cytokine_cutpoint.R")
source("R/openCyto-find_peaks_and_valleys.R")
set.seed(123)
x <- c(rnorm(100, mean = 2, sd = 0.5))

peaks <- sort(.find_peaks(x, num_peaks = 1, adjust = 1.5)[, "x"])
print(peaks)

deriv_out <- .deriv_density(x = x, adjust = 1.5, deriv = 1)
# Original code: deriv_valleys <- with(deriv_out, .find_valleys(x = x, adjust = adjust))
# with(deriv_out, ...) means it uses deriv_out$x, which is the eval.points grid, not the original data!
# But wait, .find_valleys uses density(x). So it's calculating density of the grid points. That's a uniform distribution.
# Wait, look at cytoUtils original cytokine_cutpoint
