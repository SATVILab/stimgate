source("R/cytoUtils-cytokine_cutpoint.R")
source("R/openCyto-find_peaks_and_valleys.R")
set.seed(123)
x <- c(rnorm(100, mean = 2, sd = 0.5))
peaks <- sort(.find_peaks(x, num_peaks = 1, adjust = 1.5)[, "x"])
print(peaks)
deriv_out <- .deriv_density(x = x, adjust = 1.5, deriv = 1)
# Note: we need to find valleys of deriv_out$x and deriv_out$y!
# .find_valleys does KDE, so giving it deriv_out$x just finds valleys of the distribution of grid points...
# .find_valleys uses `density(x)`, so it's expecting the raw data points!
deriv_valleys_x <- .find_valleys(x = x, adjust = 1.5)
print(deriv_valleys_x)
