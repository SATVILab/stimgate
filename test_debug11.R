source("R/cytoUtils-cytokine_cutpoint.R")
source("R/openCyto-find_peaks_and_valleys.R")
set.seed(123)
x <- c(rnorm(100, mean = 2, sd = 0.5))

peaks <- sort(.find_peaks(x, num_peaks = 1, adjust = 1.5)[, "x"])
print(peaks)

deriv_out <- .deriv_density(x = x, adjust = 1.5, deriv = 1)
# Look at original cytoUtils code:
# deriv_valleys <- with(deriv_out, .find_valleys(x = x, adjust = adjust))
# BUT deriv_out is a list with x and y.
# If .find_valleys originally accepted a DENSITY object or similar?
# If we replace `with(deriv_out, .find_valleys(x = x, adjust = adjust))` with finding valleys of the derivative:
