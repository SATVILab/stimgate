source("R/cytoUtils-cytokine_cutpoint.R")
source("R/openCyto-find_peaks_and_valleys.R")
set.seed(123)
x <- c(rnorm(100, mean = 2, sd = 0.5))

deriv_out <- .deriv_density(x = x, adjust = 1.5, deriv = 1)
str(deriv_out)
