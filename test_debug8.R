source("R/cytoUtils-cytokine_cutpoint.R")
source("R/openCyto-find_peaks_and_valleys.R")
set.seed(123)
x <- c(rnorm(100, mean = 2, sd = 0.5))
deriv_out <- .deriv_density(x = x, adjust = 1.5, deriv = 1)
str(deriv_out)

# We want the valleys of the first derivative.
# The first derivative is in `deriv_out$y`, evaluated at `deriv_out$x`.
# Local minima of `deriv_out$y` are valleys of the first derivative!
sd2 <- diff(sign(diff(deriv_out$y)))
idx <- which(sd2 == 2) + 1L
deriv_valleys <- deriv_out$x[idx]
print(deriv_valleys)
