source("R/cytoUtils-cytokine_cutpoint.R")
source("R/openCyto-find_peaks_and_valleys.R")
set.seed(123)
x <- c(rnorm(100, mean = 2, sd = 0.5))

peaks <- sort(.find_peaks(x, num_peaks = 1, adjust = 1.5)[, "x"])
print(peaks)

deriv_out <- .deriv_density(x = x, adjust = 1.5, deriv = 1)
# what if we use the same diff sign logic to find valleys of deriv_out$y directly?
sd2 <- diff(sign(diff(deriv_out$y)))
idx <- which(sd2 == 2) + 1L
deriv_valleys <- deriv_out$x[idx]
print(deriv_valleys)
deriv_valleys <- deriv_valleys[deriv_valleys > peaks[1]]
deriv_valleys <- sort(deriv_valleys)[1]

tol <- 0.01 * max(abs(deriv_out$y))
cutpoint <- with(deriv_out, x[x > deriv_valleys & abs(y) < tol])
cutpoint <- cutpoint[1]
print(cutpoint)
