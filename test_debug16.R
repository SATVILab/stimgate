source("R/cytoUtils-cytokine_cutpoint.R")
source("R/openCyto-find_peaks_and_valleys.R")
x_small <- c(1, 2, 3)
res <- tryCatch(.cytokine_cutpoint(x_small, plot = FALSE, strict=FALSE), warning=function(w) "warning", error=function(e) "error")
print(res)
# why did it throw an error in the test when wrapped in suppressWarnings?
# Oh wait, `outFunc <- ifelse(strict, stop, warning)`
# Wait! In R, `ifelse(TRUE, stop, warning)` returns a list of functions!
# No, ifelse returns the same type as the first argument, evaluating the condition.
# But `stop` and `warning` are closures! ifelse might not work correctly with closures!
