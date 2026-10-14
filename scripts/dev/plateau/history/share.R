suppressMessages(library(dplyr)); P <- Sys.getenv("P"); options(width = 200, dplyr.summarise.inform = FALSE)
f1 <- function(x, l, g) { tp <- sum(x > g & l == "gp"); fp <- sum(x > g & l == "gn"); fn <- sum(x <= g & l == "gp"); 2*tp/(2*tp+fp+fn) }
runs <- c(am_g0001 = "Gaussian low 0.01% 100k (step off)", pk2_g0001 = "Gaussian low 0.01% 100k (step on)", am_g0005_20k = "Gaussian low 0.05% 20k (step off)",
  am_s0005 = "Skew high 0.05% 100k (step off)", pk2_s0005 = "Skew high 0.05% 100k (step on)", am_g002_20k = "Gaussian low 0.2% 20k (step off)", pk2_g002_20k = "Gaussian low 0.2% 20k (step on)")
out <- list()
for (t in names(runs)) { fd <- file.path(P, paste0("dbg_", t, ".rds")); fr <- file.path(P, paste0("res_", t, ".rds")); if (!file.exists(fd) || !file.exists(fr)) next
  d <- readRDS(fd); n <- length(d); if (grepl("^pk2", t) && n %% 2 == 0 && n >= 40) d <- d[(n / 2 + 1):n]
  res <- readRDS(fr) |> filter(method == "stimgate")
  for (e in d) { ind <- as.character(e$ind); own <- e$threshold; fin <- res$threshold[as.character(res$ind) == ind][1]; s <- e$sim[[1]]
    if (!is.finite(own) || is.null(fin) || !is.finite(fin)) next
    out[[length(out) + 1]] <- tibble(run = runs[[t]], ind, own, fin, f1Own = f1(s$stim, s$labelsStim, own), f1Fin = f1(s$stim, s$labelsStim, fin)) } }
r <- bind_rows(out)
print(as.data.frame(r |> group_by(run) |> summarise(tubes = n(), sharedLower = sum(fin < own - 1e-9), sharedHigher = sum(fin > own + 1e-9), unchanged = sum(abs(fin - own) <= 1e-9),
  medMoveDown = median((fin - own)[fin < own - 1e-9]), medMoveUp = median((fin - own)[fin > own + 1e-9]),
  f1Own = median(f1Own), f1Final = median(f1Fin), worseBySharing = sum(f1Fin < f1Own - 1e-9), betterBySharing = sum(f1Fin > f1Own + 1e-9))), digits = 3)
