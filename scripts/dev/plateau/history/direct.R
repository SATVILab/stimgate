suppressMessages(devtools::load_all(quiet = TRUE)); P <- Sys.getenv("P")
A <- readRDS(file.path(P, "dbg_am_g002_20k.rds"))
mats <- list(); bl <- list()
for (i in seq_along(A)) { s <- A[[i]]$sim[[1]]
  mats[[paste0("u", i)]] <- matrix(s$uns, ncol = 1, dimnames = list(NULL, "F1"))
  mats[[paste0("s", i)]] <- matrix(s$stim, ncol = 1, dimnames = list(NULL, "F1"))
  bl[[paste0("b", i)]] <- c(paste0("u", i), paste0("s", i)) }
run <- function(on) { p <- tempfile("direct"); dir.create(p); set.seed(1); Sys.setenv(STIMGATE_DEBUG = "true", STIMGATE_INTERMEDIATE = "all")
  suppressMessages(gateStim(pathProject = p, .data = mats, batchList = bl, chnl = "F1", calcCytPosGates = FALSE, control = stimControl(locPlateauStep = on)))
  list(g = getStimGates(p), p = p) }
off <- run(FALSE); on <- run(TRUE)
a <- off$g[off$g$gateName == "loc_min", ]; b <- on$g[on$g$gateName == "loc_min", ]
cat("lowered:", sum(b$gate < a$gate - 1e-9), " raised:", sum(b$gate > a$gate + 1e-9), "\n"); print(table(b$locPlateauStep))
f <- list.files(on$p, pattern = "locPlateau", recursive = TRUE, full.names = TRUE); cat(length(f), "plateau files\n"); print(head(basename(f)))
al <- f[grepl("Align", f)]; if (length(al)) { x <- readRDS(al[1]); print(unlist(x$group)); print(head(as.data.frame(x$tbl))) }
st <- f[grepl("locPlateauStep", f)]; if (length(st)) print(table(sapply(st, function(z) readRDS(z)$reason)))
