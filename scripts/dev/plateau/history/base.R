suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
src <- readLines(file.path(P, "cmpab.R")); eval(parse(text = src[2:15]))
# Step with unstim counts scaled by a baseline stim/unstim ratio measured well below the gate.
stepR <- function(xs, xu, g, fl, bw, r0) {
  ratio <- r0 * length(xs) / length(xu); nSr <- sum(xs > g); refPur <- (nSr - sum(xu > g) * ratio) / nSr
  if (!is.finite(fl) || fl >= g || !is.finite(refPur) || refPur <= 0) return(g)
  cur <- g; hi <- g; rej <- 0
  while (rej < 2) { lo <- max(hi - bw, fl); if (lo >= hi - 1e-9) break
    nS <- sum(xs > lo & xs <= cur); nU <- sum(xu > lo & xu <= cur) * ratio
    if (nS > 0 && (nS - nU) / nS >= 0.5 * refPur) { cur <- lo; rej <- 0 } else rej <- rej + 1
    hi <- lo; if (lo <= fl) break }
  cur }
scen <- c(g002 = "Clean low 0.2% 100k", g002_20k = "Clean low 0.2% 20k", g0005_20k = "Clean low 0.05% 20k", g0001 = "Clean low 0.01% 100k",
  h002 = "Clean high 0.2% 100k", h002_20k = "Clean high 0.2% 20k", h0005_20k = "Clean high 0.05% 20k", h0001 = "Clean high 0.01% 100k", s002 = "Skew high 0.2% 100k",
  s0005 = "Skew high 0.05% 100k", s002_20k = "Skew high 0.2% 20k", s0005_20k = "Skew high 0.05% 20k",
  g002_ms02 = "Shift low 0.2% 100k", g002_20k_ms02 = "Shift low 0.2% 20k", h002_20k_ms02 = "Shift high 0.2% 20k",
  g002_sd01 = "SD low 0.2% 100k", g002_20k_sd01 = "SD low 0.2% 20k", h002_20k_sd01 = "SD high 0.2% 20k")
out <- list()
for (t in names(scen)) { A <- readRDS(file.path(P, paste0("dbg_am_", t, ".rds"))); B <- readRDS(file.path(P, paste0("dbg_pb_", t, ".rds")))
  for (i in seq_along(A)) { a <- A[[i]]; s <- a$sim[[1]]; xs <- s$stim; xu <- s$uns; l <- s$labelsStim
    g0 <- a$threshold; lb <- a$filter[[1]]$info$minProbXPos; bw <- B[[i]]$filter[[1]]$info$antimode$bw %||% NA
    if (!is.finite(g0) || !is.finite(bw) || is.null(lb)) next
    fl <- edge(sort(xs[xs >= lb]), g0)
    # baseline band: [fl - 2 bw - W, fl - 2 bw] would be scale-free; here use a band ending 1.5 units below the gate, 0.5 wide
    hiB <- g0 - 1.5; loB <- g0 - 2
    r0 <- (sum(xs > loB & xs <= hiB) / length(xs)) / (sum(xu > loB & xu <= hiB) / length(xu))
    rUse <- max(1, r0)
    truth <- mean(l == "gp") - mean(s$labelsUns == "gp")
    gS <- stepDown(xs, xu, g0, fl, bw); gR <- stepR(xs, xu, g0, fl, bw, rUse)
    out[[length(out) + 1]] <- tibble(scenario = scen[[t]], i, r0, f1Now = f1(xs, l, g0), f1Step = f1(xs, l, gS), f1Ratio = f1(xs, l, gR),
      eNow = bsf(xs, xu, g0) / truth - 1, eStep = bsf(xs, xu, gS) / truth - 1, eRatio = bsf(xs, xu, gR) / truth - 1) } }
r <- bind_rows(out); options(width = 220)
q <- function(v) sprintf("%+.0f [%+.0f,%+.0f]", 100 * median(v), 100 * quantile(v, .1), 100 * quantile(v, .9))
print(r |> group_by(scenario) |> summarise(r0 = median(r0), f1Now = median(f1Now), f1Step = median(f1Step), f1Ratio = median(f1Ratio),
  errNow = q(eNow), errStep = q(eStep), errRatio = q(eRatio)) |> as.data.frame(), digits = 3)
