suppressMessages(devtools::load_all(quiet = TRUE))
ex <- getExampleData(); str(ex, max.level = 1)
run <- function(on) { p <- tempfile("plateau"); dir.create(p)
  args <- ex; args$pathProject <- p; args$control <- stimControl(locPlateauStep = on)
  suppressMessages(do.call(gateStim, args[intersect(names(args), names(formals(gateStim)))])); getStimGates(p) }
a <- run(FALSE); b <- run(TRUE)
print(names(b)); print(head(as.data.frame(b), 4))
m <- merge(a, b, by = intersect(c("pop", "marker", "chnl", "gateName", "ind", "batch"), intersect(names(a), names(b))), suffixes = c(".off", ".on"))
print(m[, intersect(names(m), c("marker", "ind", "gate.off", "gate.on", "locPlateauStep"))])
