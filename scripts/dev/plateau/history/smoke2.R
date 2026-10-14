suppressMessages(devtools::load_all(quiet = TRUE))
ex <- getExampleData(); p <- tempfile("plateau"); dir.create(p); Sys.setenv(STIMGATE_DEBUG = "true"); set.seed(1)
rec <- NULL
trace(".gateChnlPlateauAlign", exit = quote(rec <<- returnValue()), where = asNamespace("stimgate"), print = FALSE)
suppressMessages(gateStim(pathProject = p, .data = flowWorkspace::load_gs(ex$pathGs), batchList = ex$batchList, marker = ex$marker, calcCytPosGates = FALSE, control = stimControl(locPlateauStep = TRUE)))
cat("traced:", !is.null(rec), "\n"); print(rec$locPlateauAllowed)
d <- readLines(file.path(p, "debug", "debug.txt")); print(grep("plateau", d, value = TRUE, ignore.case = TRUE)[1:5])
