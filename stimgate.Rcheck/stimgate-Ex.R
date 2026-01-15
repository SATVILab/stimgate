pkgname <- "stimgate"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "stimgate-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('stimgate')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("chnl_lab")
### * chnl_lab

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: chnl_lab
### Title: Get markers and channels
### Aliases: chnl_lab marker_lab, chnl_to_marker, marker_to_chnl,
###   get_marker, get_chnl

### ** Examples

## No test: 
# Create example flowFrame-like data structure
data(GvHD, package = "flowCore") 
fs <- GvHD[1:2]

# Get channel to marker mapping
chnl_lab(fs[[1]])
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("chnl_lab", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_stats")
### * get_stats

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_stats
### Title: Get gating statistics
### Aliases: get_stats

### ** Examples

{
# Get example dataset
example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)

# Run the stimgate pipeline
path_project <- stimgate_gate(
  path_project = file.path(tempdir(), "stimgate_example"),
  .data = gs,
  batch_list = example_data$batch_list,
  marker = example_data$marker,
  pop_gate = "root"
)

# Get statistics for the identified gates
stats <- get_stats(path_project)
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_stats", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_data_get_ex")
### * stimgate_data_get_ex

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_data_get_ex
### Title: Read saved expression data from project
### Aliases: stimgate_data_get_ex

### ** Examples

## Not run: 
##D tmp <- tempdir()
##D dir.create(file.path(tmp, "sample_data", "POP1", "ind_1"), recursive = TRUE)
##D saveRDS(c(1, 2, 3), file.path(tmp, "sample_data", "POP1", "ind_1", "chnl_BC1.rds"))
##D saveRDS(c(4, 5, 6), file.path(tmp, "sample_data", "POP1", "ind_1", "chnl_BC2.rds"))
##D stimgate_data_get_ex(tmp)
##D stimgate_data_get_ex(tmp, chnl = "BC1")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_data_get_ex", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_fcs_write")
### * stimgate_fcs_write

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_fcs_write
### Title: Write FCS files of marker-positive FCS files
### Aliases: stimgate_fcs_write

### ** Examples

## Not run: 
##D # Complete workflow example
##D library(stimgate)
##D 
##D # Load your GatingSet (gs) and define batch structure
##D # batch_list <- list(batch1 = c(1, 2, 3), batch2 = c(4, 5, 6))
##D # where the last element in each batch is the unstimulated sample
##D 
##D # First, run gating to create gates
##D path_project <- tempfile("stimgate_project")
##D # stimgate_gate(
##D #   .data = gs,
##D #   path_project = path_project,
##D #   pop_gate = "root",
##D #   batch_list = batch_list,
##D #   marker = c("IL2", "IFNg")  # your cytokine markers
##D # )
##D 
##D # Then write FCS files of cytokine-positive cells
##D path_output <- tempfile("fcs_output")
##D # stimgate_fcs_write(
##D #   path_project = path_project,
##D #   .data = gs,
##D #   ind_batch_list = batch_list,
##D #   path_dir_save = path_output,
##D #   chnl = c("IL2", "IFNg")
##D # )
##D 
##D # Alternative: provide your own gate table
##D # gate_tbl <- data.frame(
##D #   chnl = c("IL2", "IFNg"),
##D #   marker = c("IL2", "IFNg"),
##D #   batch = c(1, 1),
##D #   ind = c(1, 1),
##D #   gate = c(0.5, 0.3),
##D #   gate_name = c("gate", "gate")
##D # )
##D # stimgate_fcs_write(
##D #   path_project = path_project,
##D #   .data = gs,
##D #   ind_batch_list = batch_list,
##D #   path_dir_save = path_output,
##D #   chnl = c("IL2", "IFNg"),
##D #   gate_tbl = gate_tbl
##D # )
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_fcs_write", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_gate")
### * stimgate_gate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_gate
### Title: Identify cytokine-positive cells through automated gating
### Aliases: stimgate_gate

### ** Examples

{
example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)
path_project <- file.path(tempdir(), "demonstration")

# Run gating
stimgate::stimgate_gate(
  .data = gs,
  path_project = path_project,
  pop_gate = "root",
  batch_list = example_data$batch_list,
  marker = example_data$marker
)

# Create plots
plots <- stimgate_plot(
  ind = example_data$batch_list[[1]], # indices in `gs` to plot
  .data = gs, # GatingSet
  path_project = path_project,
  marker = example_data$marker,
  grid = TRUE
)

# Advanced usage with parameter customization

}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_gate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_gate_get")
### * stimgate_gate_get

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_gate_get
### Title: Get gates
### Aliases: stimgate_gate_get

### ** Examples

{
# Get example dataset
example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)

# Run the stimgate pipeline
path_project <- stimgate_gate(
  path_project = file.path(tempdir(), "get_gate_example"),
  .data = gs,
  batch_list = example_data$batch_list,
  marker = example_data$marker,
  pop_gate = "root"
)

# Get statistics for the identified gates
gates <- stimgate_gate_get(path_project)
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_gate_get", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_meta_read_batch_list")
### * stimgate_meta_read_batch_list

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_meta_read_batch_list
### Title: Read batch list from project
### Aliases: stimgate_meta_read_batch_list

### ** Examples

## Not run: 
##D tmp <- tempdir()
##D dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
##D saveRDS(list(batch1 = c(1,2)), file.path(tmp, "meta_data", "batch_list.rds"))
##D stimgate_meta_read_batch_list(tmp)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_meta_read_batch_list", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_meta_read_lab")
### * stimgate_meta_read_lab

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_meta_read_chnl_lab
### Title: Read channel or marker label mapping
### Aliases: stimgate_meta_read_chnl_lab stimgate_meta_read_marker_lab

### ** Examples

## Not run: 
##D tmp <- tempdir()
##D dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
##D saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "meta_data", "chnl_lab.rds"))
##D stimgate_meta_read_chnl_lab(tmp)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_meta_read_lab", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_meta_read_settings_chnl")
### * stimgate_meta_read_settings_chnl

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_meta_read_settings_chnl
### Title: Get marker settings for a single channel
### Aliases: stimgate_meta_read_settings_chnl

### ** Examples

## Not run: 
##D tmp <- tempdir()
##D dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
##D saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "meta_data", "marker_list.rds"))
##D saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "meta_data", "chnl_lab.rds"))
##D stimgate_meta_read_settings_chnl(tmp, "BC1 label")
##D stimgate_meta_read_settings_chnl(tmp, "BC1")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_meta_read_settings_chnl", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_meta_read_settings_chnls")
### * stimgate_meta_read_settings_chnls

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_meta_read_settings_chnls
### Title: Read marker settings from project
### Aliases: stimgate_meta_read_settings_chnls

### ** Examples

## Not run: 
##D tmp <- tempdir()
##D dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
##D saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "meta_data", "marker_list.rds"))
##D stimgate_meta_read_settings_markers(tmp)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_meta_read_settings_chnls", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_meta_read_settings_marker")
### * stimgate_meta_read_settings_marker

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_meta_read_settings_marker
### Title: Get settings for a named marker
### Aliases: stimgate_meta_read_settings_marker

### ** Examples

## Not run: 
##D tmp <- tempdir()
##D dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
##D saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "meta_data", "marker_list.rds"))
##D stimgate_meta_read_settings_marker(tmp, "BC1")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_meta_read_settings_marker", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_meta_read_settings_markers")
### * stimgate_meta_read_settings_markers

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_meta_read_settings_markers
### Title: Read marker list with channel labels
### Aliases: stimgate_meta_read_settings_markers

### ** Examples

## Not run: 
##D tmp <- tempdir()
##D dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
##D saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "meta_data", "marker_list.rds"))
##D saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "meta_data", "chnl_lab.rds"))
##D stimgate_meta_read_settings_chnls(tmp)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_meta_read_settings_markers", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimgate_plot")
### * stimgate_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimgate_plot
### Title: Plot stimulation gate
### Aliases: stimgate_plot

### ** Examples

# Create example data and run gating
example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)
path_project <- file.path(dirname(example_data$path_gs), "stimgate")

# Run gating
stimgate::stimgate_gate(
  .data = gs,
  path_project = path_project,
  pop_gate = "root",
  batch_list = example_data$batch_list,
  marker = example_data$marker
)

# Create plots
plots <- stimgate_plot(
  ind = example_data$batch_list[[1]], # indices in `gs` to plot
  .data = gs, # GatingSet
  path_project = path_project,
  marker = example_data$marker,
  grid = TRUE
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimgate_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
