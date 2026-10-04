# Getting Started with stimgate

``` r

library(stimgate)
```

## Introduction

The `stimgate` package provides tools to identify cells that have
possibly responded to stimulation by comparing unstimulated and
stimulated tubes from the same sample.

### Main Functions

The package provides several key functions:

- [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md):
  Main function to identify cytokine-positive cells by gating
- [`getStimStats()`](https://satvilab.github.io/stimgate/reference/getStimStats.md):
  Get statistics from gating results
- [`plotStim()`](https://satvilab.github.io/stimgate/reference/plotStim.md):
  Visualise identified gates
- [`getStimGates()`](https://satvilab.github.io/stimgate/reference/getStimGates.md):
  Extract gate information
- [`writeStimFCS()`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md):
  Write FCS files of cytokine-positive cells

### Input formats

[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md),
[`plotStim()`](https://satvilab.github.io/stimgate/reference/plotStim.md),
[`writeStimFCS()`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md)
and
[`getStimExpr()`](https://satvilab.github.io/stimgate/reference/getStimExpr.md)
accept a GatingSet, flowSet/cytoset, individual flowFrame/cytoframe, FCS
file paths or a single FCS directory, and a list of numeric matrices or
data frames (cells by channels). Lists use their names as sample names,
or `sample1`, `sample2`, etc. Columns must have matching names; their
order is aligned to the first sample. A long data frame may instead
contain a `sample` column: samples follow observed factor levels or
first appearance. Matrix column names work as channels and markers.

FCS directories are searched non-recursively, case-insensitively, with
sorted filenames. Explicit file vectors preserve order; sample names are
file basenames. Only GatingSets provide populations beyond `root`.
StimGate adds no arcsinh or logicle transformation; prepare the desired
scale before gating. FCS reading uses the defaults of
[`flowWorkspace::load_cytoset_from_fcs()`](https://rdrr.io/pkg/flowWorkspace/man/load_cytoset_from_fcs.html).

Use the same sample order when supplying data to downstream functions.

``` r

matrices <- list(control = control_matrix, stimulated = stimulated_matrix)
gateStim(
  pathProject = "matrix-results", .data = matrices,
  batchList = list(donor1 = c("control", "stimulated")),
  chnl = colnames(control_matrix)
)
```

### Specifying batches

`batchList` tells
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
which samples to compare. Each element is one batch (typically one donor
or sample) and holds its samples, as indices into `.data` or as sample
names. **The first sample in each element is the unstimulated sample**;
the remaining samples are the stimulated samples from the same batch. No
other argument identifies the unstimulated sample.

``` r

# Sample 3 is the unstimulated sample for donor 1 (stimulated: 1, 2),
# and sample 6 for donor 2 (stimulated: 4, 5)
batchList <- list(donor1 = c(3, 1, 2), donor2 = c(6, 4, 5))
```

Every batch needs an unstimulated sample and at least one stimulated
sample. An unstimulated sample may be shared by several batches,
provided it is first in each; a stimulated sample may belong to only one
batch.
[`getBatchList()`](https://satvilab.github.io/stimgate/reference/getBatchList.md)
builds a `batchList` in this order from a table of sample metadata.

### Basic Usage

``` r

# Basic gating workflow
gateStim(
  pathProject = "/path/to/project",
  .data = gs, # GatingSet object
  batchList = batchList,
  marker = c("IL2", "TNFa")
)

# Get statistics
stats <- getStimStats("/path/to/project")

# Get gate table
gates <- getStimGates("/path/to/project")

# Plot gates for the first batch
plots <- plotStim(
  ind = batchList[[1]],
  .data = gs,
  pathProject = "/path/to/project",
  marker = c("IL2", "TNFa")
)
```

For more detailed examples and advanced usage, please refer to the
function documentation.

[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
takes its data and marker arguments directly (`pathProject`, `.data`,
`batchList`, `marker` or `chnl`, `popGate`, `biasUns`, `bw`). All other
method, bandwidth and thresholding settings are supplied through its
`control` argument as a
[`stimControl()`](https://satvilab.github.io/stimgate/reference/stimControl.md)
object, and per-marker overrides through `markerControl`, keyed by
marker label or channel name.

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] stimgate_0.99.14
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        jsonlite_2.0.0      dplyr_1.2.1        
#>  [4] compiler_4.6.1      BiocManager_1.30.27 tidyselect_1.2.1   
#>  [7] jquerylib_0.1.4     systemfonts_1.3.2   scales_1.4.0       
#> [10] textshaping_1.0.5   yaml_2.3.12         fastmap_1.2.0      
#> [13] ggplot2_4.0.3       R6_2.6.1            generics_0.1.4     
#> [16] knitr_1.52          htmlwidgets_1.6.4   tibble_3.3.1       
#> [19] desc_1.4.3          bslib_0.12.0        pillar_1.11.1      
#> [22] RColorBrewer_1.1-3  rlang_1.3.0         cachem_1.1.0       
#> [25] xfun_0.61           fs_2.1.0            sass_0.4.10        
#> [28] S7_0.2.2            otel_0.2.0          cli_3.6.6          
#> [31] pkgdown_2.2.1       magrittr_2.0.5      digest_0.6.39      
#> [34] grid_4.6.1          lifecycle_1.0.5     vctrs_0.7.3        
#> [37] evaluate_1.0.5      glue_1.8.1          farver_2.1.2       
#> [40] ragg_1.5.2          rmarkdown_2.32      tools_4.6.1        
#> [43] pkgconfig_2.0.3     htmltools_0.5.9
```
