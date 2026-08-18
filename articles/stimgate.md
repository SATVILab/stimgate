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

- `stimgate_gate()`: Main function to identify cytokine-positive cells
  by gating
- `get_stats()`: Generate statistics from gating results
- `stimgate_plot()`: Visualize identified gates
- `stimgate_gate_get()`: Extract gate information
- `stimgate_fcs_write()`: Write FCS files of cytokine-positive cells

### Basic Usage

``` r

# Basic gating workflow
result <- stimgate_gate(
  pathProject = "/path/to/project",
  .data = gs, # GatingSet object
  batchList = list(batch1 = 1:10, batch2 = 11:20),
  marker = list(
    list(cut = "IL2", tol = 0.5e-8),
    list(cut = "TNFa", tol = 0.5e-8)
  )
)

# Get statistics
stats <- get_stats("/path/to/project")

# Get gate table
gates <- get_gate_tbl("/path/to/project")

# Plot gates
plots <- stimgate_plot(
  ind = 1:3,
  .data = gs,
  pathProject = "/path/to/project",
  marker = c("IL2", "TNFa")
)
```

For more detailed examples and advanced usage, please refer to the
function documentation.

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
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
#> [1] stimgate_0.104.0-2
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        jsonlite_2.0.0      dplyr_1.2.1        
#>  [4] compiler_4.6.1      BiocManager_1.30.27 tidyselect_1.2.1   
#>  [7] jquerylib_0.1.4     systemfonts_1.3.2   scales_1.4.0       
#> [10] textshaping_1.0.5   yaml_2.3.12         fastmap_1.2.0      
#> [13] ggplot2_4.0.3       R6_2.6.1            generics_0.1.4     
#> [16] knitr_1.51          htmlwidgets_1.6.4   tibble_3.3.1       
#> [19] desc_1.4.3          bslib_0.12.0        pillar_1.11.1      
#> [22] RColorBrewer_1.1-3  rlang_1.3.0         cachem_1.1.0       
#> [25] xfun_0.60           fs_2.1.0            sass_0.4.10        
#> [28] S7_0.2.2            otel_0.2.0          cli_3.6.6          
#> [31] pkgdown_2.2.1       magrittr_2.0.5      digest_0.6.39      
#> [34] grid_4.6.1          lifecycle_1.0.5     vctrs_0.7.3        
#> [37] evaluate_1.0.5      glue_1.8.1          farver_2.1.2       
#> [40] ragg_1.5.2          rmarkdown_2.31      tools_4.6.1        
#> [43] pkgconfig_2.0.3     htmltools_0.5.9
```
