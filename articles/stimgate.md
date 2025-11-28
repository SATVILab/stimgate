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

- [`stimgate_gate()`](https://satvilab.github.io/stimgate/reference/stimgate_gate.md):
  Main function to identify cytokine-positive cells by gating
- [`get_stats()`](https://satvilab.github.io/stimgate/reference/get_stats.md):
  Generate statistics from gating results
- [`stimgate_plot()`](https://satvilab.github.io/stimgate/reference/stimgate_plot.md):
  Visualize identified gates
- [`stimgate_gate_get()`](https://satvilab.github.io/stimgate/reference/stimgate_gate_get.md):
  Extract gate information
- [`stimgate_fcs_write()`](https://satvilab.github.io/stimgate/reference/stimgate_fcs_write.md):
  Write FCS files of cytokine-positive cells

### Basic Usage

``` r
# Basic gating workflow
result <- stimgate_gate(
  path_project = "/path/to/project",
  .data = gs, # GatingSet object
  batch_list = list(batch1 = 1:10, batch2 = 11:20),
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
  path_project = "/path/to/project",
  marker = c("IL2", "TNFa")
)
```

For more detailed examples and advanced usage, please refer to the
function documentation.

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
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
#> [1] stimgate_0.99.1-1
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        jsonlite_2.0.0      dplyr_1.1.4        
#>  [4] compiler_4.5.2      BiocManager_1.30.27 tidyselect_1.2.1   
#>  [7] Biobase_2.70.0      cytolib_2.22.0      cluster_2.1.8.1    
#> [10] jquerylib_0.1.4     systemfonts_1.3.1   scales_1.4.0       
#> [13] textshaping_1.0.4   yaml_2.3.10         fastmap_1.2.0      
#> [16] RProtoBufLib_2.22.0 ggplot2_4.0.1       R6_2.6.1           
#> [19] generics_0.1.4      knitr_1.50          BiocGenerics_0.56.0
#> [22] tibble_3.3.0        desc_1.4.3          bslib_0.9.0        
#> [25] pillar_1.11.1       RColorBrewer_1.1-3  rlang_1.1.6        
#> [28] cachem_1.1.0        flowCore_2.22.0     xfun_0.54          
#> [31] fs_1.6.6            sass_0.4.10         S7_0.2.1           
#> [34] cli_3.6.5           pkgdown_2.2.0       magrittr_2.0.4     
#> [37] digest_0.6.39       grid_4.5.2          lifecycle_1.0.4    
#> [40] S4Vectors_0.48.0    vctrs_0.6.5         evaluate_1.0.5     
#> [43] glue_1.8.0          farver_2.1.2        ragg_1.5.0         
#> [46] stats4_4.5.2        rmarkdown_2.30      matrixStats_1.5.0  
#> [49] tools_4.5.2         pkgconfig_2.0.3     htmltools_0.5.8.1
```
