# BUILDLOG

### v0.102.0: Miguel J Rodo (2026-06-29 11:59:11)

**Description**

Run sim transformation and bandwidth estimation scripts

**Metadata**

- Total time: 
37s
- `projr` profile: 
default

**System Resources**

- OS: Linux 5.14.0-611.55.1.el9_7.x86_64
- OS Version: #1 SMP PREEMPT_DYNAMIC Tue May 12 18:04:19 UTC 2026
- Architecture: x86_64
- Platform: x86_64-pc-linux-gnu
- CPU Cores: 48
- Total RAM: 375Gi
- Disk Space: 90T total, 21T available

**`projr` config**

```yaml
directories:
  raw-data:
    path: analysis/_raw_data
  output:
    path: analysis/_output
  docs:
    path: analysis/docs
build:
  scripts:
  - analysis/1-sim-trans.qmd
  - analysis/2-sim-bw-freq_bs.qmd

```

**Session info**

```
R version 4.4.2 (2024-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
[1] stimgate_0.102.0

loaded via a namespace (and not attached):
 [1] gtable_0.3.6         xfun_0.59            ggplot2_4.0.3       
 [4] ks_1.15.2            gh_1.6.0             processx_3.9.0      
 [7] Biobase_2.66.0       lattice_0.22-9       projr_0.5.2-10      
[10] vctrs_0.7.3          tools_4.4.2          ps_1.9.3            
[13] generics_0.1.4       parallel_4.4.2       stats4_4.4.2        
[16] flowWorkspace_4.18.1 tibble_3.3.1         cluster_2.1.8.2     
[19] pkgconfig_2.0.3      KernSmooth_2.23-26   Matrix_1.7-5        
[22] data.table_1.18.4    RColorBrewer_1.1-3   S7_0.2.2            
[25] desc_1.4.3           S4Vectors_0.44.0     graph_1.84.1        
[28] lifecycle_1.0.5      scam_1.2-22          compiler_4.4.2      
[31] farver_2.1.2         stringr_1.6.0        htmltools_0.5.9     
[34] yaml_2.3.12          pracma_2.4.6         flowCore_2.18.0     
[37] later_1.4.8          pillar_1.11.1        tidyr_1.3.2         
[40] mclust_6.1.2         nlme_3.1-169         RProtoBufLib_2.18.0 
[43] commonmark_2.0.0     gtools_3.9.5         tidyselect_1.2.1    
[46] digest_0.6.39        mvtnorm_1.4-1        stringi_1.8.7       
[49] dplyr_1.2.1          purrr_1.2.2          quarto_1.5.1        
[52] splines_4.4.2        cowplot_1.2.0        rprojroot_2.1.1     
[55] fastmap_1.2.0        grid_4.4.2           cli_3.6.6           
[58] magrittr_2.0.5       ncdfFlow_2.52.1      XML_3.99-0.23       
[61] pkgbuild_1.4.8       withr_3.0.3          scales_1.4.0        
[64] roxygen2_8.0.0       rmarkdown_2.31       httr_1.4.8          
[67] matrixStats_1.5.0    otel_0.2.0           cytolib_2.18.2      
[70] evaluate_1.0.5       knitr_1.51           mgcv_1.9-4          
[73] rlang_1.2.0          Rcpp_1.1.1-1.1       glue_1.8.1          
[76] Rgraphviz_2.50.0     BiocManager_1.30.27  xml2_1.6.0          
[79] renv_1.1.2           BiocGenerics_0.52.0  pkgload_1.5.3       
[82] rstudioapi_0.19.0    jsonlite_2.0.0       R6_2.6.1            
[85] fs_2.1.0            
```

----

