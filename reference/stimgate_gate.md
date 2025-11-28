# Identify cytokine-positive cells through automated gating

Main function for identifying cytokine-positive cells using
outlier-based gating to compare stimulated versus unstimulated samples.
This function implements a comprehensive workflow that identifies cells
responding to stimulation by detecting outliers in cytokine expression
distributions. The process includes density estimation, threshold
identification, clustering-based gate refinement, and generation of
comprehensive statistics and visualizations.

The function operates by comparing cytokine expression in stimulated
samples against corresponding unstimulated controls from the same
donor/batch to identify cells that have likely responded to stimulation.
It accounts for batch effects, background cytokine production, and
technical variability.

## Usage

``` r
stimgate_gate(
  path_project,
  .data,
  batch_list,
  marker,
  pop_gate = "root",
  bias_uns = NULL,
  bias_uns_factor = 1,
  exc_min = TRUE,
  cp_min = NULL,
  bw_min = NULL,
  min_cell = 100,
  max_pos_prob_x = Inf,
  gate_quant = c(0.25, 0.75),
  tol_clust = 1e-07,
  gate_combn = "min",
  marker_settings = NULL,
  calc_cyt_pos_gates = TRUE,
  calc_single_pos_gates = FALSE,
  debug = FALSE
)
```

## Arguments

- path_project:

  character. Path to project directory where all results will be saved.
  This directory will contain subdirectories for each marker with gate
  tables, statistics, and plots. The directory will be created if it
  doesn't exist.

- .data:

  GatingSet. A flowWorkspace GatingSet object containing the flow
  cytometry data with both stimulated and unstimulated samples. The
  GatingSet should have consistent channel names across all samples and
  include proper sample annotations.

- batch_list:

  list. Named list where each element contains indices of samples
  belonging to the same batch/donor. Names will be used for batch
  identification. Example: list(donor1 = 1:10, donor2 = 11:20). Proper
  batching is crucial for accurate background subtraction and gate
  identification.

- marker:

  list. List where each element specifies parameters for gating a
  specific marker. Each element should be a list containing at minimum
  the channel name (e.g., list(cut = "IL2")). Additional marker-specific
  parameters can override global defaults. The marker name should match
  channel names in the GatingSet.

- pop_gate:

  character vector. Population(s) within which to perform gating.
  Default is "root" to gate on all cells. Can specify other populations
  like "CD3+" or "CD4+" if these gates already exist in the GatingSet.

- bias_uns:

  numeric. Bias adjustment for unstimulated samples to account for
  background cytokine production. When NULL (default), no bias
  correction is applied. Positive values shift the unstimulated
  distribution higher, making gates more conservative. Typically ranges
  from 0.1 to 1.0 when used.

- bias_uns_factor:

  numeric. Multiplicative factor applied to bias_uns. Default is 1.
  Values \> 1 increase the bias effect, values \< 1 decrease it. This
  provides fine-tuning of the bias correction.

- exc_min:

  logical. Whether to exclude minimum expression values during analysis.
  Default is TRUE. Minimum values often represent technical artifacts or
  compensation spillover and should typically be excluded.

- cp_min:

  numeric. Minimum allowable cutpoint value. When NULL (default), no
  minimum is enforced. Useful for ensuring gates don't fall below known
  technical thresholds or background levels.

- bw_min:

  numeric. Minimum bandwidth for density estimation. When NULL
  (default), bandwidth is estimated automatically. Smaller values create
  more detailed density estimates but may be noisier. Typical range is
  0.01 to 0.1 on log-transformed data.

- min_cell:

  numeric. Minimum number of cells required for reliable gating. Default
  is 100. Samples with fewer cells will be skipped as they don't provide
  sufficient statistical power for accurate gate identification.

- max_pos_prob_x:

  numeric. Maximum x-value (expression level) to consider when
  calculating positive probabilities. Default is Inf (no limit). Can be
  used to exclude extremely high expression values that may represent
  doublets or artifacts.

- gate_quant:

  numeric vector. Quantiles used for gate combination when multiple
  gates are identified. Default is c(0.25, 0.75). The method specified
  in gate_combn determines how these quantiles are used (e.g., minimum
  of 25th percentiles).

- tol_clust:

  numeric. Convergence tolerance for clustering algorithms used in gate
  refinement. Default is 1e-7. Smaller values require more precise
  convergence but may increase computation time.

- gate_combn:

  character. Method for combining multiple gate candidates. Default is
  "min" to use the most conservative (lowest) gate. Other options may
  include "median" or "max" depending on the desired stringency.

- marker_settings:

  list. Optional list of additional marker-specific settings that
  override global defaults. Each element should be named with the marker
  name and contain parameter overrides. Default is NULL.

- calc_cyt_pos_gates:

  logical. Whether to calculate refined cytokine-positive gates using
  more sophisticated algorithms. Default is TRUE. When FALSE, only basic
  gates are calculated, which may be less accurate but faster.

- calc_single_pos_gates:

  logical. Whether to calculate single-positive gates for individual
  markers in addition to combination gates. Default is FALSE. Useful for
  detailed analysis of individual marker responses.

- debug:

  logical. Whether to enable detailed debug output and save intermediate
  results. Default is FALSE. When TRUE, additional files and verbose
  output are generated, useful for troubleshooting and method
  development.

## Value

character. Returns the path to the project directory where all results
have been saved. The directory structure created includes:

- `path_project/[marker_name]/`: Directory for each marker containing:

- `gate_tbl_init.rds`: Initial gate table with preliminary gates

- `gate_tbl.rds`: Final refined gate table

- `stats/`: Directory containing statistics files

- `plots/`: Directory containing visualization plots (if generated)

## Details

The function implements a multi-step workflow for identifying
cytokine-positive cells:

**Step 1: Data Preparation**

- Validates input parameters and GatingSet structure

- Completes marker specifications with default values

- Organizes samples by batch for proper background subtraction

**Step 2: Initial Gate Identification**

- Extracts expression data for each marker within specified populations

- Estimates probability densities for stimulated and unstimulated
  samples

- Identifies candidate cutpoints using outlier detection algorithms

- Applies clustering to refine gate positions across batches

**Step 3: Cytokine-Positive Gate Refinement (if calc_cyt_pos_gates =
TRUE)**

- Applies more sophisticated algorithms to refine initial gates

- Accounts for background cytokine production and technical variability

- Optimizes gates to minimize false positives while maintaining
  sensitivity

**Step 4: Single-Positive Gates (if calc_single_pos_gates = TRUE)**

- Calculates gates for individual markers independent of other markers

- Useful for understanding single-marker responses

**Step 5: Statistics Generation**

- Calculates comprehensive statistics including frequencies and
  combinations

- Generates cross-tabulations of cytokine-positive populations

- Saves results in structured format for downstream analysis

**Important Considerations:**

- Ensure stimulated and unstimulated samples are properly paired by
  batch

- Channel names in marker specifications must match GatingSet channels
  exactly

- Sufficient cell numbers (min_cell) are crucial for reliable gate
  identification

- Background bias correction (bias_uns) should be used cautiously and
  validated

- Debug mode generates extensive output useful for method validation

## See also

[`stimgate_gate_get`](https://satvilab.github.io/stimgate/reference/stimgate_gate_get.md)
for extracting gate information,
[`get_stats`](https://satvilab.github.io/stimgate/reference/get_stats.md)
for generating statistics from results,
[`stimgate_plot`](https://satvilab.github.io/stimgate/reference/stimgate_plot.md)
for visualizing identified gates,
[`stimgate_fcs_write`](https://satvilab.github.io/stimgate/reference/stimgate_fcs_write.md)
for exporting cytokine-positive cells,
[`GatingSet`](https://rdrr.io/pkg/flowWorkspace/man/GatingSet-class.html)
for GatingSet documentation

## Examples

``` r
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
#> see ?HDCytoData and browseVignettes('HDCytoData') for documentation
#> loading from cache
#> Done
#> To reload it, use 'load_gs' function
#> ----
#> getting base gates
#> ----
#> 
#> chnl: BC1(La139)Dd
#> getting pre-adjustment gates
#> batch 8 of 8
#> getting clustered and/or controlled gates
#> chnl: BC2(Pr141)Dd
#> getting pre-adjustment gates
#> batch 8 of 8
#> getting clustered and/or controlled gates
#> 
#> 
#> ----
#> getting single+ gates
#> ----
#> 
#> 
#> 
#> 
#> getting cyt combn frequencies
#> batch 8 of 8
```
