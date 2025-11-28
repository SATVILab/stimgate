# Plot stimulation gate

Plot bivariate hex and univariate density plots for batches of samples,
along with their gates.

## Usage

``` r
stimgate_plot(
  ind,
  .data,
  path_project,
  marker,
  ind_lab = NULL,
  marker_lab = NULL,
  exc_min = TRUE,
  limits_expand = NULL,
  limits_equal = FALSE,
  grid = TRUE,
  grid_n_col = 2,
  show_gate = TRUE,
  min_cell = 10
)
```

## Arguments

- ind:

  numeric vector. Specifies indices in `.data` to plot.

- .data:

  GatingSet. Same GatingSet passed to `stimgate_gate`.

- path_project:

  character. Path to the project directory used for `stimgate_gate`.

- marker:

  character vector of length one or two. Specifies markers (channels,
  really) to be plotted. If only one is passed, then only univariate
  plots are created.

- ind_lab:

  named character vector. Labels for `ind` used in plot. Optional.

- marker_lab:

  named character vector. Labels for `marker` used in plot. Optional.

- exc_min:

  Logical. If `TRUE`, excludes the minimum expression values when
  processing the data. Default is `TRUE`.

- limits_expand:

  list. Expand the limits of the plot axes. Default is `NULL`.

- limits_equal:

  Logical. If TRUE, forces equal lengths of the limits.

- grid:

  Logical. If TRUE, arranges the resulting plots in a grid format using
  [`cowplot::plot_grid`](https://wilkelab.org/cowplot/reference/plot_grid.html).
  Default is `TRUE`.

- grid_n_col:

  Integer. Number of columns in grid layout.

- show_gate:

  Logical. If `TRUE`, overlays gate lines on the plots.\|\> Default is
  `TRUE`.

- min_cell:

  integer. Minimum number of cells to be plotted. Will skip plots with
  fewer cells. Default is 10.

## Value

A grid of plots if `grid` is TRUE, otherwise a list of ggplot objects.

## Examples

``` r
# Create example data and run gating
example_data <- get_example_data()
#> see ?HDCytoData and browseVignettes('HDCytoData') for documentation
#> loading from cache
#> Done
#> To reload it, use 'load_gs' function
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
#> [1] "/tmp/RtmpMKL4eh/stimgate_example/stimgate"

# Create plots
plots <- stimgate_plot(
  ind = example_data$batch_list[[1]], # indices in `gs` to plot
  .data = gs, # GatingSet
  path_project = path_project,
  marker = example_data$marker,
  grid = TRUE
)
```
