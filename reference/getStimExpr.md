# Read saved expression data from project

Read channel expression vectors saved under a project's sampleData
directory and return them as a tibble with sample metadata columns.

## Usage

``` r
getStimExpr(
  pathProject,
  .data = NULL,
  pop = NULL,
  ind = NULL,
  chnl = NULL,
  marker = NULL,
  bias = FALSE,
  excMin = FALSE,
  combnExc = NULL,
  chnlGate = NULL,
  markerGate = NULL,
  gateTypeCytPos = "cyt",
  gateTypeSinglePos = "single",
  mult = FALSE,
  gateUnsMethod = "min",
  transFn = NULL,
  transChnl = NULL,
  transMarker = NULL
)
```

## Arguments

- pathProject:

  character Path to project.

- .data:

  GatingSet or NULL GatingSet object to extract expression data from.
  Default is NULL.

- pop:

  character or NULL Population name(s). Default is detected from project
  sampleData.

- ind:

  character or NULL Index/indices of samples. Default is detected from
  project sampleData.

- chnl:

  character or NULL Channel name(s) to return. Default is detected from
  project sampleData.

- marker:

  character or NULL Marker name(s) to return. Cannot be specified with
  `chnl`. Default is NULL.

- bias:

  logical Whether to add bias to unstimulated sample used in the gating.
  Default is `FALSE`.

- excMin:

  logical Whether to exclude cells with the minimum expression for any
  channels. Default is FALSE.

- combnExc:

  list or NULL Combinations of channels to exclude. Default is NULL.

- chnlGate:

  character or NULL Channel name(s) to use for gating. Cannot be
  specified with `marker_gate`. Default is NULL.

- markerGate:

  character or NULL Marker name(s) to use for gating. Cannot be
  specified with `chnl_gate`. Default is NULL.

- gateTypeCytPos:

  character Gate type to use for cytokine-positive cells. Default is
  "cyt".

- gateTypeSinglePos:

  character Gate type to use for single-positive cells. Default is
  "single".

- mult:

  logical Whether to return only multi-functional cells (positive for
  multiple markers). Default is FALSE.

- gateUnsMethod:

  character Method for gating unstimulated cells. Default is "min".

- transFn:

  function or NULL Transformation function to apply to expression
  values. Default is NULL.

- transChnl:

  character or NULL Channel name(s) to transform when using channel
  names. Default is NULL (transforms all channels).

- transMarker:

  character or NULL Marker name(s) to transform when using marker names.
  Default is NULL (transforms all markers).

## Value

A tibble with columns `pop`, `ind` and one column per requested channel.
Rows correspond to cells.

## Examples

``` r
if (FALSE) { # \dontrun{
tmp <- tempdir()
dir.create(file.path(tmp, "sampleData", "POP1", "ind_1"),
  recursive = TRUE
)
saveRDS(
  c(1, 2, 3),
  file.path(tmp, "sampleData", "POP1", "ind_1", "chnl_BC1.rds")
)
saveRDS(
  c(4, 5, 6),
  file.path(tmp, "sampleData", "POP1", "ind_1", "chnl_BC2.rds")
)
getStimExpr(tmp)
getStimExpr(tmp, chnl = "BC1")
} # }
```
