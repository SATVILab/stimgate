# Write FCS files of marker-positive FCS files

Uses the gates to write FCS files of marker-positive FCS files.

## Usage

``` r
writeStimFCS(
  pathProject,
  .data,
  pop = NULL,
  indBatchList,
  pathDirSave,
  chnl = NULL,
  gateTbl = NULL,
  transFn = NULL,
  transChnl = NULL,
  combnExc = NULL,
  gateTypeCytPos = "cyt",
  gateTypeSinglePos = "single",
  mult = FALSE,
  gateUnsMethod = "min"
)
```

## Arguments

- pathProject:

  character. Path to project directory.

- .data:

  GatingSet. GatingSet object containing the flow cytometry data.

- pop:

  character. Population that was gated on.

- indBatchList:

  list. List of indices grouped by batch.

- pathDirSave:

  character. Directory path to save the FCS files to.

- chnl:

  character vector. Specific channels to gate on.

- gateTbl:

  data.frame. Pre-computed gate table, if available.

- transFn:

  function. Transformation function to apply.

- transChnl:

  character vector. Columns to transform.

- combnExc:

  list. Combinations of channels to exclude.

- gateTypeCytPos:

  character. Gate type to use for cytokine-positive cells.

- gateTypeSinglePos:

  character. Gate type to use for single-positive cells.

- mult:

  logical. Whether cells must be multi-positive.

- gateUnsMethod:

  character. Method to calculate unstimulated thresholds.

## Details

This function processes flow cytometry data to identify and export
cytokine-positive cells to FCS files. It requires that gates have been
pre-computed using
[`gateStim`](https://satvilab.github.io/stimgate/reference/gateStim.md)
or that a complete gate table is provided.

The function will create the output directory and write FCS files for
samples that contain cytokine-positive cells. If no positive cells are
found in a sample, no FCS file will be written for that sample.

## Examples

``` r
if (FALSE) { # \dontrun{
# Complete workflow example
library(stimgate)

# Load your GatingSet (gs) and define batch structure
# batchList <- list(batch1 = c(1, 2, 3), batch2 = c(4, 5, 6))
# where the first element in each batch is the unstimulated sample

# First, run gating to create gates
pathProject <- tempfile("stimgate_project")
# gateStim(
#   .data = gs,
#   pathProject = pathProject,
#   popGate = "root",
#   batchList = batchList,
#   marker = c("IL2", "IFNg")  # your cytokine markers
# )

# Then write FCS files of cytokine-positive cells
pathOutput <- tempfile("fcs_output")
# writeStimFCS(
#   pathProject = pathProject,
#   .data = gs,
#   indBatchList = batchList,
#   pathDirSave = pathOutput,
#   chnl = c("IL2", "IFNg")
# )

# Alternative: provide your own gate table
# gateTbl <- data.frame(
#   chnl = c("IL2", "IFNg"),
#   marker = c("IL2", "IFNg"),
#   batch = c(1, 1),
#   ind = c(1, 1),
#   gate = c(0.5, 0.3),
#   gateName = c("gate", "gate")
# )
# writeStimFCS(
#   pathProject = pathProject,
#   .data = gs,
#   indBatchList = batchList,
#   pathDirSave = pathOutput,
#   chnl = c("IL2", "IFNg"),
#   gateTbl = gateTbl
# )
} # }
```
