# Write FCS files of marker-positive FCS files

Uses the gates to write FCS files of marker-positive FCS files.

## Usage

``` r
stimgate_fcs_write(
  path_project,
  .data,
  ind_batch_list,
  path_dir_save,
  chnl = NULL,
  gate_tbl = NULL,
  trans_fn = NULL,
  trans_chnl = NULL,
  combn_exc = NULL,
  gate_type_cyt_pos = "cyt",
  gate_type_single_pos = "single",
  mult = FALSE,
  gate_uns_method = "min"
)
```

## Arguments

- path_project:

  character. Path to project directory.

- .data:

  GatingSet. GatingSet object containing the flow cytometry data.

- ind_batch_list:

  list. List of indices grouped by batch.

- path_dir_save:

  character. Directory path to save the FCS files to.

- chnl:

  character vector. Specific channels to gate on.

- gate_tbl:

  data.frame. Pre-computed gate table, if available.

- trans_fn:

  function. Transformation function to apply.

- trans_chnl:

  character vector. Columns to transform.

- combn_exc:

  list. Combinations of channels to exclude.

- gate_type_cyt_pos:

  character. Gate type to use for cytokine-positive cells.

- gate_type_single_pos:

  character. Gate type to use for single-positive cells.

- mult:

  logical. Whether cells must be multi-positive.

- gate_uns_method:

  character. Method to calculate unstimulated thresholds.

## Details

This function processes flow cytometry data to identify and export
cytokine-positive cells to FCS files. It requires that gates have been
pre-computed using
[`stimgate_gate`](https://satvilab.github.io/stimgate/reference/stimgate_gate.md)
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
# batch_list <- list(batch1 = c(1, 2, 3), batch2 = c(4, 5, 6))
# where the last element in each batch is the unstimulated sample

# First, run gating to create gates
path_project <- tempfile("stimgate_project")
# stimgate_gate(
#   .data = gs,
#   path_project = path_project,
#   pop_gate = "root",
#   batch_list = batch_list,
#   marker = c("IL2", "IFNg")  # your cytokine markers
# )

# Then write FCS files of cytokine-positive cells
path_output <- tempfile("fcs_output")
# stimgate_fcs_write(
#   path_project = path_project,
#   .data = gs,
#   ind_batch_list = batch_list,
#   path_dir_save = path_output,
#   chnl = c("IL2", "IFNg")
# )

# Alternative: provide your own gate table
# gate_tbl <- data.frame(
#   chnl = c("IL2", "IFNg"),
#   marker = c("IL2", "IFNg"),
#   batch = c(1, 1),
#   ind = c(1, 1),
#   gate = c(0.5, 0.3),
#   gate_name = c("gate", "gate")
# )
# stimgate_fcs_write(
#   path_project = path_project,
#   .data = gs,
#   ind_batch_list = batch_list,
#   path_dir_save = path_output,
#   chnl = c("IL2", "IFNg"),
#   gate_tbl = gate_tbl
# )
} # }
```
