# stimgate

This is an R package intended for eventual submission to BioConductor.
We want to identify cells that have possibly responded to stimulation, by comparing the unstimulated and stimulated tubes from the same sample.

## Code standards

### Required before each commit

- Run `devtools::test()` to ensure all tests pass
- Run `styler::style_pkg()` to ensure consistent code formatting
- Run `lintr::lint_package()` to check for linting violations

### Development flow

- Build and reload: Use `devtools::load_all()` to reload the package in your R session
- Test: Run `devtools::test()` to execute all unit tests
- Coverage: Run `covr::report()` to check test coverage
- Documentation: Run `devtools::document()` to update documentation
- Full check: Run `devtools::check()` to perform a comprehensive package check

## Repository structure

- `R/`: Core R source code
  - `cp_cluster-helper.R`: Helper functions for getting thresholds lower by grouping thresholds from like distributions
  - `cp_cluster.R`: Gets thresholds lower by grouping thresholds from like distributions
  - `cp_uns_loc.R`: Gets thresholds by comparing the stim and unstim distributions to identify the point at which cells start appearing to have responded with some probability
  - `cp-sub.R`: Auxiliary functions for getting clusters
  - `cp_cut_pos_gates-helper.R`: Helper functions for getting more aggressive gates when applied to just the cytokine-positive cells.
  - `cyt_pos_gates.R`: Functions for getting more aggressive gates when applied to just the cytokine-positive cells.
  - `debug.R`: Debugging utilities
  - `ex.R`: Extract expression matrices from GatingSets
  - `fcs_write.R`: Write FCS files of just cytokine-positive cels
  - `gate_batch-helper.R`: Helper functions for gating batches of samples
  - `gate_batch.R`: Gate batches of samples
  - `gate_marker-herlper`: Helper functions for gating markers (within each marker, we gate all the batches)
  - `gate_marker.R`: Gate markers (within each marker, we gate all the batches)
  - `gate.R`: Main entry point for gating
  - `gates.R`: Extract the identified gates (thresholds)
  - `ind_batch.R`: Get the list of indices grouped by batch
  - `marker_list.R`: Get full parameters used for each marker
  - `plot_gate.R`: Plot the identified gates
  - `pos_ind.R`: Identify the indices of the cytokine-positive cells
  - `stats-helper-overall.R`: Helper functions for overall statistics
  - `stats-helper.R`: Helper functions for statistics
  - `stats.R`: Get statistics for the identified gates
- `_reference/`:
  - Reference material for the package, which at this stage is just old scripts that I didn't want to completely delete at the time.
- `.github/`: GitHub associated files
- `data-raw/`: Raw data files used for testing and examples
- `man/`: Automatically generated documentation files
- `renv/`: R package environment management files
- `tests/`: Unit tests for the package
- `_dependencies.R`: Explicitly listed dependences for `renv` to pick up via `implicit` dependencies
- `.gitignore`: Git ignore file to exclude unnecessary files from version control
- `.Rbuildignore`: R build ignore file to exclude unnecessary files from package builds
- `DESCRIPTION`: Package metadata file
- `LICENSE`: License file for the package
- `LICENSE.md`: License file in Markdown format
- `README.md`: Package overview and usage documentation
- `renv.lock`: Lock file for the `renv` package management system, capturing the state of the R package environment

## Key Guidelines

1. Begin each internal function with a `.`
2. Use `.debug` as the parameter name for debugging flags in internal functions
3. Use `.debug_msg()` for debug messages, which takes a boolean, a message and an optional value
4. Add unit tests using `testthat` for all new functionality
5. Validate inputs and provide meaningful error messages
6. Explicitly refer to all packages used, rather than using `@import` or `@importFrom`
7. Use `@export` for functions that should be available to users
