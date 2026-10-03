# Analysis test targets

Run commands from the repository root. The package is loaded from the current
checkout, and failed expectations cause a nonzero exit status.

```sh
# Existing complete analysis integration suite
Rscript analysis/tests/run_analysis_tests.R

# List the QMD targets and the test files each runs
Rscript analysis/tests/run_qmd_tests.R --list

# One QMD, a selected set, or all ten targets
Rscript analysis/tests/run_qmd_tests.R 3
Rscript analysis/tests/run_qmd_tests.R "1,3,10"
Rscript analysis/tests/run_qmd_tests.R all
```

Targets also accept document stems or paths, such as `3-sim-bw-est-base` or
`analysis/3-sim-bw-est-base.qmd`. Each target runs its existing scientific helper
and document contract tests with bounded fixtures. These tests do not render
the full analyses, run their research simulation grids, or require the external
ACS dataset and cached research outputs. Tests requiring optional dependencies
retain their existing skip behaviour.

The **analysis-qmd-tests** GitHub Actions workflow runs these same targets on
Windows with R, Bioconductor dependencies, simcyto, Python and numpy. It has only
a `workflow_dispatch` trigger. Use **Actions → analysis-qmd-tests → Run workflow**,
choose the branch, and enter `all`, one target, or a comma/space-separated set in
the `qmds` input. For example:

```sh
gh workflow run analysis-qmd-tests.yaml --ref YOUR_BRANCH -f qmds='1,3,10'
```

GitHub makes a new manual workflow available once its definition is on the
default branch. The existing automatic analysis integration workflow remains
the complete integration-suite check.
