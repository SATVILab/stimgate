# Analysis test targets

With `mode: render`, the manual-only **analysis-qmd-tests** workflow renders QMDs end to end instead (quick mode by default, optionally dev): one parallel Windows job per QMD runs the simulation render, then the plot render, and uploads the HTML and saved figures. `qmds: all` renders every simulation QMD except 5 and 6. For example: `gh workflow run analysis-qmd-tests.yaml -f mode=render -f qmds='2a,7'`.

Run commands from the repository root. The package is loaded from the current
checkout, and failed expectations cause a nonzero exit status.

```sh
# Existing complete analysis integration suite
Rscript analysis/tests/run_analysis_tests.R

# List the QMD targets and the test files each runs
Rscript analysis/tests/run_qmd_tests.R --list

# One QMD, a selected set, or all eleven targets
Rscript analysis/tests/run_qmd_tests.R 3
Rscript analysis/tests/run_qmd_tests.R "1,3,10"
Rscript analysis/tests/run_qmd_tests.R all
```

Analysis 2 is split into targets `2a` and `2b`. For example,
`Rscript analysis/tests/run_qmd_tests.R "2a,2b"` runs both.
Both targets include the shared seeded runner, resume and promotion tests.
Run state lives in `cache/sim/<analysis-key>/runs/<YYYY-MM-DD>/<run-id>/`,
beside `staging/` and `current/`. Runtime tests check this layout, exclusion of
run state from promotion, and resumption at legacy `cache/log/analysis/...`
paths recorded in staged manifests' `path_log_run` field.
Their plot fixtures check relative-error summaries averaged over cell counts
and separate outputs for each cell count without rendering the research grid.

Submit the corresponding Slurm analyses independently or together:

```sh
bash scripts/slurm/dev.sh 2a
bash scripts/slurm/dev.sh 2b
bash scripts/slurm/dev.sh 2a 2b
SIM_SIZE=final bash scripts/slurm/dev.sh 2a  # full replicate counts for reported results
```

Full-grid simulation runs default to draft: fewer replicates on the same grid,
with results and figures kept in draft folders. Analyses 7/8 keep 20 samples
in each dataset and use five datasets (final uses twenty). This is separate
from the tiny quick profile used for smoke checks.

The Slurm launcher accepts analysis IDs or launcher filenames. With no arguments,
it submits its default batch. All requested targets are validated before any job
is submitted. Each selected chunk receives the same logical run ID.
`bash analysis/tests/test-slurm-launchers.sh` checks submission and render
arguments using mock commands, without Slurm or R; both analysis CI workflows
run this check.

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

The bounded performance fixture renders five constructed-data plots (pooled
signed percentiles and ratio, dataset maximum occurrence, conditional signed
severity and ratio), plus a coverage table, without generating or gating data. It writes both
versions under `mcse_off/` and `mcse_on/`:

```sh
Rscript --no-init-file analysis/tests/render_performance_fixture.R /tmp/stimgate-performance-figures
```

The numerical regressions distinguish pooled tube percentiles from means of
within-dataset percentiles, recompute whole-dataset bootstrap statistics with
multiplicity and paired biological seeds, and keep incomplete datasets separate
from unaffected datasets. These fixtures do not validate full twenty-sample
research simulations; those require fresh HPC runs.
