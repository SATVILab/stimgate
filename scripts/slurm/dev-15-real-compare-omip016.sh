#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --job-name="dev-15-real-compare-omip016"
#SBATCH --partition=ada

# OMIP-016 (Analysis 15): StimGate, Tailgate and F-beta against the authors' gates.
# Gating and comparator runs only: the prepared data are reused unless
# RUN_PREPROCESSING=true is set. Plots are drawn in the same render.

set -euo pipefail

project_root="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
qmd_file="analysis/15-real-compare-omip016.qmd"

if [[ ! -f "$project_root/$qmd_file" ]]; then
  echo "ERROR: Could not find QMD: $project_root/$qmd_file" >&2
  exit 1
fi

# Keep native thread pools within the single-task allocation.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export RUN_SIMULATIONS="${RUN_SIMULATIONS:-true}"
export RUN_PREPROCESSING="${RUN_PREPROCESSING:-false}"
export RUN_PLOTS="${RUN_PLOTS:-true}"
export PROJECT_ROOT="$project_root"

cd "$project_root"
start_time=$(date +%s)
echo "HOSTNAME: $HOSTNAME"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-unknown}"
echo "QMD file: $qmd_file"
echo "SIM_SIZE: not used (real-data analysis)"
echo "RUN_SIMULATIONS: $RUN_SIMULATIONS"
echo "RUN_PREPROCESSING: $RUN_PREPROCESSING"
echo "RUN_PLOTS: $RUN_PLOTS"
date

r_expr="qmd_file <- '$qmd_file'; if (requireNamespace('quarto', quietly = TRUE)) { quarto::quarto_render(input = qmd_file) } else { status <- system2('quarto', c('render', qmd_file)); if (!identical(status, 0L)) quit(status = status) }"
apptainer-rscript -f stimgate -- "$r_expr"

end_time=$(date +%s)
duration=$((end_time - start_time))
printf "Elapsed time: %02d:%02d:%02d\n" \
  $((duration / 3600)) $(((duration % 3600) / 60)) $((duration % 60))
