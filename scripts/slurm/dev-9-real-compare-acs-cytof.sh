#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --job-name="dev-9-acs-cytof"
#SBATCH --partition=ada

set -euo pipefail

# Under sbatch, the script may be copied into Slurm's spool directory.
# Use the submit directory, or an explicit PROJECT_ROOT passed by a launcher.
project_root="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
qmd_file="${ACS_CYTOF_QMD_FILE:-analysis/9-real-compare-acs-cytof.qmd}"
qmd_abs="$project_root/$qmd_file"

if [[ ! -f "$qmd_abs" ]]; then
  echo "ERROR: Could not find QMD: $qmd_abs" >&2
  exit 1
fi

# Analysis 9 runs populations across two R workers. Keep native-library thread
# pools single-threaded so the two-worker Slurm allocation is not oversubscribed.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

export RUN_PREPROCESSING="${RUN_PREPROCESSING:-true}"
export RUN_METHODS="${RUN_METHODS:-true}"
export RUN_PLOTS="${RUN_PLOTS:-false}"
export ACS_N_WORKERS="${ACS_N_WORKERS:-${SLURM_NTASKS:-2}}"
export PROJECT_ROOT="$project_root"

cd "$project_root"

start_time=$(date +%s)

echo "HOSTNAME: $HOSTNAME"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-unknown}"
echo "SLURM_NTASKS: ${SLURM_NTASKS:-unknown}"
echo "QMD file: $qmd_file"
echo "PROJECT_ROOT: $project_root"
echo "RUN_PREPROCESSING: $RUN_PREPROCESSING"
echo "RUN_METHODS: $RUN_METHODS"
echo "RUN_PLOTS: $RUN_PLOTS"
echo "ACS_N_WORKERS: $ACS_N_WORKERS"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"

echo " "
echo "-------------------"
echo "Render ACS CyTOF comparison"
date

r_expr="qmd_file <- '$qmd_file'; if (requireNamespace('quarto', quietly = TRUE)) { quarto::quarto_render(input = qmd_file) } else { status <- system2('quarto', c('render', qmd_file)); if (!identical(status, 0L)) quit(status = status) }"
apptainer-rscript -f stimgate -- "$r_expr"

echo "Completed ACS CyTOF comparison"
date
echo "-------------------"
echo " "

end_time=$(date +%s)
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))

echo "--- Script Duration ---"
printf "Elapsed time: %02d:%02d:%02d\n" \
  "$hours" "$minutes" "$seconds"
echo "-----------------------"
