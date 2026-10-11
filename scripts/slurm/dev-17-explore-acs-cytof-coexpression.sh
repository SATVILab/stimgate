#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=12:00:00
#SBATCH --job-name="dev-17-acs-coexpression"
#SBATCH --partition=ada

# Draw Analysis 17's co-expression views from Analysis 9's cached gates and
# GatingSets. It gates nothing, so it renders with plots on and has no plot
# job; dev.sh starts it only after Analysis 9's job succeeds when both are
# submitted together.

set -euo pipefail

project_root="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
qmd_file="analysis/17-explore-acs-cytof-coexpression.qmd"

if [[ ! -f "$project_root/$qmd_file" ]]; then
  echo "ERROR: Could not find QMD: $project_root/$qmd_file" >&2
  exit 1
fi

# Keep native thread pools within the allocation.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export RUN_SIMULATIONS="${RUN_SIMULATIONS:-false}"
export RUN_PLOTS="${RUN_PLOTS:-true}"
export PROJECT_ROOT="$project_root"

cd "$project_root"
start_time=$(date +%s)
echo "HOSTNAME: $HOSTNAME"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-unknown}"
echo "QMD file: $qmd_file"
echo "RUN_PLOTS: $RUN_PLOTS"
date

r_expr="qmd_file <- '$qmd_file'; if (requireNamespace('quarto', quietly = TRUE)) { quarto::quarto_render(input = qmd_file) } else { status <- system2('quarto', c('render', qmd_file)); if (!identical(status, 0L)) quit(status = status) }"
apptainer-rscript -f stimgate -- "$r_expr"

end_time=$(date +%s)
duration=$((end_time - start_time))
printf "Elapsed time: %02d:%02d:%02d\n" \
  $((duration / 3600)) $(((duration % 3600) / 60)) $((duration % 60))
