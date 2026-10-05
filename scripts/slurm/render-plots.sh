#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --job-name="dev-plots"
#SBATCH --partition=ada

# Render full analysis reports from saved results: simulations off, plots on.
# dev.sh submits this after an analysis's simulation jobs finish, with
# PLOT_QMD_FILES listing the QMDs to render (separated by ':', relative to the
# project root).

set -euo pipefail

project_root="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
qmd_files="${PLOT_QMD_FILES:-}"

if [[ -z "$qmd_files" ]]; then
  echo "ERROR: PLOT_QMD_FILES must name at least one QMD." >&2
  exit 1
fi

cd "$project_root"
export PROJECT_ROOT="$project_root"
export RUN_SIMULATIONS=false
export RUN_PLOTS=true

start_time=$(date +%s)
echo "HOSTNAME: $HOSTNAME"
echo "SIM_SIZE: ${SIM_SIZE:-final}"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-unknown}"
echo "PROJECT_ROOT: $project_root"
echo "QMD files: $qmd_files"

IFS=':' read -r -a qmd_file_vec <<< "$qmd_files"
for qmd_file in "${qmd_file_vec[@]}"; do
  if [[ ! -f "$qmd_file" ]]; then
    echo "ERROR: Could not find QMD: $project_root/$qmd_file" >&2
    exit 1
  fi
  echo "-------------------"
  echo "Render plots: $qmd_file"
  date
  r_expr="qmd_file <- '$qmd_file'; if (requireNamespace('quarto', quietly = TRUE)) { quarto::quarto_render(input = qmd_file) } else { status <- system2('quarto', c('render', qmd_file)); if (!identical(status, 0L)) quit(status = status) }"
  apptainer-rscript -f stimgate -- "$r_expr"
  echo "Completed rendering $qmd_file"
  date
done

end_time=$(date +%s)
duration=$((end_time - start_time))
printf "Elapsed time: %02d:%02d:%02d\n" \
  $((duration / 3600)) $(((duration % 3600) / 60)) $((duration % 60))
