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
export RUN_PREPROCESSING=false
export RUN_STIMGATE=false
export RUN_COMPARATORS=false
export RUN_METHODS=false
# Match the simulation launchers: single-threaded BLAS keeps reruns (e.g. the
# OMIP-111 diagnostic pages) bit-identical to the saved results.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
# Empty for manual renders; submitted reports must use this submission's run.
export ANALYSIS_EXPECTED_RUN_ID="${ANALYSIS_RUN_ID:-}"

start_time=$(date +%s)
echo "HOSTNAME: $HOSTNAME"
echo "SIM_SIZE: ${SIM_SIZE:-parameters.sim_size in _projr.yml}"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-unknown}"
echo "PROJECT_ROOT: $project_root"
echo "QMD files: $qmd_files"

IFS=':' read -r -a qmd_file_vec <<< "$qmd_files"
for qmd_file in "${qmd_file_vec[@]}"; do
  if [[ ! -f "$qmd_file" ]]; then
    echo "ERROR: Could not find QMD: $project_root/$qmd_file" >&2
    exit 1
  fi
  # Render directly, without projr's build/clean step or the isolated chunk
  # helper (which deletes its HTML on exit). All reports embed resources;
  # each HTML is self-contained even if Quarto reuses input-stem _files.
  qmd_stem=$(basename -- "$qmd_file" .qmd)
  # QMDs with Monte Carlo intervals render with them on only; set SHOW_MCSE=off
  # in a manual render for the version without them. QMDs without intervals
  # render under their usual name.
  if grep -q "show_mcse" "$qmd_file"; then mcse_modes=(on); else mcse_modes=(none); fi
  for mcse_mode in "${mcse_modes[@]}"; do
    if [[ "$mcse_mode" == "none" ]]; then
      unset SHOW_MCSE
      output_file="${qmd_stem}.html"
    else
      export SHOW_MCSE="$mcse_mode"
      output_file="${qmd_stem}-mcse_${mcse_mode}.html"
    fi
    echo "Render plots: $qmd_file ($mcse_mode) -> $output_file"
    date
    # Render under Quarto's default name and rename afterwards: --output with
    # embed-resources makes Quarto look for its support files in the wrong place.
    r_expr="qmd_file <- '$qmd_file'; if (requireNamespace('quarto', quietly = TRUE)) { quarto::quarto_render(input = qmd_file) } else { status <- system2('quarto', c('render', qmd_file)); if (!identical(status, 0L)) quit(status = status) }"
    apptainer-rscript -f stimgate -- "$r_expr"
    default_html="$(dirname -- "$qmd_file")/${qmd_stem}.html"
    if [[ "$output_file" != "${qmd_stem}.html" ]]; then
      mv -f -- "$default_html" "$(dirname -- "$qmd_file")/$output_file"
    fi
    echo "Completed rendering $output_file"
    date
  done
done

end_time=$(date +%s)
duration=$((end_time - start_time))
printf "Elapsed time: %02d:%02d:%02d\n" \
  $((duration / 3600)) $(((duration % 3600) / 60)) $((duration % 60))
