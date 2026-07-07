#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --job-name="dev-6-freq_bs"
#SBATCH --partition=ada

set -euo pipefail

chunk_index="${1:-${SIM_GRID_CHUNK_INDEX:-1}}"
n_chunks="${SIM_GRID_N_CHUNKS:-4}"

# Under sbatch, the script may be copied into Slurm's spool directory.
# So do not infer the project root from ${BASH_SOURCE[0]}. Use the submit
# directory, or an explicit PROJECT_ROOT passed from the launcher.
project_root="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
qmd_file="${SIM_GRID_QMD_FILE:-analysis/split/6-sim-bw-freq_bs-adaptive-${chunk_index}.qmd}"
qmd_abs="$project_root/$qmd_file"

if [[ ! "$n_chunks" =~ ^[0-9]+$ ]] || (( n_chunks < 1 )); then
  echo "ERROR: SIM_GRID_N_CHUNKS must be a positive integer. Got: $n_chunks" >&2
  exit 1
fi

if [[ ! "$chunk_index" =~ ^[0-9]+$ ]] || (( chunk_index < 1 || chunk_index > n_chunks )); then
  echo "ERROR: chunk_index must be an integer between 1 and $n_chunks. Got: $chunk_index" >&2
  exit 1
fi

if [[ ! -f "$qmd_abs" ]]; then
  echo "ERROR: Could not find split QMD: $qmd_abs" >&2
  exit 1
fi

export SIM_GRID_CHUNK_INDEX="$chunk_index"
export SIM_GRID_N_CHUNKS="$n_chunks"
export RUN_SIMULATIONS="${RUN_SIMULATIONS:-true}"
export RUN_PLOTS="${RUN_PLOTS:-false}"
export PROJECT_ROOT="$project_root"

cd "$project_root"

# Record the start time
start_time=$(date +%s)

echo "HOSTNAME: $HOSTNAME"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-unknown}"
echo "SIM_GRID_CHUNK_INDEX: $SIM_GRID_CHUNK_INDEX"
echo "SIM_GRID_N_CHUNKS: $SIM_GRID_N_CHUNKS"
echo "QMD file: $qmd_file"
echo "PROJECT_ROOT: $project_root"

echo " "
echo " "
echo " "

echo "-------------------"
echo "Render Quarto directly"
date
r_expr="qmd_file <- '$qmd_file'; if (requireNamespace('quarto', quietly = TRUE)) { quarto::quarto_render(input = qmd_file) } else { status <- system2('quarto', c('render', qmd_file)); if (!identical(status, 0L)) quit(status = status) }"
apptainer-rscript -f stimgate -- "$r_expr"
echo "Completed rendering Quarto"
date
echo "-------------------"
echo " "

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Convert duration to human-readable format
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

# Append the duration to the Slurm standard output log
echo "--- Script Duration ---"
printf "Elapsed time: %02d:%02d:%02d\n" $hours $minutes $seconds
echo "-----------------------"
