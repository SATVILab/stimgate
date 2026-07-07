#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --job-name="dev-2-bw-freq_bs"
#SBATCH --partition=ada

set -euo pipefail

chunk_index="${1:-${SIM_GRID_CHUNK_INDEX:-1}}"
n_chunks="${SIM_GRID_N_CHUNKS:-4}"

# Under sbatch, the script may be copied into Slurm's spool directory.
# So do not infer the project root from ${BASH_SOURCE[0]}. Use the submit
# directory, or an explicit PROJECT_ROOT passed from the launcher.
project_root="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
qmd_file="${SIM_GRID_QMD_FILE:-analysis/split/2-sim-bw-freq_bs-global-${chunk_index}.qmd}"
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
echo "Run projr"
date
apptainer-rscript -f stimgate -- "projr::projr_build_dev(profile = '', file = '$qmd_file')"
echo "Completed running projr"
date
echo "-------------------"
echo " "

end_time=$(date +%s)
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

echo "--- Script Duration ---"
printf "Elapsed time: %02d:%02d:%02d\n" $hours $minutes $seconds
echo "-----------------------"
