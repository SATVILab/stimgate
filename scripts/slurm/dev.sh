#!/usr/bin/env bash
set -euo pipefail

# get location of script
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
project_root=$(cd -- "$script_dir/../.." &> /dev/null && pwd)

scripts=(
  # "dev-1-sim-trans.sh"
  # "dev-2-stim-bw-freq_bs-global.sh"
  # "dev-3-sim-bw-est-base.sh"
  # "dev-4-sim-bw-est-norm.sh"
  # "dev-5-sim-bw-est-adaptive.sh"
  # "dev-6-sim-bw-freq_bs-adaptive.sh"
  # "dev-7-sim-compare-freq_bs.sh"
  # "dev-8-sim-compare-freq_bs-batch.sh"
  "dev-9-real-compare-acs-cytof.sh"
)

poll_seconds="${POLL_SECONDS:-5}"
sim_grid_n_chunks="${SIM_GRID_N_CHUNKS:-4}"
sim_grid_shuffle_seed="${SIM_GRID_SHUFFLE_SEED:-20260707}"
analysis_run_id="${ANALYSIS_RUN_ID:-analysis-slurm-$(date -u +%Y%m%dT%H%M%S)-$$}"

install_script="$script_dir/install.sh"

chunked_qmd_stem_for_script() {
  case "$1" in
    dev-2-stim-bw-freq_bs-global.sh)
      echo "2-sim-bw-freq_bs-global"
      ;;
    dev-3-sim-bw-est-base.sh)
      echo "3-sim-bw-est-base"
      ;;
    dev-4-sim-bw-est-norm.sh)
      echo "4-sim-bw-est-norm"
      ;;
    dev-5-sim-bw-est-adaptive.sh)
      echo "5-sim-bw-est-adaptive"
      ;;
    dev-6-sim-bw-freq_bs-adaptive.sh)
      echo "6-sim-bw-freq_bs-adaptive"
      ;;
    *)
      echo ""
      ;;
  esac
}

tmp_before="$(mktemp)"
trap 'rm -f "$tmp_before"' EXIT

echo "Recording currently active Slurm jobs"
squeue -h -u "$USER" -o "%A" | sort -u > "$tmp_before"

echo "Submitting install job"
slurm-sbatch "$install_script"

echo "Finding newly submitted install job"

install_jobid=""

for attempt in {1..30}; do
  install_jobid="$(
    awk '
      FNR == NR {
        before[$1] = 1
        next
      }

      $2 == "install" && !($1 in before) {
        print $1
      }
    ' "$tmp_before" <(squeue -h -u "$USER" -o "%A %j") | tail -n 1
  )"

  if [[ -n "$install_jobid" ]]; then
    break
  fi

  sleep 2
done

if [[ -z "$install_jobid" ]]; then
  echo "ERROR: Could not identify the new install job."
  echo "Current jobs:"
  squeue -u "$USER"
  exit 1
fi

echo "Install job ID: $install_jobid"
echo "Waiting for install job to finish"

while squeue -h -j "$install_jobid" 2>/dev/null | grep -q .; do
  date
  squeue -j "$install_jobid"
  sleep "$poll_seconds"
done

echo "Install job has left the queue"
echo "Checking final state"

install_state=""

if command -v sacct >/dev/null 2>&1; then
  for attempt in {1..20}; do
    install_state="$(
      sacct -j "$install_jobid" -n -X -o State 2>/dev/null |
        awk 'NF { print $1; exit }'
    )"

    if [[ -n "$install_state" ]]; then
      break
    fi

    sleep 3
  done
else
  echo "WARNING: sacct not available, so final install status cannot be checked."
fi

if [[ -n "$install_state" && "$install_state" != "COMPLETED" ]]; then
  echo "ERROR: install job did not complete successfully."
  echo "Final state: $install_state"
  echo
  echo "sacct output:"
  sacct -j "$install_jobid" -o JobID,JobName,State,ExitCode,Elapsed
  exit 1
fi

if [[ "$install_state" == "COMPLETED" ]]; then
  echo "Install job completed successfully"
fi

echo "Submitting downstream jobs"

for script in "${scripts[@]}"; do
  qmd_stem="$(chunked_qmd_stem_for_script "$script")"

  if [[ -n "$qmd_stem" ]]; then
    echo "Submitting transactional chunks for $script"
    echo "ANALYSIS_RUN_ID: $analysis_run_id"

    for chunk_index in $(seq 1 "$sim_grid_n_chunks"); do
      log_dir="_tmp/log/sbatch/${qmd_stem}/chunk-${chunk_index}"
      job_name="${qmd_stem}-${chunk_index}"

      echo "Submitting $script chunk $chunk_index of $sim_grid_n_chunks"
      echo "Log directory: $log_dir"

      slurm-sbatch -l "$log_dir" -n "$script_dir/$script" -- \
        --job-name="$job_name" \
        --export=ALL,PROJECT_ROOT="$project_root",ANALYSIS_RUN_ID="$analysis_run_id",SIM_GRID_CHUNK_INDEX="$chunk_index",SIM_GRID_N_CHUNKS="$sim_grid_n_chunks",SIM_GRID_SHUFFLE_SEED="$sim_grid_shuffle_seed",RUN_SIMULATIONS=true,RUN_PLOTS=false
    done
  else
    echo "Submitting $script"
    slurm-sbatch "$script_dir/$script"
  fi
done

echo "All downstream jobs submitted"
