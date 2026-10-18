#!/usr/bin/env bash
set -euo pipefail

install_stimgate=FALSE

# get location of script
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
project_root=$(cd -- "$script_dir/../.." &> /dev/null && pwd)

scripts=(
  "dev-1-sim-trans.sh"
  "dev-2a-stim-bw-freq_bs-global.sh"
  "dev-2b-stim-bias_uns-freq_bs.sh"
  "dev-3-sim-bw-est-base.sh"
  "dev-4-sim-bw-est-norm.sh"
  "dev-6-sim-tune-comparators.sh"
  "dev-7-sim-compare-freq_bs.sh"
  "dev-7b-sim-compare-null.sh"
  "dev-8-sim-compare-freq_bs-batch.sh"
  "dev-9a-real-tune-acs-cytof-tailgate.sh"
  "dev-9b-real-compare-acs-cytof.sh"
  "dev-11a-sim-low-separation-cyt-pos.sh"
  "dev-11b-sim-low-separation-batch.sh"
  "dev-12-sim-cluster-gates.sh"
  "dev-13-real-debug-acs-cytof.sh"
  "dev-14-real-compare-omip111.sh"
  "dev-15-real-compare-omip016.sh"
  "dev-17-explore-acs-cytof-coexpression.sh"
)

# With arguments, submit only those analysis IDs or launcher filenames.
# Validate the entire selection before submitting any jobs.
if [[ "${1:-}" == "--help" ]]; then
  echo "Usage: bash scripts/slurm/dev.sh [analysis ID or launcher filename ...]"
  echo "Examples: dev.sh 2a; dev.sh 2b; dev.sh 2a 2b; dev.sh 11a 11b"
  echo "Simulation size: parameters.sim_size in _projr.yml (draft or final)."
  echo "SIM_SIZE overrides it for this submission. Draft results are kept separately."
  echo "Example: SIM_SIZE=final bash scripts/slurm/dev.sh 2a"
  exit 0
fi
# Without an explicit override, each render reads the active projr parameter.
sim_size="${SIM_SIZE:-}"
if [[ -n "$sim_size" && "$sim_size" != "final" && "$sim_size" != "draft" ]]; then
  echo "ERROR: SIM_SIZE must be final or draft. Got: $sim_size" >&2
  exit 1
fi
sim_size_export=""
if [[ -n "$sim_size" ]]; then
  export SIM_SIZE="$sim_size"
  sim_size_export=",SIM_SIZE=$sim_size"
fi
if (( $# > 0 )); then
  scripts=()
  for target in "$@"; do
    matched_script=""
    for launcher_path in "$script_dir"/dev-*.sh; do
      launcher_name="${launcher_path##*/}"
      analysis_id="${launcher_name#dev-}"
      analysis_id="${analysis_id%%-*}"
      if [[ "$target" == "$analysis_id" || "$target" == "$launcher_name" ]]; then
        matched_script="$launcher_name"
        break
      fi
    done
    if [[ -z "$matched_script" ]]; then
      echo "ERROR: Unknown analysis target: $target. Use an ID such as 2a or 2b." >&2
      exit 1
    fi
    if [[ " ${scripts[*]} " != *" $matched_script "* ]]; then
      scripts+=("$matched_script")
    fi
  done
fi

# Some analyses read another's results, so submit the one they read first
# when both run: 9b reads 9a, and 13, 17 and 18 read 9b (and 9a's
# selection). The OMIP
# analyses and 17 also wait for Analysis 1's projr build, which clears the
# figure folder.
move_after() {
  local first="$1" second="$2" idx_first=-1 idx_second=-1 i moved
  for i in "${!scripts[@]}"; do
    [[ "${scripts[$i]}" == "$first" ]] && idx_first=$i
    [[ "${scripts[$i]}" == "$second" ]] && idx_second=$i
  done
  if (( idx_first >= 0 && idx_second >= 0 && idx_second < idx_first )); then
    moved="${scripts[$idx_second]}"
    unset 'scripts[idx_second]'
    scripts=("${scripts[@]}" "$moved")
  fi
}
for omip_script in dev-14-real-compare-omip111.sh dev-15-real-compare-omip016.sh \
  dev-17-explore-acs-cytof-coexpression.sh; do
  move_after dev-1-sim-trans.sh "$omip_script"
done
move_after dev-9a-real-tune-acs-cytof-tailgate.sh dev-9b-real-compare-acs-cytof.sh
move_after dev-9b-real-compare-acs-cytof.sh dev-13-real-debug-acs-cytof.sh
move_after dev-9b-real-compare-acs-cytof.sh dev-17-explore-acs-cytof-coexpression.sh
move_after dev-9b-real-compare-acs-cytof.sh dev-18-explore-acs-cytof-coexpression-gates.sh

poll_seconds="${POLL_SECONDS:-5}"
sim_grid_n_chunks_set="${SIM_GRID_N_CHUNKS:-}"
sim_grid_n_chunks="${SIM_GRID_N_CHUNKS:-4}"
if [[ ! "$sim_grid_n_chunks" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: SIM_GRID_N_CHUNKS must be a positive integer. Got: $sim_grid_n_chunks" >&2
  exit 1
fi

# Chunks per chunked analysis: SIM_GRID_N_CHUNKS when the caller sets it,
# otherwise 4, except Analysis 8, whose crossed grid of 42 baselines needs 40
# chunks of about an hour each.
n_chunks_for_script() {
  if [[ -n "$sim_grid_n_chunks_set" ]]; then
    echo "$sim_grid_n_chunks"
    return
  fi
  case "$1" in
    dev-8-sim-compare-freq_bs-batch.sh)
      echo 40
      ;;
    *)
      echo "$sim_grid_n_chunks"
      ;;
  esac
}
sim_grid_shuffle_seed="${SIM_GRID_SHUFFLE_SEED:-20260707}"
analysis_run_id="${ANALYSIS_RUN_ID:-analysis-slurm-$(date -u +%Y%m%dT%H%M%S)-$$}"

install_script="$script_dir/install.sh"

chunked_qmd_stem_for_script() {
  case "$1" in
    dev-2a-stim-bw-freq_bs-global.sh)
      echo "2a-sim-bw-freq_bs-global"
      ;;
    dev-2b-stim-bias_uns-freq_bs.sh)
      echo "2b-sim-bias_uns-freq_bs"
      ;;
    dev-3-sim-bw-est-base.sh)
      echo "3-sim-bw-est-base"
      ;;
    dev-4-sim-bw-est-norm.sh)
      echo "4-sim-bw-est-norm"
      ;;
    dev-7-sim-compare-freq_bs.sh)
      echo "7-sim-compare-freq_bs"
      ;;
    dev-7b-sim-compare-null.sh)
      echo "7b-sim-compare-null"
      ;;
    dev-8-sim-compare-freq_bs-batch.sh)
      echo "8-sim-compare-freq_bs-batch"
      ;;
    dev-11a-sim-low-separation-cyt-pos.sh)
      echo "11a-sim-low-separation-cyt-pos"
      ;;
    dev-11b-sim-low-separation-batch.sh)
      echo "11b-sim-low-separation-batch"
      ;;
    *)
      echo ""
      ;;
  esac
}

# QMDs rendered with plots (simulations off) once a launcher's jobs finish,
# separated by ':' so the list passes through --export unchanged.
plot_qmds_for_script() {
  case "$1" in
    dev-1-sim-trans.sh) echo "analysis/1-sim-trans.qmd" ;;
    dev-6-sim-tune-comparators.sh) echo "analysis/6-sim-tune-comparators.qmd" ;;
    dev-7-sim-compare-freq_bs.sh) echo "analysis/7-sim-compare-freq_bs.qmd" ;;
    dev-8-sim-compare-freq_bs-batch.sh)
      echo "analysis/8-sim-compare-freq_bs-batch.qmd"
      ;;
    dev-9a-real-tune-acs-cytof-tailgate.sh)
      echo "analysis/9a-real-tune-acs-cytof-tailgate.qmd"
      ;;
    # Analysis 10 only presents Analysis 9b's saved results.
    dev-9b-real-compare-acs-cytof.sh)
      echo "analysis/9b-real-compare-acs-cytof.qmd:analysis/10-real-compare-acs-cytof-validation.qmd"
      ;;
    dev-12-sim-cluster-gates.sh)
      echo "analysis/12-sim-cluster-gates.qmd"
      ;;
    dev-13-real-debug-acs-cytof.sh)
      echo "analysis/13-real-debug-acs-cytof.qmd"
      ;;
    # The OMIP analyses (14, 15) and 17 draw their plots in their
    # own job (plots_in_sim_job), so they have no plot job.
    dev-14-real-compare-omip111.sh | dev-15-real-compare-omip016.sh | \
      dev-17-explore-acs-cytof-coexpression.sh)
      ;;
    *)
      qmd_stem="$(chunked_qmd_stem_for_script "$1")"
      if [[ -n "$qmd_stem" ]]; then
        echo "analysis/${qmd_stem}.qmd"
      fi
      ;;
  esac
}

# Analyses quick enough to draw their plots in the simulation render, and
# Analysis 17, which only draws plots.
plots_in_sim_job() {
  case "$1" in
    dev-14-real-compare-omip111.sh | dev-15-real-compare-omip016.sh | \
      dev-17-explore-acs-cytof-coexpression.sh | dev-18-explore-acs-cytof-coexpression-gates.sh)
      return 0
      ;;
    *) return 1 ;;
  esac
}

# Simulation jobs that must wait for another launcher's jobs to succeed, as
# ':'-prefixed job IDs (empty when that launcher is not in this submission).
# Analysis 9b reads Analysis 9a's GatingSets and Tailgate selection. Analyses
# 13 and 18 read Analysis 9b's caches and results, and 17 reads its gates.
sim_dependency_for_script() {
  case "$1" in
    dev-9b-real-compare-acs-cytof.sh)
      echo "${script_job_ids[dev-9a-real-tune-acs-cytof-tailgate.sh]:-}"
      ;;
    dev-13-real-debug-acs-cytof.sh | dev-17-explore-acs-cytof-coexpression.sh | \
      dev-18-explore-acs-cytof-coexpression-gates.sh)
      echo "${script_job_ids[dev-9b-real-compare-acs-cytof.sh]:-}"
      ;;
    *) echo "" ;;
  esac
}

# Submit through slurm-sbatch, show its output and set `submitted_job_id`.
submit_job() {
  local out
  out="$(slurm-sbatch "$@")"
  printf '%s\n' "$out"
  submitted_job_id="$(
    printf '%s\n' "$out" | sed -n 's/.*Submitted batch job \([0-9][0-9]*\).*/\1/p' | tail -n 1
  )"
}

tmp_before="$(mktemp)"
trap 'rm -f "$tmp_before"' EXIT

if [[ "$install_stimgate" == "TRUE" ]]; then
  echo "Installing stimgate"
  slurm-sbatch "$install_script"

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
else
  echo "Skipping stimgate installation"
fi

echo "Submitting downstream jobs"
echo "SIM_SIZE: ${sim_size:-parameters.sim_size in _projr.yml}"

declare -A script_job_ids=()
all_sim_job_ids=()

for script in "${scripts[@]}"; do
  qmd_stem="$(chunked_qmd_stem_for_script "$script")"
  script_job_ids["$script"]=""

  if [[ -n "$qmd_stem" ]]; then
    echo "Submitting transactional chunks for $script"
    echo "ANALYSIS_RUN_ID: $analysis_run_id"

    script_n_chunks="$(n_chunks_for_script "$script")"
    for chunk_index in $(seq 1 "$script_n_chunks"); do
      log_dir="_tmp/log/sbatch/${qmd_stem}/chunk-${chunk_index}"
      job_name="${qmd_stem}-${chunk_index}"

      echo "Submitting $script chunk $chunk_index of $script_n_chunks"
      echo "Log directory: $log_dir"

      submit_job -l "$log_dir" -n "$script_dir/$script" -- \
        --job-name="$job_name" \
        --export=ALL,PROJECT_ROOT="$project_root",ANALYSIS_RUN_ID="$analysis_run_id",SIM_GRID_CHUNK_INDEX="$chunk_index",SIM_GRID_N_CHUNKS="$script_n_chunks",SIM_GRID_SHUFFLE_SEED="$sim_grid_shuffle_seed",RUN_SIMULATIONS=true,RUN_PLOTS=false"$sim_size_export"
      script_job_ids["$script"]+=":$submitted_job_id"
      all_sim_job_ids+=("$submitted_job_id")
    done
  else
    echo "Submitting $script"
    dependency_args=()
    dependency=""
    dependency_ids="$(sim_dependency_for_script "$script")"
    if [[ -n "$dependency_ids" ]]; then
      dependency="afterok${dependency_ids}"
    fi
    run_plots=false
    if plots_in_sim_job "$script"; then
      run_plots=true
      # Analysis 1's projr build clears projr's output folder, where figures go.
      projr_ids="${script_job_ids[dev-1-sim-trans.sh]:-}"
      if [[ -n "$projr_ids" ]]; then
        dependency+="${dependency:+,}afterany${projr_ids}"
      fi
    fi
    if [[ -n "$dependency" ]]; then
      dependency_args=(--dependency="$dependency")
    fi
    submit_job "$script_dir/$script" -- \
      ${dependency_args[@]+"${dependency_args[@]}"} \
      --export=ALL,PROJECT_ROOT="$project_root",ANALYSIS_RUN_ID="$analysis_run_id",RUN_SIMULATIONS=true,RUN_PLOTS="$run_plots""$sim_size_export"
    script_job_ids["$script"]=":$submitted_job_id"
    all_sim_job_ids+=("$submitted_job_id")
  fi
done

# Render each analysis's full report from its saved results (simulations off,
# plots on). A plot job starts after its own simulation jobs succeed and every
# other simulation job of this submission has finished: the projr builds clear
# projr's output folder, so plots must not be written while one may still run.
echo "Submitting plot renders"
for script in "${scripts[@]}"; do
  plot_qmds="$(plot_qmds_for_script "$script")"
  own_ids="${script_job_ids[$script]}"
  if [[ -z "$plot_qmds" ]]; then
    continue
  fi
  if [[ "$own_ids" == *"::"* || "$own_ids" == ":" || -z "$own_ids" ]]; then
    echo "WARNING: Could not read job IDs for $script; skipping its plot render." >&2
    continue
  fi
  other_ids=""
  for job_id in "${all_sim_job_ids[@]}"; do
    if [[ "$own_ids:" != *":$job_id:"* ]]; then
      other_ids+=":$job_id"
    fi
  done
  dependency="afterok${own_ids}"
  if [[ -n "$other_ids" ]]; then
    dependency+=",afterany${other_ids}"
  fi
  plot_stem="${script#dev-}"
  plot_stem="${plot_stem%.sh}"
  echo "Submitting plot render for $script ($plot_qmds)"
  submit_job -l "_tmp/log/sbatch/plots/${plot_stem}" -n "$script_dir/render-plots.sh" -- \
    --job-name="plots-${plot_stem}" \
    --dependency="$dependency" \
    --export=ALL,PROJECT_ROOT="$project_root",PLOT_QMD_FILES="$plot_qmds",ANALYSIS_RUN_ID="$analysis_run_id",RUN_SIMULATIONS=false,RUN_PLOTS=true"$sim_size_export"
done

echo "All downstream jobs submitted"
