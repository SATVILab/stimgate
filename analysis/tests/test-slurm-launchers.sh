#!/usr/bin/env bash
# Exercise launch selection and rendering without submitting jobs or running R.
set -euo pipefail
project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT
mkdir "$test_dir/bin"
export SLURM_TEST_LOG="$test_dir/calls"
cat > "$test_dir/bin/slurm-sbatch" <<'EOF'
#!/usr/bin/env bash
printf '%s|' "$@" >> "$SLURM_TEST_LOG"
printf '\n' >> "$SLURM_TEST_LOG"
EOF
cat > "$test_dir/bin/apptainer-rscript" <<'EOF'
#!/usr/bin/env bash
printf '%s|' "$@" >> "$SLURM_TEST_LOG"
printf 'run=%s|chunk=%s|n_chunks=%s\n' "$ANALYSIS_RUN_ID" "$SIM_GRID_CHUNK_INDEX" "$SIM_GRID_N_CHUNKS" >> "$SLURM_TEST_LOG"
EOF
chmod +x "$test_dir/bin/"*
export PATH="$test_dir/bin:$PATH"
export PROJECT_ROOT="$project_root"
export SIM_GRID_N_CHUNKS=2 ANALYSIS_RUN_ID=slurm-test-run
export SIM_GRID_SHUFFLE_SEED=123
unset SIM_GRID_QMD_FILE SIM_GRID_CHUNK_INDEX SLURM_ARRAY_TASK_ID

run_selection() {
  : > "$SLURM_TEST_LOG"
  bash "$project_root/scripts/slurm/dev.sh" "$@" > "$test_dir/output" 2>&1
}

run_selection
[[ $(grep -c 'dev-2b-' "$SLURM_TEST_LOG") -eq 2 ]]
for chunk in 1 2; do
  bias_job=$(grep -F "2b-sim-bias_uns-freq_bs/chunk-${chunk}|" "$SLURM_TEST_LOG")
  [[ "$bias_job" == *"PROJECT_ROOT=$project_root,ANALYSIS_RUN_ID=slurm-test-run,SIM_GRID_CHUNK_INDEX=$chunk,SIM_GRID_N_CHUNKS=2,SIM_GRID_SHUFFLE_SEED=123,RUN_SIMULATIONS=true,RUN_PLOTS=false"* ]]
done

for analysis_id in 2a 2b; do
  run_selection "$analysis_id"
  [[ $(wc -l < "$SLURM_TEST_LOG") -eq 2 ]]
  [[ $(grep -c "dev-${analysis_id}-" "$SLURM_TEST_LOG") -eq 2 ]]
  grep -Fq -- 'ANALYSIS_RUN_ID=slurm-test-run' "$SLURM_TEST_LOG"
  grep -Fq -- 'SIM_GRID_CHUNK_INDEX=1,SIM_GRID_N_CHUNKS=2' "$SLURM_TEST_LOG"
  grep -Fq -- 'SIM_GRID_CHUNK_INDEX=2,SIM_GRID_N_CHUNKS=2' "$SLURM_TEST_LOG"
done

run_selection 2a 2b 2a
[[ $(wc -l < "$SLURM_TEST_LOG") -eq 4 ]]
[[ $(grep -c 'dev-2a-' "$SLURM_TEST_LOG") -eq 2 ]]
[[ $(grep -c 'dev-2b-' "$SLURM_TEST_LOG") -eq 2 ]]
run_selection dev-2b-stim-bias_uns-freq_bs.sh
[[ $(wc -l < "$SLURM_TEST_LOG") -eq 2 ]]

: > "$SLURM_TEST_LOG"
if bash "$project_root/scripts/slurm/dev.sh" 2a missing > "$test_dir/output" 2>&1; then
  echo 'Unknown selections must fail.' >&2
  exit 1
fi
[[ ! -s "$SLURM_TEST_LOG" ]]
grep -Fq -- 'Unknown analysis target: missing' "$test_dir/output"

for stem in 2a-stim-bw-freq_bs-global 2b-stim-bias_uns-freq_bs; do
  : > "$SLURM_TEST_LOG"
  bash "$project_root/scripts/slurm/dev-${stem}.sh" 2 > "$test_dir/output" 2>&1
  # Launcher filenames use 'stim'; the corresponding QMD uses 'sim'.
  qmd_stem="${stem/-stim-/-sim-}"
  grep -Fq -- "analysis/${qmd_stem}.qmd" "$SLURM_TEST_LOG"
  grep -Fq -- 'run=slurm-test-run|chunk=2|n_chunks=2' "$SLURM_TEST_LOG"
done
echo 'Slurm 2a/2b selection and render contracts passed.'
