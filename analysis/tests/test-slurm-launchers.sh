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
n=$(( $(cat "$SLURM_TEST_LOG.n" 2>/dev/null || echo 100) + 1 ))
echo "$n" > "$SLURM_TEST_LOG.n"
echo "Submitted batch job $n"
EOF
cat > "$test_dir/bin/apptainer-rscript" <<'EOF'
#!/usr/bin/env bash
printf '%s|' "$@" >> "$SLURM_TEST_LOG"
printf 'run=%s|chunk=%s|n_chunks=%s\n' "$ANALYSIS_RUN_ID" "$SIM_GRID_CHUNK_INDEX" "$SIM_GRID_N_CHUNKS" >> "$SLURM_TEST_LOG"
# Record whether the QMD being rendered exists at render time.
qmd=$(printf '%s' "${!#}" | sed -n "s/^qmd_file <- '\([^']*\)'.*/\1/p")
if [[ -n "$qmd" && -f "$qmd" ]]; then
  printf 'render_file_exists=%s\n' "$qmd" >> "$SLURM_TEST_LOG"
fi
EOF
chmod +x "$test_dir/bin/"*
export PATH="$test_dir/bin:$PATH"
export PROJECT_ROOT="$project_root"
export SIM_GRID_N_CHUNKS=2 ANALYSIS_RUN_ID=slurm-test-run
export SIM_GRID_SHUFFLE_SEED=123
unset SIM_GRID_QMD_FILE SIM_GRID_CHUNK_INDEX SLURM_ARRAY_TASK_ID

run_selection() {
  : > "$SLURM_TEST_LOG"
  echo 100 > "$SLURM_TEST_LOG.n"
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
  # Two chunk jobs and one plot render.
  [[ $(wc -l < "$SLURM_TEST_LOG") -eq 3 ]]
  [[ $(grep -c "dev-${analysis_id}-" "$SLURM_TEST_LOG") -eq 2 ]]
  grep -Fq -- 'ANALYSIS_RUN_ID=slurm-test-run' "$SLURM_TEST_LOG"
  grep -Fq -- 'SIM_GRID_CHUNK_INDEX=1,SIM_GRID_N_CHUNKS=2' "$SLURM_TEST_LOG"
  grep -Fq -- 'SIM_GRID_CHUNK_INDEX=2,SIM_GRID_N_CHUNKS=2' "$SLURM_TEST_LOG"
done

run_selection 2a 2b 2a
[[ $(wc -l < "$SLURM_TEST_LOG") -eq 6 ]]
[[ $(grep -c 'dev-2a-' "$SLURM_TEST_LOG") -eq 2 ]]
[[ $(grep -c 'dev-2b-' "$SLURM_TEST_LOG") -eq 2 ]]
run_selection dev-2b-stim-bias_uns-freq_bs.sh
[[ $(wc -l < "$SLURM_TEST_LOG") -eq 3 ]]

# Plot renders wait for their own simulation jobs to succeed and for every
# other simulation job of the submission to finish.
run_selection 2a 7
plot_2a=$(grep -F 'render-plots.sh' "$SLURM_TEST_LOG" | grep -F 'plots-2a-')
plot_7=$(grep -F 'render-plots.sh' "$SLURM_TEST_LOG" | grep -F 'plots-7-')
[[ "$plot_2a" == *"--dependency=afterok:101:102,afterany:103|"* ]]
[[ "$plot_7" == *"--dependency=afterok:103,afterany:101:102|"* ]]
[[ "$plot_2a" == *"PLOT_QMD_FILES=analysis/2a-sim-bw-freq_bs-global.qmd,RUN_SIMULATIONS=false,RUN_PLOTS=true"* ]]
[[ "$plot_7" == *"PLOT_QMD_FILES=analysis/7-sim-compare-freq_bs.qmd,"* ]]
run_selection 9
grep -Fq -- 'PLOT_QMD_FILES=analysis/9-real-compare-acs-cytof.qmd:analysis/10-real-compare-acs-cytof-validation.qmd,' "$SLURM_TEST_LOG"

# The plot launcher renders each listed QMD itself with simulations off.
: > "$SLURM_TEST_LOG"
PLOT_QMD_FILES="analysis/9-real-compare-acs-cytof.qmd:analysis/10-real-compare-acs-cytof-validation.qmd" \
  bash "$project_root/scripts/slurm/render-plots.sh" > "$test_dir/output" 2>&1
grep -Fq -- "render_file_exists=analysis/9-real-compare-acs-cytof.qmd" "$SLURM_TEST_LOG"
grep -Fq -- "render_file_exists=analysis/10-real-compare-acs-cytof-validation.qmd" "$SLURM_TEST_LOG"

: > "$SLURM_TEST_LOG"
if bash "$project_root/scripts/slurm/dev.sh" 2a missing > "$test_dir/output" 2>&1; then
  echo 'Unknown selections must fail.' >&2
  exit 1
fi
[[ ! -s "$SLURM_TEST_LOG" ]]
grep -Fq -- 'Unknown analysis target: missing' "$test_dir/output"

for stem in 2a-stim-bw-freq_bs-global 2b-stim-bias_uns-freq_bs 3-sim-bw-est-base 4-sim-bw-est-norm; do
  : > "$SLURM_TEST_LOG"
  bash "$project_root/scripts/slurm/dev-${stem}.sh" 2 > "$test_dir/output" 2>&1
  # Launcher filenames use 'stim'; the corresponding QMD uses 'sim'.
  qmd_stem="${stem/-stim-/-sim-}"
  # Each chunk renders its own copy of the QMD, which exists during the
  # render and is removed afterwards.
  grep -Fq -- "render_file_exists=analysis/${qmd_stem}--chunk2-job" "$SLURM_TEST_LOG"
  grep -Fq -- 'run=slurm-test-run|chunk=2|n_chunks=2' "$SLURM_TEST_LOG"
  if compgen -G "$project_root/analysis/${qmd_stem}--chunk*" > /dev/null; then
    echo "Temporary render copy of ${qmd_stem} was not removed." >&2
    exit 1
  fi
done
echo 'Slurm 2a/2b selection and render contracts passed.'
