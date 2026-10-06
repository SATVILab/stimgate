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
printf 'run=%s|chunk=%s|n_chunks=%s\n' "${ANALYSIS_RUN_ID:-}" "${SIM_GRID_CHUNK_INDEX:-}" "${SIM_GRID_N_CHUNKS:-}" >> "$SLURM_TEST_LOG"
printf 'expected_run=%s|preprocessing=%s|stimgate=%s|comparators=%s|simulations=%s|plots=%s\n' "${ANALYSIS_EXPECTED_RUN_ID:-}" "${RUN_PREPROCESSING:-}" "${RUN_STIMGATE:-}" "${RUN_COMPARATORS:-}" "${RUN_SIMULATIONS:-}" "${RUN_PLOTS:-}" >> "$SLURM_TEST_LOG"
printf 'sim_size=%s|mcse=%s\n' "${SIM_SIZE:-}" "${SHOW_MCSE:-}" >> "$SLURM_TEST_LOG"
# Record whether the QMD being rendered exists at render time.
qmd=$(printf '%s' "${!#}" | sed -n "s/^qmd_file <- '\([^']*\)'.*/\1/p")
if [[ -n "$qmd" && -f "$qmd" ]]; then
  printf 'render_file_exists=%s\n' "$qmd" >> "$SLURM_TEST_LOG"
fi
if [[ "${MOCK_RENDER_ARTIFACTS:-false}" == true ]]; then
  printf '%s\n' "${SHOW_MCSE:-none}" > "$(dirname -- "$qmd")/$(basename -- "$qmd" .qmd).html"
fi
exit "${MOCK_RENDER_STATUS:-0}"
EOF
chmod +x "$test_dir/bin/"*
export PATH="$test_dir/bin:$PATH"
export PROJECT_ROOT="$project_root"
export SIM_GRID_N_CHUNKS=2 ANALYSIS_RUN_ID=slurm-test-run
export SIM_GRID_SHUFFLE_SEED=123
unset SIM_GRID_QMD_FILE SIM_GRID_CHUNK_INDEX SLURM_ARRAY_TASK_ID SIM_SIZE

run_selection() {
  : > "$SLURM_TEST_LOG"
  echo 100 > "$SLURM_TEST_LOG.n"
  bash "$project_root/scripts/slurm/dev.sh" "$@" > "$test_dir/output" 2>&1
}

run_selection
[[ $(grep -c 'dev-2b-' "$SLURM_TEST_LOG") -eq 2 ]]
for chunk in 1 2; do
  bias_job=$(grep -F "2b-sim-bias_uns-freq_bs/chunk-${chunk}|" "$SLURM_TEST_LOG")
  [[ "$bias_job" == *"PROJECT_ROOT=$project_root,ANALYSIS_RUN_ID=slurm-test-run,SIM_GRID_CHUNK_INDEX=$chunk,SIM_GRID_N_CHUNKS=2,SIM_GRID_SHUFFLE_SEED=123,RUN_SIMULATIONS=true,RUN_PLOTS=false|"* ]]
done
# No override: renderers consult the active project parameter.
! grep -Fq -- 'SIM_SIZE=' "$SLURM_TEST_LOG"
grep -Fq -- 'SIM_SIZE: parameters.sim_size in _projr.yml' "$test_dir/output"

for analysis_id in 2a 2b 7 8; do
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

# Reject invalid chunk counts before submitting even a nonchunked analysis.
for invalid_chunks in 0 invalid 1.5; do
  : > "$SLURM_TEST_LOG"
  if SIM_GRID_N_CHUNKS="$invalid_chunks" run_selection 1 2a; then
    echo "Accepted invalid chunk count: $invalid_chunks" >&2
    exit 1
  fi
  [[ ! -s "$SLURM_TEST_LOG" ]]
  grep -Fq 'SIM_GRID_N_CHUNKS must be a positive integer' "$test_dir/output"
done

# Plot renders wait for their own simulation jobs to succeed and for every
# other simulation job of the submission to finish.
run_selection 2a 7
plot_2a=$(grep -F 'render-plots.sh' "$SLURM_TEST_LOG" | grep -F 'plots-2a-')
plot_7=$(grep -F 'render-plots.sh' "$SLURM_TEST_LOG" | grep -F 'plots-7-')
[[ "$plot_2a" == *"--dependency=afterok:101:102,afterany:103:104|"* ]]
[[ "$plot_7" == *"--dependency=afterok:103:104,afterany:101:102|"* ]]
[[ "$plot_2a" == *"PLOT_QMD_FILES=analysis/2a-sim-bw-freq_bs-global.qmd,ANALYSIS_RUN_ID=slurm-test-run,RUN_SIMULATIONS=false,RUN_PLOTS=true"* ]]
[[ "$plot_7" == *"PLOT_QMD_FILES=analysis/7-sim-compare-freq_bs.qmd,"* ]]
[[ "$plot_2a" == *",RUN_PLOTS=true|"* ]]

# Explicit final still reaches every simulation and plot job.
export SIM_SIZE=final
run_selection 2a 7
grep -Fq -- 'SIM_SIZE: final' "$test_dir/output"
[[ $(grep -c ',RUN_SIMULATIONS=true,RUN_PLOTS=false,SIM_SIZE=final|' "$SLURM_TEST_LOG") -eq 4 ]]
[[ $(grep -c ',RUN_SIMULATIONS=false,RUN_PLOTS=true,SIM_SIZE=final|' "$SLURM_TEST_LOG") -eq 2 ]]

# SIM_SIZE=draft reaches chunk, projr and plot jobs; invalid sizes submit nothing.
export SIM_SIZE=draft
run_selection 2a 7
grep -Fq -- 'SIM_SIZE: draft' "$test_dir/output"
[[ $(grep -c ',RUN_SIMULATIONS=true,RUN_PLOTS=false,SIM_SIZE=draft|' "$SLURM_TEST_LOG") -eq 4 ]]
[[ $(grep -c ',RUN_SIMULATIONS=false,RUN_PLOTS=true,SIM_SIZE=draft|' "$SLURM_TEST_LOG") -eq 2 ]]
: > "$SLURM_TEST_LOG"
bash "$project_root/scripts/slurm/dev-2a-stim-bw-freq_bs-global.sh" 2 > "$test_dir/output" 2>&1
grep -Fq -- 'sim_size=draft' "$SLURM_TEST_LOG"
grep -Fq -- 'SIM_SIZE: draft' "$test_dir/output"
for launcher in dev-1-sim-trans.sh dev-7-sim-compare-freq_bs.sh dev-8-sim-compare-freq_bs-batch.sh; do
  : > "$SLURM_TEST_LOG"
  bash "$project_root/scripts/slurm/$launcher" > "$test_dir/output" 2>&1
  grep -Fq -- 'sim_size=draft' "$SLURM_TEST_LOG"
  grep -Fq -- 'SIM_SIZE: draft' "$test_dir/output"
done
: > "$SLURM_TEST_LOG"
# A temporary fixture, so the renamed mock reports stay out of analysis/.
mkdir -p "$test_dir/draft"
printf '%s\n' '---' 'params:' '  show_mcse: on' '---' > "$test_dir/draft/fixture.qmd"
MOCK_RENDER_ARTIFACTS=true PLOT_QMD_FILES="$test_dir/draft/fixture.qmd" \
  bash "$project_root/scripts/slurm/render-plots.sh" > "$test_dir/output" 2>&1
grep -Fq -- 'sim_size=draft' "$SLURM_TEST_LOG"
grep -Fq -- 'SIM_SIZE: draft' "$test_dir/output"
export SIM_SIZE=huge
: > "$SLURM_TEST_LOG"
if run_selection 2a; then
  echo 'Invalid SIM_SIZE values must fail.' >&2
  exit 1
fi
[[ ! -s "$SLURM_TEST_LOG" ]]
grep -Fq -- 'SIM_SIZE must be final or draft' "$test_dir/output"
unset SIM_SIZE

run_selection 9
grep -Fq -- 'PLOT_QMD_FILES=analysis/9-real-compare-acs-cytof.qmd:analysis/10-real-compare-acs-cytof-validation.qmd,' "$SLURM_TEST_LOG"

# Non-chunked simulation jobs share the plot job's submission run ID.
run_selection 1 9
[[ $(grep -c 'ANALYSIS_RUN_ID=slurm-test-run' "$SLURM_TEST_LOG") -eq 4 ]]

# The plot launcher overrides inherited ACS stage controls.
export RUN_PREPROCESSING=true RUN_STIMGATE=true RUN_COMPARATORS=true
# The plot launcher renders each listed QMD itself with simulations off.
: > "$SLURM_TEST_LOG"
PLOT_QMD_FILES="analysis/9-real-compare-acs-cytof.qmd:analysis/10-real-compare-acs-cytof-validation.qmd" \
  bash "$project_root/scripts/slurm/render-plots.sh" > "$test_dir/output" 2>&1
grep -Fq -- "render_file_exists=analysis/9-real-compare-acs-cytof.qmd" "$SLURM_TEST_LOG"
grep -Fq -- "render_file_exists=analysis/10-real-compare-acs-cytof-validation.qmd" "$SLURM_TEST_LOG"
grep -Fq -- 'expected_run=slurm-test-run|preprocessing=false|stimgate=false|comparators=false|simulations=false|plots=true' "$SLURM_TEST_LOG"
[[ $(grep -c 'render_file_exists=' "$SLURM_TEST_LOG") -eq 2 ]]
# QMDs 9 and 10 have no Monte Carlo intervals, so each renders once.
[[ $(grep -c 'mcse=off' "$SLURM_TEST_LOG") -eq 0 ]]
[[ $(grep -c 'mcse=on' "$SLURM_TEST_LOG") -eq 0 ]]
[[ $(grep -c 'sim_size=.*|mcse=$' "$SLURM_TEST_LOG") -eq 2 ]]
# Reports are renamed after rendering; --output breaks embed-resources.
! grep -Fq -- '--output' "$SLURM_TEST_LOG"

unset RUN_PREPROCESSING RUN_STIMGATE RUN_COMPARATORS

# Persistent output names survive the second render (the isolated chunk helper
# must not be used for reports, since its exit trap removes HTML outputs).
mkdir "$test_dir/reports"
printf '%s\n' '---' 'format: html' 'params:' '  show_mcse: on' '---' > "$test_dir/reports/fixture.qmd"
printf '%s\n' '---' 'format: html' '---' > "$test_dir/reports/plain.qmd"
MOCK_RENDER_ARTIFACTS=true PLOT_QMD_FILES="$test_dir/reports/fixture.qmd:$test_dir/reports/plain.qmd" \
  bash "$project_root/scripts/slurm/render-plots.sh" > "$test_dir/output" 2>&1
[[ $(cat "$test_dir/reports/fixture-mcse_off.html") == off ]]
[[ $(cat "$test_dir/reports/fixture-mcse_on.html") == on ]]
[[ $(cat "$test_dir/reports/plain.html") == none ]]
[[ ! -e "$test_dir/reports/plain-mcse_off.html" ]]


: > "$SLURM_TEST_LOG"
if bash "$project_root/scripts/slurm/dev.sh" 2a missing > "$test_dir/output" 2>&1; then
  echo 'Unknown selections must fail.' >&2
  exit 1
fi
[[ ! -s "$SLURM_TEST_LOG" ]]
grep -Fq -- 'Unknown analysis target: missing' "$test_dir/output"

for stem in 2a-stim-bw-freq_bs-global 2b-stim-bias_uns-freq_bs 3-sim-bw-est-base 4-sim-bw-est-norm 7-sim-compare-freq_bs 8-sim-compare-freq_bs-batch; do
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
# Every simulation launcher and the plot launcher propagate renderer failures.
export MOCK_RENDER_STATUS=23
for launcher_path in "$project_root"/scripts/slurm/dev-*.sh; do
  if bash "$launcher_path" > "$test_dir/output" 2>&1; then
    echo "Renderer failure must fail $launcher_path." >&2
    exit 1
  else
    [[ $? -eq 23 ]]
  fi
  grep -Fq -- 'SIM_SIZE:' "$test_dir/output" || [[ "$launcher_path" == *dev-9-* ]]
done
if PLOT_QMD_FILES="analysis/1-sim-trans.qmd" bash "$project_root/scripts/slurm/render-plots.sh" > "$test_dir/output" 2>&1; then
  echo 'Renderer failure must fail the plot launcher.' >&2
  exit 1
else
  [[ $? -eq 23 ]]
fi
unset MOCK_RENDER_STATUS
echo 'Slurm selection, run provenance, stage controls and render failure contracts passed.'
