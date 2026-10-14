#!/bin/bash
# Batch-shift reruns, at most 4 at a time: current rule (am_) and plateau step before the cap (pb_).
# Count only processes named R (not shells whose command text mentions the script).
P=$1; cd /scratch/01463748/projects/2024/stimgate/stimgate
nR() { ps -C R -o args= | grep -c "rerun[.]R"; }
MS='list(stim_mean_shift = 0.2, stim_mean_shift_clusters = "gn", mismatch_type = "mean_shift_negative")'
SD='list(stim_sd_multiplier = 1.1, stim_sd_multiplier_clusters = "gn", mismatch_type = "sd_inflation")'
for base in "low 0.002 2e4 g002_20k" "high 0.002 2e4 h002_20k" "low 0.002 1e5 g002"; do set -- $base
  for sh in ms02 sd01; do ov=$MS; [ $sh = sd01 ] && ov=$SD
    for v in am pb; do
      while [ $(nR) -ge 4 ]; do sleep 15; done
      pl=0; [ $v = pb ] && pl=1
      (P=$P TRANS=gaussian SEP=$1 PLATEAU=$pl nohup Rscript $P/rerun.R $2 $3 ${v}_$4_$sh "list()" "$ov" > $P/log_${v}_$4_$sh.txt 2>&1 &)
      sleep 5
    done; done; done
while [ $(nR) -gt 0 ]; do sleep 20; done
echo "not done:"; grep -L "^done" $P/log_*_ms02.txt $P/log_*_sd01.txt
echo batch-finished
