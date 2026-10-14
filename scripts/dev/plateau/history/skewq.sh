#!/bin/bash
# Skew high-separation reruns, at most 4 at a time: current rule (am_) and plateau step before the cap (pb_).
P=$1; cd /scratch/01463748/projects/2024/stimgate/stimgate
nR() { ps -C R -o args= | grep -c "rerun[.]R"; }
for base in "5e-4 1e5 s0005" "1e-4 1e5 s0001" "0.002 2e4 s002_20k" "5e-4 2e4 s0005_20k"; do set -- $base
  for v in am pb; do
    while [ $(nR) -ge 4 ]; do sleep 15; done
    pl=0; [ $v = pb ] && pl=1
    (P=$P TRANS=skew SEP=high PLATEAU=$pl nohup Rscript $P/rerun.R $1 $2 ${v}_$3 > $P/log_${v}_$3.txt 2>&1 &)
    sleep 5
  done; done
while [ $(nR) -gt 0 ]; do sleep 20; done
echo "not done:"; grep -L "^done" $P/log_*_s0005.txt $P/log_*_s0001.txt $P/log_*_s002_20k.txt $P/log_*_s0005_20k.txt
echo skew-finished
