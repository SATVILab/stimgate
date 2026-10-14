#!/bin/bash
# Run the Gaussian scenarios with the plateau step off, at most 4 R processes at a time.
P=$1; cd /scratch/01463748/projects/2024/stimgate/stimgate-plateau
nR() { ps -C R -o args= | grep -c "rerun[.]R"; }
for s in "high 1e-4 1e5" "high 5e-4 2e4" "high 5e-4 1e5" "high 2e-3 2e4" "high 2e-3 1e5" "low 1e-4 1e5" "low 5e-4 2e4" "low 5e-4 1e5" "low 2e-3 5e3" "low 2e-3 2e4" "low 2e-3 1e5" "low 1e-2 5e3" "low 1e-2 2e4" "low 1e-2 1e5"; do
  set -- $s; tag="off_${1}_${2}_${3}"
  while [ $(nR) -ge 4 ]; do sleep 15; done
  (P=$P TRANS=gaussian SEP=$1 nohup Rscript $P/rerun.R $2 $3 $tag "list(locPlateauStep = FALSE)" > $P/log_$tag.txt 2>&1 &); sleep 5
done
while [ $(nR) -gt 0 ]; do sleep 20; done
echo "not done:"; grep -L "^done" $P/log_off_*.txt
echo queue-finished
