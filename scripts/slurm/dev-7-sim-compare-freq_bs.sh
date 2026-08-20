#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --job-name="dev-7-compare"
#SBATCH --partition=ada

# Record the start time
start_time=$(date +%s)

echo "HOSTNAME: $HOSTNAME"

echo " "
echo " "
echo " "

echo "-------------------"
echo "Run projr"
date

apptainer-rscript -f stimgate -- \
  'projr::projr_build_dev(profile = "", file = "analysis/7-sim-compare-freq_bs.qmd")'

echo "Completed running projr"
date
echo "-------------------"
echo " "

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Convert the duration to human-readable format
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))

echo "--- Script Duration ---"
printf "Elapsed time: %02d:%02d:%02d\n" \
  "$hours" "$minutes" "$seconds"
echo "-----------------------"
