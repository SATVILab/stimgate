#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --job-name="dev-1-trans"
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
apptainer-rscript -f stimgate -- 'projr::projr_build_dev(profile = "", file = "analysis/1-sim-trans.qmd")'
echo "Completed running projr"
date
echo "-------------------"
echo " "

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Convert duration to human-readable format
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

# Append the duration to the Slurm standard output log
echo "--- Script Duration ---"
printf "Elapsed time: %02d:%02d:%02d\n" $hours $minutes $seconds
echo "-----------------------"
