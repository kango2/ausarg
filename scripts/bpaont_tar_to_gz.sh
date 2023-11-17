#!/bin/bash
#PBS -N fast5tar
#PBS -P xl04
#PBS -q hugemem
#PBS -l walltime=10:00:00
#PBS -l mem=32GB
#PBS -l ncpus=8
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=1400GB
#PBS -M kirat.alreja@anu.edu.au

# Exit script on any error
set -e

# Variables for input files
fast5pass=${passfile}
fast5fail=${failfile}

# Change to job specific file system and check for success
cd ${PBS_JOBFS} || { echo "Failed to change directory to ${PBS_JOBFS}"; exit 1; }

# Extract pass and fail files and check for success
tar -xvf ${fast5pass} || { echo "Failed to extract ${fast5pass}"; exit 1; }

# Process 'pass' files
cd *pass* || { echo "Failed to change directory to $(basename "$fast5pass" .tar)"; exit 1; }
cat *.gz > ${sample}_pass.gz || { echo "Failed to concatenate pass files"; exit 1; }
rsync -av ${sample}_pass.gz ${output} || { echo "Failed to rsync pass files"; exit 1; }

rm *.gz || { echo "Failed to delete pass files"; exit 1; }

# Process 'fail' files
cd ..
tar -xvf ${fast5fail} || { echo "Failed to extract ${fast5fail}"; exit 1; }
cd *fail* || { echo "Failed to change directory to $(basename "$fast5fail" .tar)"; exit 1; }
cat *.gz > ${sample}_fail.gz || { echo "Failed to concatenate fail files"; exit 1; }
rsync -av ${sample}_fail.gz ${output} || { echo "Failed to rsync fail files"; exit 1; }

rm *.gz || { echo "Failed to delete pass files"; exit 1; }

# End of script
echo "Processing complete."
