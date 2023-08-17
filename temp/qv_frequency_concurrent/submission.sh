#!/bin/bash

# Directory containing your FASTQ chunk files
CHUNKS_DIR="/g/data/xl04/ka6418/QCflow/qv_concurrent/chunks"

# Iterate through each chunk in the directory
for CHUNK in $CHUNKS_DIR/*.gz; do

    # Create an individual PBS script for the chunk
    JOB_FILE="job_$$.pbs"
    cat > $JOB_FILE <<EOL
#!/bin/bash
#PBS -N quality_value
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=0:30:00
#PBS -l mem=32GB
#PBS -l ncpus=16
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

# Change to the directory the script is run from

module load biopython # Load necessary modules. Adjust as necessary

python3 /g/data/xl04/ka6418/QCflow/qv_concurrent/quality.py $CHUNK /g/data/xl04/ka6418/QCflow/qv_concurrent/csvs

EOL

    # Submit the PBS script
    qsub $JOB_FILE

    # Optional: Remove the job file after submission
    # rm $JOB_FILE

done
