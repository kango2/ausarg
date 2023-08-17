#!/bin/bash

# Get the paths from command-line arguments
CHUNKS_DIR="$1"
CSV_DIR="$2"

# Check if arguments are provided
if [ -z "$CHUNKS_DIR" ] || [ -z "$CSV_DIR" ]; then
    echo "Usage: $0 <path_to_chunks_directory> <path_to_csv_directory>"
    exit 1
fi

# Iterate through each chunk in the directory
for CHUNK in "$CHUNKS_DIR"/*.fq; do

    # Create an individual PBS script for the chunk
    JOB_FILE="job_$$.pbs"
    cat > "$JOB_FILE" <<EOL
#!/bin/bash
#PBS -N quality_value
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=0:05:00
#PBS -l mem=8GB
#PBS -l ncpus=4
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -o /g/data/xl04/ka6418/logs
#PBS -e /g/data/xl04/ka6418/logs
#PBS -M kirat.alreja@anu.edu.au

# Change to the directory the script is run from

module load biopython # Load necessary modules. Adjust as necessary

python3 /g/data/xl04/ka6418/ausarg/temp/qv_frequency_concurrent/quality.py "$CHUNK" "$CSV_DIR"

EOL

    # Submit the PBS script
    qsub "$JOB_FILE"

    # Optional: Remove the job file after submission
    # rm "$JOB_FILE"

done
