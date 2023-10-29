#!/bin/bash

#PBS -P te53
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o matching.log

# Directory containing .fastq.gz files
FASTQ_DIR="/g/data/te53/rdna/Bam-to-Fastq-conversion"

# Directory containing .done files
DONE_DIR="/g/data/te53/ontsv/hrp561/cnv"

# File to save the matching names
OUTPUT_FILE="/g/data/te53/nj8315/matching-names.txt"

# Debugging logs
DEBUG_LOG="/g/data/te53/nj8315/debug.log"
> "$DEBUG_LOG"

# Clear the contents of the output file (if it exists) before writing new matches
> "$OUTPUT_FILE"

for fastq in $FASTQ_DIR/*.fastq.gz; do
    # Extract base name from the fastq file
    base_name=$(basename "$fastq" .fastq.gz)

    # Construct the .done file name
    done_file="${DONE_DIR}/${base_name}_pass.1kbp.bedcov.done"

    # Debugging: Log what is being checked
    echo "Checking for: $done_file" >> "$DEBUG_LOG"

    # Check if the .done file exists
    if [[ -f "$done_file" ]]; then
        echo "$fastq has been previously analyzed."
        echo "$base_name" >> "$OUTPUT_FILE"
    else
        echo "$fastq has NOT been previously analyzed."
        # Debugging: Log the names that don't match
        echo "No match for: $fastq" >> "$DEBUG_LOG"
    fi
done
