#!/bin/bash
#PBS -P te53
#PBS -l ncpus=240
#PBS -l mem=190GB
#PBS -l walltime=24:00:00
#PBS -l iointensive=10
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o analysis.log

# Load necessary modules
module load samtools 

# Define paths and directories
REFERENCE="/g/data/te53/rdna/reference/KY962517.1.fasta"
BAM_DIR="/g/data/te53/nj8315/rdnafastq-aligned-KY962517"
DEPTH_DIR="${BAM_DIR}/depth_files"
STATS_DIR="${BAM_DIR}/alignment_stats"

mkdir -p $DEPTH_DIR $STATS_DIR

# Loop through each .sorted.bam file in the directory
for BAM_FILE in ${BAM_DIR}/*sorted.bam; do
    BASE_NAME=$(basename $BAM_FILE .sorted.bam)
    
    # Extract alignment scores, mismatches, and read length
    echo "Processing $BAM_FILE..."
    echo "Extracting alignment scores, mismatches, and read length for $BASE_NAME..."
    samtools view $BAM_FILE | 
    awk '{for(i=12;i<=NF;++i) {if ($i ~ /^AS:i:/) {score=$i} else if ($i ~ /^NM:i:/) {mismatches=$i}}; print score, mismatches, length($10)}' > $STATS_DIR/${BASE_NAME}_alignment_info.txt

    # Compute the coverage for each position
    echo "Computing coverage for $BASE_NAME..."
    DEPTH_FILE="${DEPTH_DIR}/${BASE_NAME}_depth.txt"
    samtools depth $BAM_FILE > $DEPTH_FILE

    # Extract start positions of all reads
    echo "Extracting start positions of reads for $BASE_NAME..."
    samtools view $BAM_FILE | awk '{print $4}' > $STATS_DIR/${BASE_NAME}_start_positions.txt
done

echo "Analysis completed."

