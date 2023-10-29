#!/bin/bash

#PBS -P te53
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o aligned.log

# Load the required modules
module load samtools/1.12
module load bwa/0.7.17

# Directories and reference file
FASTQ_DIR="/g/data/te53/rdna/converted-illumin-fastq"
REFERENCE="/g/data/te53/rdna/reference/KY962517.1.fasta"
OUTPUT_DIR="/g/data/te53/kh3349/rDNA/rdnafastq-aligned-KY962517"
LIST_FILE="/g/data/te53/kh3349/rDNA/list.txt"  # Added this line

# Loop through the fastq files listed in list.txt and align them to the reference
while read -r fastq_file_name; do
    # Construct the input and output paths
    input_file="${FASTQ_DIR}/${fastq_file_name}"
    sam_output="${OUTPUT_DIR}/${fastq_file_name%.fastq.gz}.sam"
    bam_output="${OUTPUT_DIR}/${fastq_file_name%.fastq.gz}.bam"
    sorted_bam_output="${OUTPUT_DIR}/${fastq_file_name%.fastq.gz}.sorted.bam"

    # Align the fastq file to the reference and output as SAM
    bwa mem $REFERENCE $input_file > $sam_output

    # Convert SAM to BAM
    samtools view -Sb $sam_output > $bam_output

    # Sort the BAM file
    samtools sort $bam_output -o $sorted_bam_output

    # Optional: Remove the unsorted BAM file to save space
    rm $bam_output

    echo "Aligned and sorted ${fastq_file_name}. Output saved to ${sorted_bam_output}"

done < $LIST_FILE  # Modified this line to use the absolute path

echo "Alignment and sorting completed."
