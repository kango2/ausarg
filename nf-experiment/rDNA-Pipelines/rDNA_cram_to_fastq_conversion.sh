#!/bin/bash
#PBS -P te53
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o conversion.log

# Load necessary modules
module load samtools
module load bedtools
module load parallel

# Define input directories and file paths
BAM_DIR="/g/data/te53/chm13_phase3_20220817/"
REFERENCE_BED="/g/data/te53/ontsv/hrp561/scripts/chm13-t2t-rdnamodels.ucsc.20230611.bed"
REFERENCE_FASTA="/g/data/te53/ontsv/references/chm13_reference_dir/chm13.draft_v1.1.fasta" 
# Specify output directory for the converted files
output_dir="/g/data/te53/rdna/converted-illumin-fastq"

# Number of threads for samtools
THREADS=4

# Define function to convert CRAM files specific to rDNA regions into FASTQ format
convert_cram_to_fastq_rDNA() {
    cram_file="$1"
    # Extract the sample ID (directory name) from the CRAM file path
    SAMPLE_ID=$(basename $(dirname ${cram_file}))
    # Define temporary file paths
    TMP_BAM="${output_dir}/${SAMPLE_ID}.tmp.bam"
    TMP_SORTED_BAM="${output_dir}/${SAMPLE_ID}.tmp_sorted.bam"
    OUTPUT_FASTQ="${output_dir}/${SAMPLE_ID}.fastq"
    # Convert CRAM to BAM format
    samtools view -b -@ $THREADS -T $REFERENCE_FASTA -o $TMP_BAM $cram_file
    # Extract rDNA regions from the BAM file using bedtools
    bedtools intersect -abam $TMP_BAM -b $REFERENCE_BED -wa > $TMP_BAM
    # Sort the temporary BAM file using threads
    samtools sort $TMP_BAM -o $TMP_SORTED_BAM -@ $THREADS
    # Convert the sorted BAM file to FASTQ format using threads
    samtools fastq -@ $THREADS $TMP_SORTED_BAM > $OUTPUT_FASTQ
    # Cleanup: remove temporary files
    rm $TMP_BAM $TMP_SORTED_BAM
    echo "rDNA regions from $cram_file have been extracted and saved as FASTQ file at $OUTPUT_FASTQ"
}

# Export the function and variables
export -f convert_cram_to_fastq_rDNA
export BAM_DIR REFERENCE_BED REFERENCE_FASTA output_dir THREADS

# Discover BAM files in the subdirectories
BAM_FILES=$(find ${BAM_DIR} -type f -name "*.cram")

# Use GNU parallel to process multiple BAM files concurrently
parallel -j 10 convert_cram_to_fastq_rDNA ::: ${BAM_FILES}

echo "Conversion completed."
