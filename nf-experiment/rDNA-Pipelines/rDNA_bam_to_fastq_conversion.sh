#!/bin/bash
#PBS -P te53
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o conversion.log

# Load necessary modules
module load samtools
module load bedtools
module load parallel

# Define input directories and file paths
BAM_DIR="/g/data/te53/ontsv/data/alignments/minimap2/v2.22-r1101"
REFERENCE_BED="/g/data/te53/ontsv/hrp561/scripts/chm13-t2t-rdnamodels.ucsc.20230611.bed"
# Specify output directory for the converted files
output_dir="/Path/to/Bam-to-Fastq-conversion/"

# Define function to convert BAM files specific to rDNA regions into FASTQ format
convert_bam_to_fastq_rDNA() {
    bam_file="$1"
    # Extract the RunID from the BAM file path
    RUN_ID=$(basename $(dirname $(dirname ${bam_file})))
    # Define temporary file paths
    TMP_BAM="${output_dir}/${RUN_ID}.tmp.bam"
    TMP_SORTED_BAM="${output_dir}/${RUN_ID}.tmp_sorted.bam"
    OUTPUT_FASTQ="${output_dir}/${RUN_ID}.fastq"
    # Extract rDNA regions from the BAM file using bedtools
    bedtools intersect -abam $bam_file -b $REFERENCE_BED -wa > $TMP_BAM
    # Sort the temporary BAM file
    samtools sort $TMP_BAM -o $TMP_SORTED_BAM
    # Convert the sorted BAM file to FASTQ format
    samtools fastq $TMP_SORTED_BAM > $OUTPUT_FASTQ
    # Cleanup: remove temporary files
    rm $TMP_BAM $TMP_SORTED_BAM
    echo "rDNA regions from $bam_file have been extracted and saved as FASTQ file at $OUTPUT_FASTQ"
}

# Export relevant variables and functions for parallel processing
export -f convert_bam_to_fastq_rDNA
export BAM_DIR
export REFERENCE_BED
export output_dir

# List of specific subdirectory patterns to look for BAM files
SUBDIRS=("PBXP*" "PBYP*" "PGXX*" "PLPN*" "PUXP*")

# Loop through each subdirectory pattern
for SUBDIR in "${SUBDIRS[@]}"; do
  # Discover BAM files in the subdirectories
  BAM_FILES=$(find ${BAM_DIR} -type f -path "${BAM_DIR}/${SUBDIR}/*/*_pass.bam")
  # Use GNU parallel to process multiple BAM files concurrently
  parallel -j 4 convert_bam_to_fastq_rDNA ::: ${BAM_FILES}
done
echo "Conversion completed."
