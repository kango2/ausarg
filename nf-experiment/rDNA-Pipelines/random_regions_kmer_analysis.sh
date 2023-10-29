#!/bin/bash
#PBS -P te53
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=48:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o randomregions.log

# Load necessary modules
module load samtools
module load parallel
module load jellyfish

# Define input directories and file paths
BAM_DIR="/g/data/te53/chm13_phase3_20220817/"
REFERENCE_BED="/g/data/te53/rdna/tmp/randomregions.100.bed"
REFERENCE_FASTA="/g/data/te53/ontsv/references/chm13_reference_dir/chm13.draft_v1.1.fasta" 

# Specify output directory for the randomregion files
output_dir="/g/data/te53/nj8315/randomregions"

# Number of threads for samtools
THREADS=4

# Define function to convert CRAM files specific to rDNA regions into FASTQ format and perform jellyfish analysis
randomregions() {
    cram_file="$1"
    # Extract the sample ID (directory name) from the CRAM file path
    SAMPLE_ID=$(basename $(dirname ${cram_file}))
    # Define BAM file path
    BAM="${output_dir}/${SAMPLE_ID}.bam"
    OUTPUT_FASTQ="${output_dir}/${SAMPLE_ID}.random.fastq"
    
    samtools view -b -@ $THREADS -o $BAM -T $REFERENCE_FASTA -L $REFERENCE_BED $cram_file
    samtools fastq -@ $THREADS -o $OUTPUT_FASTQ $BAM
    
    # Run jellyfish analysis for k-mer sizes 9, 11, 13, and 17
    for kmer_size in 9 11 13 17; do
    OUTPUT_JF="${output_dir}/${SAMPLE_ID}.k${kmer_size}.jf"
    OUTPUT_HISTO="${output_dir}/${SAMPLE_ID}.k${kmer_size}.histo"
    jellyfish count --min-quality=20 --quality-start=33 -s 100M -m $kmer_size -C -o $OUTPUT_JF $OUTPUT_FASTQ
    jellyfish histo $OUTPUT_JF > $OUTPUT_HISTO
done

    echo "rDNA regions from $cram_file have been extracted, saved as FASTQ file at $OUTPUT_FASTQ, and analyzed using jellyfish."
}

# Export the function and variables
export -f randomregions
export BAM_DIR REFERENCE_BED REFERENCE_FASTA output_dir THREADS

# Discover CRAM files in the directory
CRAM_FILES=$(find ${BAM_DIR} -type f -name "*.cram")

# Use GNU parallel to process multiple CRAM files concurrently
echo "$CRAM_FILES" | parallel -j 10 randomregions

echo "Conversion and analysis completed."
