#!/bin/bash

# DESCRIPTION:
# This script performs the following operations:
# LiftoverVcf: This step transitions variants within the VCF files from hg38 to hg19 standards using Picard.
# To perform the liftover, you'll need a reference file matching the chromosome naming of your VCF files. 
# In our context, hg19.fa sourced from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) was the most suitable. 
# It's imperative to ensure hg19.fa.fai and hg19.fa.dict are also accessible in the same directory.
# Additionally, the chain file hg38ToHg19.over.chain is essential, and it's available for download at [USCS](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/).
# Subsetting, Filtering & Sorting: Managed through Bcftools, this step involves the viewing, subsetting, filtering, and sorting of genomic data.
# For the subsetting process, a list of SNPs is necessary. These SNPs correspond to markers listed on the SNP array. Fetch this file from [screening-array] (https://sapac.support.illumina.com/downloads/infinium-global-screening-array-v3-0-support-files.html).
# Authors: Kosar Hooshmand

#PBS -P te53
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o merge.log

# Load required modules
module load java bcftools

# Define input, output directory, and references
inputdir="/path/to/input_directory"
outputdir="/path/to/output_directory"
reference="/path/to/hg19.fa"
input_list="/path/to/SNP_list.txt"

#############
# LiftoverVcf step
#############
# Transition variants within the VCF files from hg38 to hg19 standards

input_files=("FS28687806.norm.vcf" "FS28687807.norm.vcf" "FS28687808.norm.vcf" "FS28687809.norm.vcf" "FS28687810.norm.vcf")

for file in "${input_files[@]}"
do
  java -jar -Xmx80g -XX:ParallelGCThreads=$PBS_NCPUS \
  /g/data/if89/apps/picard/2.27.4/picard.jar LiftoverVcf \
  I="${input_dir}/${file}" \
  O="${output_dir}/${file%.*}_liftover.vcf.gz" \
  CHAIN="/path/to/hg38ToHg19.over.chain" \
  REJECT="${output_dir}/${file%.*}_rejected_variants.vcf" \
  WARN_ON_MISSING_CONTIG=true R=${reference}
done

#############
# bcftools view step for subsetting
#############
# Subset the lifted over VCFs based on SNP list

input_files=("FS28687806.norm_liftover.vcf" "FS28687807.norm_liftover.vcf" "FS28687808.norm_liftover.vcf" "FS28687809.norm_liftover.vcf" "FS28687810.norm_liftover.vcf")

for file in "${input_files[@]}"
do
bcftools view --include ID==@"${input_list}" "${input_dir}/${file}" -O z -o ${output_dir}/${file%.*}_subset.vcf.gz 
done

#############
# Filtering and sorting step
#############
# Remove unwanted chromosomes and sort

input_files=("FS28687806.norm_liftover_subset.vcf" "FS28687807.norm_liftover_subset.vcf" "FS28687808.norm_liftover_subset.vcf" "FS28687809.norm_liftover_subset.vcf" "FS28687810.norm_liftover_subset.vcf")

for file in "${input_files[@]}"
do
    bcftools view "${input_dir}/${file}" | grep -v -e '^chrX' -e '^chrY' -e '^chrM' > "${output_dir}/${file%.*}_autosome.vcf" && \
    bcftools view "${output_dir}/${file%.*}_autosome.vcf" | grep -v | sort -k1,1V -k2,2n > "${output_dir}/${file%.*}_autosome_sorted.vcf"
done
