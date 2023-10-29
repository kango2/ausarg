#!/bin/bash

# DESCRIPTION:
# This script processes genomic data from raw FASTQ files, performs alignment, variant calling, and 
# finally converts the variants into the 23andMe format.
# Overview of the pipeline:
# Alignment & Indexing: Reads are aligned to the GRCh38 reference genome using BWA. The resulting BAM files are subsequently indexed.
# Post-processing: Implemented with the bwa-postalt.js JavaScript program located within the bwakit directory. Operated through k8, this process facilitates 
# variant calling across all alternate contigs, notably including leukocyte antigen (HLA) alternate contigs.
# Sorting & Deduplication: Sequences are sorted, and duplicates are marked using Samtools.
# Genotyping, Filtering & Annotation: This step involves the genotyping, subsequent filtering, and annotation of variants, all managed using Bcftools.
# 23andMe Conversion: VCF files are transformed into the 23andMe format, an operation combinedly performed by Bcftools, plink2, and plink1.9.
# Authors: Kosar Hooshmand

#PBS -j oe
#PBS -l ncpus=16
#PBS -l mem=63GB
#PBS -l walltime=24:00:00
#PBS -P te53
#PBS -l jobfs=390GB
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89

# Load required modules
module load bwa/0.7.17 samtools/1.12 k8/0.2.5 bcftools plink2

BWAKIT=/g/data/te53/software/bwa/0.7.17-r1198-dirty/bwa-postalt.js

# Define input, output directory, and values for RG headers
inputdir="/path/to/input_directory"
outputdir="/path/to/output_directory"
reference="/path/to/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.fa"
dbSNP="/path/to/dbsnp154/GCF_000001405.38.seqidmod.vcf.gz"
pl="illumina"
pm="hiseq"

# Error handling settings
set -e
set -o pipefail
set -u
set -x
set -o functrace

# Declare an associative array to hold unique sample IDs
declare -A unique_sampleIds

# Loop through each FASTQ file to extract unique sample IDs
for fastq in ${inputdir}/*R1.fq.gz; do
    sampleId=$(basename "$fastq" | awk -F'_' '{print $4}')
    unique_sampleIds["$sampleId"]=1
done

# Convert the associative array keys into an indexed array
sampleIds=("${!unique_sampleIds[@]}")

for sampleId in "${sampleIds[@]}"; do
    sm="${sampleId}"
    
    # Find the associated R1 and R2 files for this sampleId
    r1_files=($(ls "${inputdir}"/*_${sampleId}_*R1.fq.gz))
    r2_files=($(ls "${inputdir}"/*_${sampleId}_*R2.fq.gz))

    # Loop through each set of R1-R2 files for the given sampleId
    for index in "${!r1_files[@]}"; do
        r1="${r1_files[$index]}"
        r2="${r2_files[$index]}"

        bwa mem -t ${PBS_NCPUS} -K 100000000 -Y -R "@RG\\tID:${sampleId}\\tPL:${pl}\\tPM:${pm}\\tSM:${sm}" ${reference} ${r1} ${r2} \
        | k8 $BWAKIT -p ${outputdir}/${sampleId} ${reference}.alt \
        | samtools fixmate -m --threads ${PBS_NCPUS} --reference ${reference} - - \
        | samtools sort --threads ${PBS_NCPUS} -m 1536M --reference ${reference} -T ${PBS_JOBFS}/${sampleId}.sorttmp \
        | samtools markdup -T ${PBS_JOBFS}/${sampleId}.deduptmp --reference ${reference} --threads ${PBS_NCPUS} --write-index - ${outputdir}/${sampleId}.fmsortdedup.bam
    done
    
    # Merge BAM files for the current sample
    merged_bam_file="${outputdir}/${sampleId}_merged.bam"
    bam_files=($(ls "${outputdir}"/*_${sampleId}*.fmsortdedup.bam))
    samtools merge "${merged_bam_file}" "${bam_files[@]}"
    samtools index "${merged_bam_file}"
done

#############
# Begin Variant Calling
#############

# Create sample list file with full path to merged bam files
echo "Creating sample list file with full path to merged bam files..."
find "${outputdir}" -name "*_merged.bam" > "${outputdir}/merged_sample_list.txt"

# Define function to call variants for a given chromosome
call_variants () {
  chr=$1
  echo "Calling variants for chromosome ${chr}..."
  bcftools mpileup \
    -Ou -f "${reference}" -b "${outputdir}/merged_sample_list.txt" -r "${chr}" \
    --min-MQ 20 --min-BQ 20 --adjust-MQ 50 \
    --threads $PBS_NCPUS \
    --output-type b \
    --gvcf 10,20,30,60 \
    | bcftools call --ploidy 2 -vmO z -f GQ -o "${outputdir}/chr${chr}.g.vcf.gz"
  
  # Index the resulting g.vcf.gz file
  echo "Indexing the resulting g.vcf.gz file for chromosome ${chr}..."
  bcftools index "${outputdir}/chr${chr}.g.vcf.gz"
}

# Call variants from bam files for each chromosome separately in parallel
echo "Calling variants from bam files for each chromosome separately in parallel..."
for chr in $(seq 1 22) X Y M; do
  call_variants "${chr}" &
done

# Wait for all background processes to complete
wait

echo "Variant calling completed!"

#############
# Begin Post-processing
#############

# Concatenate produced VCF files
echo "Concatenating VCF files..."
input_list="${outputdir}/chr{1..22,X,Y,M}.g.vcf.gz"
bcftools concat -a -D -o "${outputdir}/bcftools_merged.vcf.gz" -Ov -f "${input_list}" --threads $PBS_NCPUS

# Index the merged VCF file
echo "Indexing the merged VCF file..."
bcftools index -t -f "${outputdir}/bcftools_merged.vcf.gz"

# Annotate SNPs
echo "Annotating SNPs..."
bcftools annotate --columns ID  --annotations ${dbSNP} --threads $PBS_NCPUS -Oz -o "${outputdir}/bcftools_annotate.vcf.gz" "${outputdir}/bcftools_merged.vcf.gz"

# Filtering SNPs
echo "Filtering SNPs..."
bcftools filter -e "QUAL<30 || FORMAT/DP>100" -s "LOW_QUAL" "${outputdir}/bcftools_annotate.vcf.gz" --threads $PBS_NCPUS -Oz -o "${outputdir}/bcftools_filtered.vcf.gz" 

# Extract the list of sample names from the VCF header
echo "Extracting sample names..."
samples=$(bcftools query -l "${outputdir}/bcftools_filtered.vcf.gz")

# Split the multi-sample VCF file into individual VCF files per sample
echo "Splitting VCF file per sample..."
for sample in $samples; do
    bcftools view -i "TYPE=='snp'" -f "PASS" -m2 -M2 --threads $PBS_NCPUS -Oz -s "$sample" -o "${outputdir}/${sample}.vcf.gz" "${outputdir}/bcftools_filtered.vcf.gz"
done

echo "Post-processing completed!"

#############
# Begin conversion to 23andMe
#############

for sampleId in "${sampleIds[@]}"; do
    echo "Processing sample: ${sample_id}"
    
    # Normalize the VCF file
    bcftools norm --multiallelics - ${inputdir}/${sampleId}.vcf.gz -o ${outputdir}/${sampleId}.norm.vcf
    
    # Convert the normalized VCF to binary ped format (BED)
    plink2 --vcf ${outputdir}/${sampleId}.norm.vcf --threads $PBS_NCPUS --keep-allele-order --make-bed --out "${outputdir}/${sampleId}"
    
    # Convert the BED format to 23andMe format
    /g/data/te53/kh3349/apps/plink --bfile "${outputdir}/${sampleId}" --threads $PBS_NCPUS --snps-only --recode 23 --out "${outputdir}/${sampleId}_23andme"
done

echo "conversion to 23andMe completed!"
