#!/bin/bash
# Hap.py-based Haplotype Comparison Script
# This script is designed for comparing haplotypes using Hap.py.
# It requires the Hap.py tool installed using "singularity pull docker://pkrusche/hap.py".
# Authors: Kosar Hooshmand

#PBS -l ncpus=48
#PBS -l walltime=48:00:00
#PBS -l mem=190GB
#PBS -N HapComparison
#PBS -j oe
#PBS -l jobfs=400GB
#PBS -o hap.log
#PBS -l storage=scratch/te53+gdata/te53+gdata/xy86+gdata/if89
#PBS -P te53

module load singularity

# singularity pull docker://pkrusche/hap.py

# Define the paths and parameters
SINGULARITY_IMAGE="/g/data/te53/kh3349/apps/hap.py_latest.sif"
HAPPY_PATH="/opt/hap.py/bin/hap.py"
TRUESET_VCF="/g/data/te53/kh3349/Genome_Reference_Materials/giab_data/NA12878_HG001/HG001_GRCh38_1_22_v4.2.1_benchmark_annotated.vcf.gz"
BED_FILE="/g/data/te53/kh3349/Genome_Reference_Materials/giab_data/NA12878_HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
OUTPUT_DIR="/g/data/te53/kh3349/varaints-comparison/hap.py_output"
REFERENCE="/g/data/te53/humanreference/GRCh38.p13/refresource/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.fa"

# SNP VCF
ANALYSED_SNP_VCF="/g/data/xy86/nextflow-outputs/HG001/variantannotator_SNP/SNP_recal.annotated.vcf.gz"
export HGREF=$REFERENCE
singularity exec $SINGULARITY_IMAGE $HAPPY_PATH $TRUESET_VCF $ANALYSED_SNP_VCF -o $OUTPUT_DIR/SNP -f $BED_FILE --bcftools-norm --force-interactive --window-size 30 --threads ${PBS_NCPUS}

# INDEL VCF
ANALYSED_INDEL_VCF="/g/data/xy86/nextflow-outputs/HG001/variantannotator_INDEL/INDEL_recal.annotated.vcf.gz"
export HGREF=$REFERENCE
singularity exec $SINGULARITY_IMAGE $HAPPY_PATH $TRUESET_VCF $ANALYSED_INDEL_VCF -o $OUTPUT_DIR/INDEL -f $BED_FILE --bcftools-norm --force-interactive --window-size 30 --threads ${PBS_NCPUS}
