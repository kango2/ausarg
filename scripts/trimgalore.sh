#!/bin/bash
#PBS -N Illumina_Trimming_QC
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32GB
#PBS -l ncpus=8
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=100GB
#PBS -M kirat.alreja@anu.edu.au

source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate trim_env
module load fastqc/0.11.7

trim_galore --fastqc_args "--ourdir ${outdir_qc}" --outdir ${outdir_trim} --cores ${PBS_NCPUS} ${fastq}

