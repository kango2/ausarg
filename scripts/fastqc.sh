#!/bin/bash
#PBS -N Illumina_FastQC
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32GB
#PBS -l ncpus=8
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=100GB
#PBS -M kirat.alreja@anu.edu.au

module load fastqc/0.11.7
fastqc --threads 8 --outdir ${outdir} ${fastq}

