#!/bin/bash
#PBS -N cutadapt
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd

set -ex 

module load cutadapt/3.7

cutadapt --cores ${PBS_NCPUS} --anywhere file:${PACBIOADAPTERS} \
 --error-rate 0.1 --overlap 25 --match-read-wildcards --revcomp --discard-trimmed \
 --json ${sample}.cutadapt.json \
 -o "${sample}.fastq" \
 ${fastq}