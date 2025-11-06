#!/bin/bash
#PBS -N indexasm
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=4:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd


module load samtools bwa-mem2

bwa-mem2 index ${fasta}
samtools faidx ${fasta}
