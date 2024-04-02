#!/bin/bash
#PBS -N alignment
#PBS -P xl04
#PBS -q normal
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l walltime=48:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l wd

module load minimap2 samtools

platform=$platform
rawreads=$rawreads
reference=$reference
output=$output

