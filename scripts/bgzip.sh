#!/bin/bash
#PBS -N blow5merge
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au


module load htslib

bgzip -@ ${PBS_NCPUS} ${file}
