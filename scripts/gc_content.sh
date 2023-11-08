#!/bin/bash
#PBS -N gc_calculation
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=05:00:00
#PBS -l mem=32GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=100GB
#PBS -M kirat.alreja@anu.edu.au

#TODO : Add usage and flags for public use

module load parallel/20191022 biopython/1.79 kentutils/0.0

inputfile="$input"
outputdir="$output"

faSplit sequence $inputfile 10000 ${PBS_JOBFS}/chunk

cd ${PBS_JOBFS}
filelist=$(ls ${PBS_JOBFS}/chunk*)
printf "%s\n" "${filelist[@]}" | parallel -I{} --jobs ${PBS_NCPUS} python3 /g/data/xl04/ka6418/github/ausarg/scripts/gc_content_helper.py {} 

awk 'NR==1{print; next} FNR>1' *.csv > "${outputdir}/$(basename "$inputfile")_GC.csv"
