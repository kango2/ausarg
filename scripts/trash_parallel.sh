#!/bin/bash

#PBS -N TRASH_parallel
#PBS -P xl04
#PBS -q express
#PBS -l walltime=00:30:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=100GB
#PBS -M kirat.alreja@anu.edu.au

module load parallel/20191022 kentutils/0.0 

inputfile="$inputfasta"
output="$outputdir"
seqt="$template"

source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate trash

faSplit sequence "$inputfile" 10000 "${PBS_JOBFS}/chunk"

cd "${PBS_JOBFS}"
filelist=$(ls ${PBS_JOBFS}/chunk*)

num_jobs=${PBS_NCPUS}

if [ -z "$seqt" ]; then
    printf "%s\n" "${filelist[@]}" | parallel  -I{} --jobs "$num_jobs" /g/data/xl04/ka6418/TRASH/TRASH_run.sh  {} --o "${PBS_JOBFS}" --N.max.div 5
else
    printf "%s\n" "${filelist[@]}" | parallel  -I{} --jobs "$num_jobs" /g/data/xl04/ka6418/TRASH/TRASH_run.sh  {} --o "${PBS_JOBFS}" --N.max.div 5 --seqt "$seqt"
fi

awk 'NR==1{print; next} FNR>1' Summary*.csv > "${output}/$(basename "$inputfile")_centromeres.csv"
