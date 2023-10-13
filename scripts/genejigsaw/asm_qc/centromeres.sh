#!/bin/bash

#PBS -N TRASH_parallel
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -j oe
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=100GB
#PBS -M kirat.alreja@anu.edu.au

inputfile="$inputfasta"
output="$outputdir"
seqt="$template"

source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate trash

if [ -z "$seqt" ]; then
    echo "No sequence template"
    /g/data/xl04/ka6418/TRASH/TRASH_run.sh ${inputfile} --par ${PBS_NCPUS} --o ${output} --N.max.div 5
else
    echo "Running with template"
    /g/data/xl04/ka6418/TRASH/TRASH_run.sh ${inputfile} --par ${PBS_NCPUS} --o ${output} --N.max.div 5 --seqt "$seqt"
fi


