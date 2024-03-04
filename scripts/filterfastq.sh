#!/bin/bash

#PBS -N filterfastq
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=16GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=100GB
#PBS -M kirat.alreja@anu.edu.au


module load pythonlib

python3 /g/data/xl04/ka6418/github/ausarg/scripts/python/filterfastq.py --input-seq ${input} --out ${outputdir} --filter-qual 8