#!/bin/bash
#PBS -N qc
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=96GB
#PBS -l ncpus=1
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

module load biopython

python3 /g/data/xl04/ka6418/github/ausarg/scripts/long-read-qv.py --i ${input} --o ${output} --f ${flowcell} --p ${platform} --s ${sample}

