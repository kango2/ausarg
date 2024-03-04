#!/bin/bash
#PBS -N greenhill
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=400GB
#PBS -M kirat.alreja@anu.edu.au


source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate greenhill

cd ${outputdir}

greenhill \
-c ${fasta} \
-HIC ${r1} ${r2} \
-t ${PBS_NCPUS} \
2>3D.log