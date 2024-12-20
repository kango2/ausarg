#!/bin/bash
#PBS -N olfactory
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l jobfs=400GB
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate orfinder
module load hmmer
cp ${inputgenome} ${PBS_JOBFS}
perl /g/data/xl04/ka6418/annotation/koala/ORFinder.pl -reffasta ${PBS_JOBFS}/$(basename ${inputgenome}) -speciescode $(basename ${inputgenome} .fasta) -outputdir ${outputdir} -threads ${PBS_NCPUS}