#!/bin/bash
#PBS -N teliprocator
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=5:00:00
#PBS -l mem=8GB
#PBS -l ncpus=1
#PBS -l storage=gdata/xl04+gdata/if89+gdata/te53
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

module load biopython

cd ${outdir}
prefix=$(basename "${fasta}" .fasta)

python3 /g/data/xl04/ka6418/temp/chromsyn_testrun/temp/telociraptor/code/telociraptor.py -seqin ${fasta} -basefile ${prefix} test i=-1 tweak=F telonull=T > log.txt

