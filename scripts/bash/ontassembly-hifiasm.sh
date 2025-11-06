#!/bin/bash
#PBS -N hifiasm
#PBS -P xl04
#PBS -q hugemem
#PBS -l walltime=48:00:00
#PBS -l mem=1470GB
#PBS -l ncpus=48
#PBS -l jobfs=1400GB
#PBS -l storage=gdata/if89+gdata/xl04
#PBS -l wd
#PBS -j oe 


#Usage
#outdir : output folder
#outbase : basename for hifiasm output files 
#inputfastq : ":" colon separated fastq files

set -ex 

module load hifiasm/0.24.0

cat ${inputfastq//:/ } > ${PBS_JOBFS}/${outbase}.fastq.gz

hifiasm -t ${PBS_NCPUS} --ont -o ${outdir}/${outbase} ${PBS_JOBFS}/${outbase}.fastq.gz