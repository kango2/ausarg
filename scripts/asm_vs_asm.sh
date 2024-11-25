#!/bin/bash
#PBS -N asm_vs_asm
#PBS -P xl04
#PBS -q express
#PBS -l walltime=2:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -j oe

#usage : qsub -o /path/to/logs -v outputdir=/path/to/output,target=/path/to/target.fa,query=path/to/query.fa path/to/asm_vs_asm.sh

module load minimap2

# Assuming outputdir, ref, and target are passed as environment variables or defined elsewhere in the script
output_file="${outputdir}/$(basename ${target})_$(basename ${query}).paf"

minimap2 -t ${PBS_NCPUS} -o ${output_file} -x asm20 ${target} ${query}  
