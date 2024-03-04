#!/bin/bash
#PBS -N asm_vs_asm
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=4:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48 
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

module load minimap2

# Assuming outputdir, ref, and target are passed as environment variables or defined elsewhere in the script
output_file="${outputdir}/$(basename ${ref})_$(basename ${target}).paf"

minimap2 -t 48 -o ${output_file} -x asm20 ${ref} ${target}

#sort -k1,1 -k4,4n
