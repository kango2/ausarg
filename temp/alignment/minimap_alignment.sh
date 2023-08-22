#!/bin/bash
#PBS -N minimap alignment
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=5:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=400GB
#PBS -M kirat.alreja@anu.edu.au

# Run minimap2 and samtools
minimap2 -ax $minimap_preset -t $PBS_NCPUS $reference_assembly $query_fastq | \
samtools view -u - | \
samtools sort -u -m 4G -T $PBS_JOBFS/tmp -O BAM --reference $reference_assembly --threads 48 -o $output_filepath

# Index the output BAM file
samtools index -c $output_filepath