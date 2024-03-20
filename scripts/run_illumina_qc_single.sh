#!/bin/bash
#PBS -N IlluminaQC
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=96GB
#PBS -l ncpus=24
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

# Extract the base names without extensions
R1_basename=$(basename ${R1} | rev | cut -d. -f2- | rev)

output_name=${output}/${R1_basename}_QC

module load python3
python3 /g/data/xl04/ka6418//github/ausarg/scripts/illumina_fastq_stats_parallel.py -R1 $R1 -o $output_name --format csv --cores ${PBS_NCPUS}

#TODO: add update to status, make it 3 keyword system, implement "running"
