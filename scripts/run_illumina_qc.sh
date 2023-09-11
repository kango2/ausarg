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

R1_file=$(echo ${filepair} | cut -f1 -d ':')
R2_file=$(echo ${filepair} | cut -f2 -d ':')
# Extract the base names without extensions
R1_basename=$(basename ${R1_file} | rev | cut -d. -f2- | rev)
R2_basename=$(basename ${R2_file} | rev | cut -d. -f2- | rev)

output_name=${output}/${R1_basename}_${R2_basename}_QC

module load python3
python3 /g/data/xl04/ka6418/ausarg/scripts/illumina_fastq_stats_parallel.py -R1 $R1_file -R2 $R2_file -o $output_name --format csv --cores ${PBS_NCPUS}

tail -n +2 ${output_name}.csv | sqlite3 /g/data/xl04/ka6418/ausarg/database/ausarg.db ".mode csv" ".import /dev/stdin illumina_metrics"
