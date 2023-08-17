#!/bin/bash

#PBS -N master_qv_concurrent
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=16GB
#PBS -l ncpus=4
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=400GB
#PBS -M kirat.alreja@anu.edu.au


input_file="$1"  # Path to the input gzipped FASTQ file
csv_dir="$2"

mkdir -p "$csv_dir/chunks"
mkdir -p "$csv_dir/csvs"

# Split the input file into chunks within the temporary directory
zcat "$input_file" | split -l 40000 -d - "${PBS_JOBFS}/input_chunk_"

./submission.sh "$csv_dir/chunks" "$csv_dir/csvs"





