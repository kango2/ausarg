#!/bin/bash

#PBS -N master_qv_concurrent
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l mem=32GB
#PBS -l ncpus=16
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=400GB
#PBS -M kirat.alreja@anu.edu.au


input_file="$INPUT_FILE"  # Path to the input gzipped FASTQ file
temp_dir="$TEMP_DIR"

mkdir -p "$temp_dir/chunks"
mkdir -p "$temp_dir/csvs"

input_base=$(basename "$input_file" .fq.gz)  # Remove .fq.gz extension
output_prefix="${temp_dir}/chunks/${input_base}_chunk_"

zcat "$input_file" | split -a 3 -l 40000 -d - "$output_prefix"
for file in "${output_prefix}"*; do
    mv "$file" "$file.fq"
done


./submission.sh "$temp_dir/chunks" "$temp_dir/csvs"

