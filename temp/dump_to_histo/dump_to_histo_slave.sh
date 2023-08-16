#!/bin/bash
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=0:30:00
#PBS -l mem=32GB
#PBS -l ncpus=16
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

input_file="$INPUT_FILE"

if [ ! -f "$input_file" ]; then
    echo "Input file not found."
    exit 1
fi

directory=$(dirname "$input_file")
filename=$(basename "$input_file")
base_name="${filename%.*}"
output_file="$directory/$base_name.histo"

jellyfish histo -t 16 "$input_file" > "$output_file"

echo "Histogram generated and saved as: $output_file"
