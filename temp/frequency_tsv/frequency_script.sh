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

# Construct output directory path
output_dir="$(dirname "$input_file")/../evaluation/length_frequency/tsv"
mkdir -p "$output_dir"  # Create the directory if it doesn't exist

id=$(basename "$input_file" .fq.gz | rev | cut -d "_" -f 1 | rev)
sequencing_technology=$(basename "$input_file" .fq.gz | rev | cut -d "_" -f 2 | rev)

# Construct output file name
output_file_name=${id}_${sequencing_technology}_length_freq.tsv"

echo -e "Read_Length\tRead_Numbers\tFlowcell_ID\tPlatform" > "$output_file_name"

if [ ! -r "$input_file" ]; then
  echo "Input file '$input_file' not found or not readable."
  exit 1
fi

zcat "$input_file" | sed -n '2~4p' | perl -lne '$l{int(length($_)/100)*100}++; END { map { print "$_\t$l{$_}\t'"$id"'\t'"$sequencing_technology"'"} sort {$a<=>$b} keys %l }' >> "$output_file_name"

echo "Processing completed. Output file saved as: $output_file_name"
