#!/bin/bash

# Variables for r1, r2, and outputdir
r1="your_r1_variable"  # Replace with your actual r1 variable
r2="your_r2_variable"  # Replace with your actual r2 variable
outputdir="your_outputdir"  # Replace with your actual outputdir

# Directory containing the fasta files
fasta_dir="/g/data/xl04/genomeprojects/Pogona_vitticeps/fasta"

# Output CSV file
output_csv="output.csv"

# Create the CSV header
echo "sample,hap,fasta,r1,r2,outputdir" > "$output_csv"

# Loop through each fasta file in the directory
for fasta_file in "$fasta_dir"/*.fasta; do
    # Extract the filename without the path
    filename=$(basename "$fasta_file")
    
    # Extract the sample (POGVITdef or POGVITscaf)
    sample="${filename%%.*}"

    # Extract the hap (h1, h2, or p)
    hap="${filename#*.}"
    hap="${hap%%.*}"

    # Write the row to the CSV file
    echo "${sample},${hap},${fasta_file},${r1},${r2},${outputdir}" >> "$output_csv"
done

echo "CSV file generated: $output_csv"
