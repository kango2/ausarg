#!/bin/bash
#PBS -N kat-experimental
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=400GB
#PBS -M kirat.alreja@anu.edu.au

module load kat
module load pythonlib

# Change this to the path of your fasta file
FASTA_FILE="/g/data/xl04/ka6418/species/Tiliqua_Rugosa/rTilRug0.5/assembly/rTilRug0.5.asm.hp1.fasta"

# Set the output directory
OUTPUT_DIR="/g/data/xl04/ka6418/ausarg-data"

# K-mer size to use
KMER_SIZE=75

# Run KAT comp with the fasta file as input, comparing the file to itself
kat comp -t $PBS_NCPUS -o "${OUTPUT_DIR}/kat_comp" -m $KMER_SIZE -H 1000000 "${FASTA_FILE}" "${FASTA_FILE}"

# Extract the main diagonal value, representing homozygous k-mers
MAIN_DIAGONAL=$(grep main_diagonal "${OUTPUT_DIR}/kat_comp-main.mx" | awk '{print $2}')

# Extract the total number of k-mers
TOTAL_KMERS=$(grep total "${OUTPUT_DIR}/kat_comp-main.mx" | awk '{print $2}')

# Calculate the heterozygosity rate
heterozygosity_rate=$(echo "scale=6; 1 - (${MAIN_DIAGONAL} / ${TOTAL_KMERS})" | bc)

# Print the heterozygosity rate
echo "Heterozygosity rate: ${heterozygosity_rate}"
