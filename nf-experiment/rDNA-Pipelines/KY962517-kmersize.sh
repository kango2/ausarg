#!/bin/bash
#PBS -P te53
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o KY962517-kmer.log

module load jellyfish

# Define the reference file path
REFERENCE="/g/data/te53/rdna/reference/KY962517.1.fasta"

# Define k-mer sizes
kmer_sizes=(9 11 13 17)

# Loop through each k-mer size and run jellyfish
for k in "${kmer_sizes[@]}"; do
    # Count the k-mers
    jellyfish count -m $k -s 100M -t 10 -C -o "KY962517.1-${k}mers.jf" $REFERENCE

    # Dump k-mer counts to a human-readable file
    jellyfish dump -c "KY962517.1-${k}mers.jf" > "KY962517.1-${k}mers_counts.txt"

    # Generate a histogram of k-mer occurrences
    jellyfish histo "KY962517.1-${k}mers.jf" > "KY962517.1-${k}mers_histo.txt"
done

echo "Jellyfish processing complete!"

