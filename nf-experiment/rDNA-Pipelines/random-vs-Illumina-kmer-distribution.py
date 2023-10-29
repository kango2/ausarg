import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Create a directory for saving outputs
output_directory = "/g/data/te53/kh3349/rDNA/random-vs-rDNA/visualizations"
os.makedirs(output_directory, exist_ok=True)

# Paths to directories with k-mer results
randomregions_paths = {
    9: "/g/data/te53/nj8315/randomregions/random-illumina-k9",
    11: "/g/data/te53/nj8315/randomregions/random-illumina-k11",
    13: "/g/data/te53/nj8315/randomregions/random-illumina-k13",
    17: "/g/data/te53/nj8315/randomregions/random-illumina-k17"
}

rdna_paths = {
    9: "/g/data/te53/nj8315/updated-illumina/Illumina-k9",
    11: "/g/data/te53/nj8315/updated-illumina/Illumina-k11",
    13: "/g/data/te53/nj8315/updated-illumina/Illumina-k13",
    17: "/g/data/te53/nj8315/updated-illumina/Illumina-k17"
}

# Modified script to add a grid overlay on both x and y axes

def load_kmer_counts(directory, suffix):
    """Load k-mer counts from Jellyfish .histo files in a directory."""
    files = glob.glob(f"{directory}/*{suffix}")
    aggregated_counts = {}
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                count, num_kmers = map(int, line.strip().split())
                aggregated_counts[count] = aggregated_counts.get(count, 0) + num_kmers
    return aggregated_counts

# Visualization section with the modified axes and grid overlay on both x and y axes
for k in [9, 11, 13, 17]:
    # Load data
    random_counts = load_kmer_counts(randomregions_paths[k], ".histo")
    rdna_counts = load_kmer_counts(rdna_paths[k], ".fastq.gz.histo")  # Adjusted suffix for rdna_paths
    
    # Create a new figure
    plt.figure(figsize=(20, 8))
    
    # Bar Plot for Random Regions
    plt.bar(random_counts.keys(), random_counts.values(), width=0.4, align='center', label='Random Regions', alpha=0.7, color='blue')
    
    # Bar Plot for rDNA Regions (shifted slightly to avoid overlap)
    plt.bar([k+0.5 for k in rdna_counts.keys()], rdna_counts.values(), width=0.4, align='center', label='rDNA Regions', alpha=0.7, color='red')
    
    plt.yscale('log')
    plt.title(f'K-mer Count Distribution for k={k}')
    plt.xlabel('Occurrences of K-mer')
    plt.ylabel('Number of Unique K-mers (log scale)')
    plt.legend()
    plt.grid(axis='both', which='both', linestyle='--', linewidth=0.5, alpha=0.5)  # Added grid overlay to both axes
    
    # Save the figure
    plt.tight_layout()
    # Uncomment the below line to save the figure when you run the script in your environment
    plt.savefig(f"{output_directory}/corrected_combined_k{k}.png")
    plt.show()
