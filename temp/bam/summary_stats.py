import csv

def calculate_summary_stats(depth_values):
    """Calculate summary statistics for a list of depth values."""
    stats = {
        'mean': sum(depth_values) / len(depth_values),
        'median': sorted(depth_values)[len(depth_values) // 2],
        'min': min(depth_values),
        'max': max(depth_values),
    }
    return stats

def process_depth_file(file_path, output_csv):
    haplotype_stats = {}

    with open(file_path, 'r') as file:
        for line in file:
            # Assuming the depth file has three columns: chromosome, position, and depth
            chromosome, position, depth = line.strip().split('\t')

            # Extract the haplotype information from the chromosome name (modify as per your naming convention)
            haplotype = chromosome.split(':')[0]

            # If haplotype not encountered before, initialize an empty list for depth values
            if haplotype not in haplotype_stats:
                haplotype_stats[haplotype] = []

            # Add depth value to the haplotype's list
            haplotype_stats[haplotype].append(int(depth))

    # Calculate summary statistics for each haplotype
    haplotype_summary_stats = {
        haplotype: calculate_summary_stats(depth_values)
        for haplotype, depth_values in haplotype_stats.items()
    }

    # Get the fieldnames from the first haplotype's stats
    fieldnames = ['Haplotype'] + list(next(iter(haplotype_summary_stats.values())).keys())

    # Write the summary statistics to a CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for haplotype, stats in haplotype_summary_stats.items():
            writer.writerow({'Haplotype': haplotype, **stats})


# Usage example
depth_file = '/g/data/xl04/ka6418/sequence_alignment/bam_outputs/tiliqua_rugosa/pacbio/merged/assembly.haplotype1/coverage.txt'  # Replace with the path to your depth file
output_csv = 'summary_stats.csv'  # Replace with the desired path for the output CSV file
process_depth_file(depth_file, output_csv)
