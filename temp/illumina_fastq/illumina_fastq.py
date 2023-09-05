import csv
import gzip

def parse_fastq(fastq_file):
    """Parse the FASTQ file and return a list of sequences."""
    sequences = []
    with gzip.open(fastq_file, 'r') as f:
        for i, line in enumerate(f):
            # Sequence lines are the second line in every set of 4 lines
            if i % 4 == 1:
                sequences.append(line.strip())
    return sequences

def count_reads(sequences):
    """Count the number of reads from the sequences."""
    return len(sequences)

def count_total_bases(sequences):
    """Count the total number of bases from the sequences."""
    return sum(len(seq) for seq in sequences)

def average_read_length(num_reads, total_bases):
    """Calculate the average read length."""
    return total_bases / num_reads if num_reads > 0 else 0

def save_to_csv(output_file, num_reads, total_bases, avg_read_length):
    """Save the results to a CSV file."""
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Number of reads', 'Total number of bases', 'Average read length'])
        csvwriter.writerow([num_reads, total_bases, avg_read_length])

def main():
    fastq_file = "/g/data/xl04/ka6418/ausarg/temp/illumina_fastq/alternate_test.fastq.gz"
    output_csv = "statistics.csv"
    
    sequences = parse_fastq(fastq_file)
    num_reads = count_reads(sequences)
    total_bases = count_total_bases(sequences)
    avg_read_length_val = average_read_length(num_reads, total_bases)
    
    save_to_csv(output_csv, num_reads, total_bases, avg_read_length_val)
    print(f"Statistics saved to {output_csv}")

if __name__ == "__main__":
    main()
