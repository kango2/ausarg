import csv

def count_reads(fastq_file):
    """Count the number of reads in a FASTQ file."""
    with open(fastq_file, 'r') as f:
        return sum(1 for line in f if line.startswith('@')) // 4

def count_total_bases(fastq_file):
    """Count the total number of bases in a FASTQ file."""
    total_bases = 0
    with open(fastq_file, 'r') as f:
        for i, line in enumerate(f):
            # Sequence lines are the second line in every set of 4 lines
            if i % 4 == 1:
                total_bases += len(line.strip())
    return total_bases

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
    fastq_file = "path_to_your_fastq_file.fastq"
    output_csv = "statistics.csv"

    num_reads = count_reads(fastq_file)
    total_bases = count_total_bases(fastq_file)
    avg_read_length = average_read_length(num_reads, total_bases)
    
    save_to_csv(output_csv, num_reads, total_bases, avg_read_length)
    print(f"Statistics saved to {output_csv}")

if __name__ == "__main__":
    main()
