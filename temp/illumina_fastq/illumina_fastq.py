import csv
import gzip
from collections import defaultdict

def parse_fastq(fastq_file):
    """Parse the FASTQ file and return a list of sequences and quality scores."""
    sequences = []
    quality_scores = []
    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            # Sequence lines are the second line in every set of 4 lines
            if i % 4 == 1:
                sequences.append(line.strip())
            # Quality scores are the fourth line in every set of 4 lines
            elif i % 4 == 3:
                quality_scores.append(line.strip())
    return sequences, quality_scores

def count_reads(sequences):
    """Count the number of reads from the sequences."""
    return len(sequences)

def count_total_bases(sequences):
    """Count the total number of bases from the sequences."""
    return sum(len(seq) for seq in sequences)

def average_read_length(num_reads, total_bases):
    """Calculate the average read length."""
    return total_bases / num_reads if num_reads > 0 else 0

def phred_to_quality(phred_string):
    """Convert Phred encoded quality string to a list of rounded quality values."""
    return [ord(char) - 33 for char in phred_string]

def average_per_base_quality(quality_scores):
    """Calculate average per base quality values across all reads."""
    total_scores = [0] * max(len(qs) for qs in quality_scores)
    total_reads = len(quality_scores)
    
    for qs in quality_scores:
        scores = phred_to_quality(qs)
        for i, score in enumerate(scores):
            total_scores[i] += score

    return [round(score / total_reads) for score in total_scores]

def average_nucleotide_frequencies(sequences):
    """Calculate average per base nucleotide frequencies (A, C, G, T) across all reads."""
    bases = ['A', 'C', 'G', 'T']
    total_bases = [defaultdict(int) for _ in range(max(len(seq) for seq in sequences))]
    total_reads = len(sequences)
    
    for seq in sequences:
        for i, base in enumerate(seq):
            total_bases[i][base] += 1

    average_frequencies = []
    for base_counts in total_bases:
        avg_freq = [base_counts.get(base, 0) / total_reads for base in bases]
        average_frequencies.append(avg_freq)

    return average_frequencies

def overall_nucleotide_content(sequences):
    """Calculate overall nucleotide content (A, C, G, T) for the entire FASTQ file."""
    bases = ['A', 'C', 'G', 'T']
    total_counts = defaultdict(int)
    total_length = sum(len(seq) for seq in sequences)
    
    for seq in sequences:
        for base in seq:
            total_counts[base] += 1

    return [total_counts.get(base, 0) / total_length for base in bases]

def main():
    fastq_file = "/g/data/xl04/ka6418/ausarg/temp/illumina_fastq/alternate_test.fastq.gz"
    output_csv = "statistics.csv"
    
    sequences, quality_scores = parse_fastq(fastq_file)
    num_reads = count_reads(sequences)
    total_bases = count_total_bases(sequences)
    avg_read_length_val = average_read_length(num_reads, total_bases)
    avg_quality_values = average_per_base_quality(quality_scores)
    nucleotide_freq = average_nucleotide_frequencies(sequences)
    overall_content = overall_nucleotide_content(sequences)

    # Format nucleotide frequencies for CSV
    nucleotide_freq_str = [",".join(f"{freq:.2f}" for freq in freq_list) for freq_list in nucleotide_freq]
    overall_content_str = ",".join(f"{freq:.2f}" for freq in overall_content)

    # Save to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        headers = ['Number of reads', 'Total number of bases', 'Average read length', 
                   'Average per-base quality values', 'Average per-base nucleotide frequencies', 
                   'Overall nucleotide content (A,C,G,T)']
        csvwriter.writerow(headers)
        csvwriter.writerow([num_reads, total_bases, avg_read_length_val, 
                            ",".join(map(str, avg_quality_values)), 
                            "|".join(nucleotide_freq_str), overall_content_str])
    
    print(f"Statistics saved to {output_csv}")

if __name__ == "__main__":
    main()
