import csv
from collections import defaultdict
import heapq

class StreamingMedian:
    def __init__(self):
        self.min_heap = []  # to store the higher half of the list
        self.max_heap = []  # to store the lower half of the list

    def add_number(self, num):
        if not self.max_heap or num < -self.max_heap[0]:
            heapq.heappush(self.max_heap, -num)
        else:
            heapq.heappush(self.min_heap, num)

        # Balance the heaps
        if len(self.max_heap) < len(self.min_heap):
            heapq.heappush(self.max_heap, -heapq.heappop(self.min_heap))
        elif len(self.max_heap) > len(self.min_heap) + 1:
            heapq.heappush(self.min_heap, -heapq.heappop(self.max_heap))

    def find_median(self):
        if len(self.max_heap) == len(self.min_heap):
            return (-self.max_heap[0] + self.min_heap[0]) / 2
        else:
            return -self.max_heap[0]

# Initialize defaultdict for mean, min, max, and median
means = defaultdict(int)
counts = defaultdict(int)
mins = defaultdict(int)
maxs = defaultdict(int)
medians = defaultdict(StreamingMedian)

# Open the csv file and read it line by line
with open('/g/data/xl04/ka6418/sequence_alignment/bam_outputs/tiliqua_rugosa/pacbio/merged/assembly.haplotype1/coverage.txt', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for chrom, _, depth in reader:
        depth = int(depth)
        counts[chrom] += 1
        means[chrom] += (depth - means[chrom]) / counts[chrom]  # Online mean calculation
        mins[chrom] = min(mins[chrom], depth) if chrom in mins else depth
        maxs[chrom] = max(maxs[chrom], depth) if chrom in maxs else depth
        medians[chrom].add_number(depth)

# Open a new csv file to write the results
with open('coverage_stats.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Haplotype", "Mean", "Median", "Min", "Max"])  # Write the header

    # For each haplotype, write the results to the csv file
    for chrom in means.keys():
        writer.writerow([chrom, means[chrom], medians[chrom].find_median(), mins[chrom], maxs[chrom]])
