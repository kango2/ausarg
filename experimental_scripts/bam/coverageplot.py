import numpy as np
import matplotlib.pyplot as plt

with open('/g/data/xl04/ka6418/sequence_alignment/bam_outputs/tiliqua_rugosa/pacbio/merged/assembly.haplotype1/coverage.txt', 'r') as f:
    data = np.loadtxt(f, dtype=int, usecols=2)

plt.hist(data, bins=50, range=(0, 100))
plt.title('Coverage Histogram')
plt.xlabel('Coverage Depth')
plt.ylabel('Frequency')
plt.savefig('histogram.png')
