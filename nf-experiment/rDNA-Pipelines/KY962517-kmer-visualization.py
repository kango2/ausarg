import pandas as pd
import matplotlib.pyplot as plt

# Updated file paths with the provided directory path
file_paths = {
    9: "/g/data/te53/kh3349/rDNA/KY962517-kmer/KY962517.1-9mers_histo.txt",
    11: "/g/data/te53/kh3349/rDNA/KY962517-kmer/KY962517.1-11mers_histo.txt",
    13: "/g/data/te53/kh3349/rDNA/KY962517-kmer/KY962517.1-13mers_histo.txt",
    17: "/g/data/te53/kh3349/rDNA/KY962517-kmer/KY962517.1-17mers_histo.txt"
}

plt.figure(figsize=(15, 10))

for k, file_path in file_paths.items():
    data = pd.read_csv(file_path, sep=" ", header=None, names=["K-mer Count", "Frequency"])
    plt.loglog(data["K-mer Count"], data["Frequency"], label=f"k={k}")

plt.title("K-mer Histograms for Different Sizes")
plt.xlabel("K-mer Count")
plt.ylabel("Frequency")
plt.legend()
plt.grid(True, which="both", ls="--", c='0.65')
plt.tight_layout()

# Save the plot as a PDF
output_pdf_path = "/g/data/te53/kh3349/rDNA/KY962517-kmer/kmer_histograms.pdf"
plt.savefig(output_pdf_path, format='pdf')

print(f"Saved k-mer histograms as {output_pdf_path}")

