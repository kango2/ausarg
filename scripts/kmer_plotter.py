# todo make it user-friendly and add arguments

import os
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# Load data from file
def load_data(file_path):
    """Load data from the given file path into lists of x and y values."""
    x_values, y_values = [], []
    with open(file_path, 'r') as f:
        for line in f:
            x, y = line.strip().split()
            x_values.append(int(x))
            y_values.append(int(y))
    return x_values, y_values

# Extract info from filename
def extract_info_from_filename(filename):
    """Extract specie name, sequencing tech, and kmer size from the filename."""
    base_name = os.path.basename(filename).replace('.histo', '')
    specie_name, sequencing_tech, kmer_size = base_name.split("_")[:3]
    kmer_size = int(kmer_size)  # Convert kmer size to integer
    return specie_name, sequencing_tech, kmer_size

# Filter data based on x-axis limit
def filter_data(x_data, y_data, limit=500):
    filtered_x = [x for x in x_data if x <= limit]
    filtered_y = y_data[:len(filtered_x)]
    return filtered_x, filtered_y
def automatic_plot_combined(files,save_path=None):
    file_info = [extract_info_from_filename(file) for file in files]
    files = [os.path.join(directory_path, filename) for filename in os.listdir(directory_path) if filename.endswith('.histo')]

    # Group files by sequencing tech
    grouped_files = defaultdict(list)
    for file, (specie_name, tech, kmer_size) in zip(files, file_info):
        grouped_files[tech].append((file, specie_name, kmer_size))

    total_rows = len(grouped_files)
    max_subplots_per_row = max(len(tech_files) for tech_files in grouped_files.values())

    # Create a single figure with multiple subplots arranged in rows by sequencing tech
    fig, ax = plt.subplots(total_rows, max_subplots_per_row, figsize=(5 * max_subplots_per_row, 5 * total_rows))

    # To handle cases where there's only one sequencing tech
    if total_rows == 1:
        ax = [ax]

    row_idx = 0
    for tech, tech_files in grouped_files.items():
        tech_files = sorted(tech_files, key=lambda x: x[2])  # Sort by k-mer size

        for col_idx, (file, specie_name, kmer_size) in enumerate(tech_files):
            x_data, y_data = filter_data(*load_data(file), limit=200)  # Set default limit to 200

            ax[row_idx][col_idx].bar(x_data, y_data, width=1, align='center', color=sns.color_palette("flare")[col_idx])
            ax[row_idx][col_idx].set_title(f"{specie_name} ({tech}, k={kmer_size})")
            ax[row_idx][col_idx].set_xlabel("Kmer Coverage")
            ax[row_idx][col_idx].set_yscale('log')

        row_idx += 1

    # Common y-axis label
    fig.text(0.00, 0.5, 'Frequency', va='center', rotation='vertical')
    
    # Set the mega title
    mega_title = "QV = 10"
    plt.suptitle(mega_title, fontsize=16)

    # Adjust layout and show plot
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust top spacing for mega title
    
    fig.set_facecolor("white")
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
    plt.show()


directory_path = "/g/data/xl04/bpadata/Pogona_vitticeps/raw/evaluation/kmer"

save_path = "/g/data/xl04/ka6418/gc/PogVit_Kmers_10.png"

# For demonstration purposes using the previously loaded files
automatic_plot_combined(files,save_path)
