import matplotlib.pyplot as plt

# List of sample names
samples = ['ELD100106','ELD102637','ELD102679','ELD110586']

# Data structures to hold all sample data
all_depths = []
all_scores_mismatches = []
all_positions = []

# Load all the data first
for sample in samples:
    # Paths for the current sample
    depth_file = f"/g/data/te53/kh3349/rDNA/multi-analysis/depth_files/{sample}_depth.txt"
    alignment_file = f"/g/data/te53/kh3349/rDNA/multi-analysis/alignment_stats/{sample}_alignment_info.txt"
    start_positions_file = f"/g/data/te53/kh3349/rDNA/multi-analysis/alignment_stats/{sample}_start_positions.txt"

    # --- Coverage Visualization ---
    positions, depths = [], []
    with open(depth_file, 'r') as f:
        for line in f:
            pos, depth = line.strip().split()[1:3]
            positions.append(int(pos))
            depths.append(int(depth))
    all_depths.append(depths)

    # --- Alignment Quality Visualization ---
    scores, mismatches = [], []
    with open(alignment_file, 'r') as f:
        for line in f:
            score, mismatch, _ = line.strip().split()
            scores.append(int(score.split(":")[-1]))
            mismatches.append(int(mismatch.split(":")[-1]))
    all_scores_mismatches.append((scores, mismatches))

    # --- Start Positions Visualization ---
    pos = []
    with open(start_positions_file, 'r') as f:
        for line in f:
            pos.append(int(line.strip()))
    all_positions.append(pos)

# Now, plot combined figures
fig, axs = plt.subplots(len(samples), 3, figsize=(20, 6 * len(samples)))

# Ensure axs is always a 2D array, even if len(samples) is 1
if len(samples) == 1:
    axs = [axs]

for idx, sample in enumerate(samples):
    axs[idx][0].plot(positions, all_depths[idx], color="blue", label="Coverage")
    axs[idx][0].set_title(f"Coverage for {sample}")
    axs[idx][0].set_xlabel("Position in reference")
    axs[idx][0].set_ylabel("Coverage depth")
    axs[idx][0].legend()

    axs[idx][1].scatter(all_scores_mismatches[idx][0], all_scores_mismatches[idx][1], color="red", alpha=0.5)
    axs[idx][1].set_title(f"Alignment Quality for {sample}")
    axs[idx][1].set_xlabel("Alignment Score")
    axs[idx][1].set_ylabel("Number of Mismatches")

    axs[idx][2].hist(all_positions[idx], bins=100, color="green", alpha=0.7)
    axs[idx][2].set_title(f"Start Positions for {sample}")
    axs[idx][2].set_xlabel("Position in reference")
    axs[idx][2].set_ylabel("Number of Reads")

plt.tight_layout()
plt.savefig("/g/data/te53/kh3349/rDNA/multi-analysis/combined_plots.png")
plt.show()
