import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Load the CSV data
file_path = '/g/data/xl04/ka6418/ausarg/temp/qv_frequency_concurrent/tiliqua_ont.csv'
data = pd.read_csv(file_path)

# Pivot the table to get a 2D array for z
pivot_data = data.pivot_table(index='QV', columns='Read_Length', values='Read_Numbers', fill_value=0)

# Define the cutoff for QV
QV_CUTOFF = (80, 100)  # Adjust these values as needed

# Filter the data based on the cutoff
filtered_data = pivot_data[(pivot_data.index >= QV_CUTOFF[0]) & (pivot_data.index <= QV_CUTOFF[1])]

# Extract the x, y, and z values for contour plot after filtering
X_filtered, Y_filtered = np.meshgrid(filtered_data.columns, filtered_data.index)
Z_filtered = filtered_data.values

# Set up the figure and axis for the filtered range
fig, ax = plt.subplots(figsize=(12, 8))

# Create the contour fill
contourf = ax.contourf(X_filtered, Y_filtered, Z_filtered, levels=50, cmap='plasma')

# Add colorbar with label
cbar = plt.colorbar(contourf, ax=ax)
cbar.set_label('Read Numbers', fontsize=14)
cbar.ax.tick_params(labelsize=12)

# Set axis labels, title, and grid
ax.set_xlabel('Read Lengths', fontsize=16)
ax.set_ylabel('QV', fontsize=16)
title = f'Contour Plot of Read Numbers vs. Read Lengths and QV ({QV_CUTOFF[0]} ≤ QV ≤ {QV_CUTOFF[1]})'
ax.set_title(title, fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)

# Tight layout for better spacing
plt.tight_layout()

# Determine the output filename
output_filename = os.path.splitext(os.path.basename(file_path))[0] + '.png'

# Save the figure in high resolution
plt.savefig(output_filename, dpi=300)

# Close the plot
plt.close()
