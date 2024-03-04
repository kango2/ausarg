import pandas as pd
import plotly.subplots as sp
import plotly.graph_objects as go
import math
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate plots for chromosome data.')
parser.add_argument('-i', '--input_file', required=True, help='Input file path')
parser.add_argument('-o', '--output_folder', required=True, help='Output folder path')
args = parser.parse_args()

# Get the base name of the input file and construct the output file name
input_file_base = os.path.basename(args.input_file)
output_file_name = os.path.splitext(input_file_base)[0] + '.html'
output_file = os.path.join(args.output_folder, output_file_name)

# Read the CSV file into a Pandas DataFrame
data = pd.read_csv(args.input_file)

# Get the unique chromosomes
chromosomes = data['name'].unique()

# Define layout
cols = 3  # Number of columns
rows = math.ceil(len(chromosomes) / cols)  # Calculate number of rows needed

# Create subplots
fig = sp.make_subplots(rows=rows, cols=cols, subplot_titles=chromosomes)

# Generate individual graphs per chromosome
for i, chromosome in enumerate(chromosomes):
    row = i // cols + 1
    col = i % cols + 1
    chromosome_data = data[data['name'] == chromosome]

    fig.add_trace(go.Scatter(x=chromosome_data['start'], y=chromosome_data['width'],
                             mode='markers', name=chromosome, hovertext=chromosome_data['most.freq.value.N']),
                  row=row, col=col)

    fig.update_xaxes(title_text='Start Coordinate', row=row, col=col)
    fig.update_yaxes(title_text='Width of Repeat', row=row, col=col)

fig.update_layout(height=200 * rows, width=1200, showlegend=False)  # Adjust height and width as needed

# Display the plot
fig.show()
fig.write_html(output_file)
