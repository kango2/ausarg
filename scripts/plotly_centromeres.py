import pandas as pd
import plotly.subplots as sp
import plotly.graph_objects as go
import math

# Read the CSV file into a Pandas DataFrame
data = pd.read_csv('/g/data/xl04/ka6418/greenhill/bassiana_greenhill/centromeres/Summary.of.repetitive.regions.out_ConsensusOutput.fa.csv')

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
fig.write_html('chromosome_plot_bassiana_kirat.html')
