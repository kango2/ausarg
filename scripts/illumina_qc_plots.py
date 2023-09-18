import os
def generate_qc_plot(csv_path, output_png_path):
    # Load the CSV file into a DataFrame
    qc_metrics_df = pd.read_csv(csv_path)
    
    # Splitting columns to extract data for the two paired-end files
    qc_metrics_df['Mean_QV_at_read_position_1'], qc_metrics_df['Mean_QV_at_read_position_2'] = qc_metrics_df['Mean_QV_at_read_position'].str.split(';', 1).str
    
    # Extracting Mean QV at read position for both paired-end files
    mean_qv_1 = list(map(float, qc_metrics_df['Mean_QV_at_read_position_1'][0].split(',')))
    mean_qv_2 = list(map(float, qc_metrics_df['Mean_QV_at_read_position_2'][0].split(',')))
    
    # Extract base names for titles
    file_path_content = qc_metrics_df['File_path'][0]
    base_names = [os.path.basename(path) for path in file_path_content.split(';')]
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8), sharey=True)
    
    # Plotting data for Paired-end 1
    ax1.plot(mean_qv_1, color='blue', alpha=0.8)
    ax1.set_title(base_names[0])
    ax1.set_xlabel('Read Position')
    ax1.set_ylabel('Mean QV')
    ax1.set_ylim(0, 50)
    
    # Plotting data for Paired-end 2
    ax2.plot(mean_qv_2, color='green', alpha=0.8)
    ax2.set_title(base_names[1])
    ax2.set_xlabel('Read Position')
    ax2.set_ylim(0, 50)
    
    # Adding main title to the plot
    fig.suptitle('Mean Quality Value (QV) at Read Position')
    
    # Save the plot to a .png file
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.savefig(output_png_path)
    plt.close()

def generate_nucleotide_frequency_plot(csv_path, output_png_path):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(csv_path)
    
    # Split the 'Nucleotide_count_at_read_position' column to get data for R1 and R2
    r1_nucleotide_data, r2_nucleotide_data = df['Nucleotide_count_at_read_position'].iloc[0].split(';')
    
    # Extract DataFrames for R1 and R2
    r1_df = parse_nucleotide_data(r1_nucleotide_data)
    r2_df = parse_nucleotide_data(r2_nucleotide_data)
    
    # Calculate frequencies for R1 and R2
    r1_freq_df = calculate_frequencies(r1_df)
    r2_freq_df = calculate_frequencies(r2_df)
    
    # Update the color for 'N' to pastel orange
    colors['N'] = '#FFD6A5'
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # R1 Stacked Bar Chart
    r1_freq_df.plot(kind='bar', stacked=True, color=[colors[nucleotide] for nucleotide in r1_freq_df.columns], ax=ax1, legend=False)
    ax1.set_title("R1 Nucleotide Frequencies")
    ax1.set_ylabel("Frequency")
    ax1.set_xlabel("Position")
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xticks, rotation=45)
    ax1.set_ylim(0, 1)

    # R2 Stacked Bar Chart
    r2_freq_df.plot(kind='bar', stacked=True, color=[colors[nucleotide] for nucleotide in r2_freq_df.columns], ax=ax2)
    ax2.set_title("R2 Nucleotide Frequencies")
    ax2.set_ylabel("Frequency")
    ax2.set_xlabel("Position")
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticks, rotation=45)
    ax2.set_ylim(0, 1)
    
    # Layout adjustment
    plt.tight_layout()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    # Save the plot to the provided output path
    plt.savefig(output_png_path)
    plt.close()

def generate_qv_histogram(csv_path, output_png_path):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(csv_path)
    
    # Split the 'Mean_QV_per_read' column to get data for R1 and R2
    r1_qv_data, r2_qv_data = df['Mean_QV_per_read'].iloc[0].split(';')
    
    # Get Series for R1 and R2 QV frequencies
    r1_qv_series = parse_qv_data(r1_qv_data)
    r2_qv_series = parse_qv_data(r2_qv_data)
    
    # Plotting histograms for QV frequencies of R1 and R2
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # R1 Histogram
    ax1.bar(r1_qv_series.index, r1_qv_series.values, color='blue', alpha=0.7)
    ax1.set_title("R1 QV Frequencies")
    ax1.set_xlabel("Quality Value")
    ax1.set_ylabel("Frequency")
    
    # R2 Histogram
    ax2.bar(r2_qv_series.index, r2_qv_series.values, color='green', alpha=0.7)
    ax2.set_title("R2 QV Frequencies")
    ax2.set_xlabel("Quality Value")
    ax2.set_ylabel("Frequency")
    
    # Save the plot to the provided output path
    plt.tight_layout()
    plt.savefig(output_png_path)
    plt.close()


# Test the function
output_path = "/mnt/data/quality_plot.png"
generate_qc_plot('/mnt/data/Pogona_combined_ILM.1.fastq_Pogona_combined_ILM.2.fastq_QC.csv', output_path)

output_path
