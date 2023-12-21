import pandas as pd
import matplotlib.pyplot as plt
import os

def generate_plots_from_csv(file_path,output_dir):
    # Load the CSV file
    data = pd.read_csv(file_path)

    # Mean QV at Read Position
    data[['Mean_QV_at_read_position_1', 'Mean_QV_at_read_position_2']] = data['Mean_QV_at_read_position'].str.split(';', n=1, expand=True)
    mean_qv_1 = list(map(float, data['Mean_QV_at_read_position_1'][0].split(',')))
    mean_qv_2 = list(map(float, data['Mean_QV_at_read_position_2'][0].split(',')))

    # Nucleotide Counts
    nucleotide_counts = data['Nucleotide_count_at_read_position'].iloc[0]
    r1_counts, r2_counts = nucleotide_counts.split(';')
    r1_counts_list = [list(map(int, position.split(':'))) for position in r1_counts.split(',')]
    r2_counts_list = [list(map(int, position.split(':'))) for position in r2_counts.split(',')]
    r1_df = pd.DataFrame(r1_counts_list, columns=['A', 'C', 'G', 'T', 'N'])
    r2_df = pd.DataFrame(r2_counts_list, columns=['A', 'C', 'G', 'T', 'N'])

    # Total Nucleotide Content
    nucleotide_content = data['Nucleotide_content'].iloc[0]
    r1_content, r2_content = nucleotide_content.split(';')
    r1_content_list = [list(map(int, position.split(':'))) for position in r1_content.split(',')]
    r2_content_list = [list(map(int, position.split(':'))) for position in r2_content.split(',')]
    r1_content_df = pd.DataFrame(r1_content_list, columns=['A', 'C', 'G', 'T', 'N'])
    r2_content_df = pd.DataFrame(r2_content_list, columns=['A', 'C', 'G', 'T', 'N'])
    r1_totals = r1_content_df.sum()
    r2_totals = r2_content_df.sum()

    # Mean QV Per Read
    mean_qv_per_read = data['Mean_QV_per_read'].iloc[0]
    r1_qv, r2_qv = mean_qv_per_read.split(';')
    r1_qv_list = list(map(int, r1_qv.split(',')))
    r2_qv_list = list(map(int, r2_qv.split(',')))

    # Plotting
    fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(18, 24))

    # Mean QV at Read Position plots
    axes[0, 0].plot(mean_qv_1, color='cornflowerblue', alpha=0.8)
    axes[0, 0].set_xlabel('Read Position')
    axes[0, 0].set_ylabel('Mean QV')
    axes[0, 0].set_ylim(25, 40)
    axes[0, 0].set_title('Mean QV at Read Position (R1)')

    axes[0, 1].plot(mean_qv_2, color='cornflowerblue', alpha=0.8)
    axes[0, 1].set_xlabel('Read Position')
    axes[0, 1].set_ylim(25, 40)
    axes[0, 1].set_title('Mean QV at Read Position (R2)')

    # Nucleotide Counts plots
    for nucleotide in ['A', 'C', 'G', 'T', 'N']:
        axes[1, 0].plot(r1_df[nucleotide], label=nucleotide, alpha=0.6)
    axes[1, 0].set_title('Nucleotide Counts for R1 File')
    axes[1, 0].set_xlabel('Read Position')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].legend()
    axes[1, 0].grid(True)

    for nucleotide in ['A', 'C', 'G', 'T', 'N']:
        axes[1, 1].plot(r2_df[nucleotide], label=nucleotide, alpha=0.6)
    axes[1, 1].set_title('Nucleotide Counts for R2 File')
    axes[1, 1].set_xlabel('Read Position')
    axes[1, 1].legend()
    axes[1, 1].grid(True)

    # Total Nucleotide Content plots
    colors = {'A': '#add8e6', 'C': '#90ee90', 'G': '#ffcccb', 'T': '#dda0dd', 'N': '#ffd700'}
    
    bottom_r1 = 0
    bottom_r2 = 0
    for nucleotide, color in colors.items():
        axes[2, 0].bar(1, r1_totals[nucleotide], bottom=bottom_r1, color=color, edgecolor='white', label=nucleotide)
        bottom_r1 += r1_totals[nucleotide]

        axes[2, 1].bar(1, r2_totals[nucleotide], bottom=bottom_r2, color=color, edgecolor='white', label=nucleotide)
        bottom_r2 += r2_totals[nucleotide]

    axes[2, 0].set_title('Total Nucleotide Content for R1 File')
    axes[2, 0].set_ylabel('Total Nucleotide Content')
    axes[2, 0].legend()
    axes[2, 0].get_xaxis().set_visible(False)
    axes[2, 1].get_xaxis().set_visible(False)

    axes[2, 1].set_title('Total Nucleotide Content for R2 File')
    axes[2, 1].legend()

    # Mean QV Per Read plots
    axes[3, 0].bar(range(61), r1_qv_list[:61], color='cornflowerblue')
    axes[3, 0].set_title('Mean QV Per Read for R1 File (QV 0 to 60)')
    axes[3, 0].set_xlabel('QV')
    axes[3, 0].set_ylabel('Count of Reads')
    axes[3, 0].set_xlim(0, 60)

    axes[3, 1].bar(range(61), r2_qv_list[:61], color='cornflowerblue')
    axes[3, 1].set_title('Mean QV Per Read for R2 File (QV 0 to 60)')
    axes[3, 1].set_xlabel('QV')
    axes[3, 1].set_ylabel('Count of Reads')
    axes[3, 1].set_xlim(0, 60)

    plt.tight_layout()
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Saving the plot as a PDF in the specified output directory
    base_name = os.path.basename(file_path)  # Extract base name of the file
    pdf_filename = os.path.splitext(base_name)[0] + '.pdf'  # Replace extension with .pdf
    output_path = os.path.join(output_dir, pdf_filename)  # Full path for the output file
    plt.savefig(output_path, format='pdf')  # Save the figure as a PDF

    plt.show()


# Path to the CSV file (replace with actual file path if needed)
csv_file_path = "/g/data/xl04/bpadata/Bassiana_duperreyi/raw/evaluation/illumina_qc/350747_AusARG_UNSW_HTYH7DRXX_GTATTCCACC-TTGTCTACAT_S2_L001_R1_001.fastq_350747_AusARG_UNSW_HTYH7DRXX_GTATTCCACC-TTGTCTACAT_S2_L001_R2_001.fastq_QC.csv"  # Replac

generate_plots_from_csv(csv_file_path,"/g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental")