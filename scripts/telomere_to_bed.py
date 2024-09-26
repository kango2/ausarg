import os
import pandas as pd
import argparse

def main(input_file, output_dir):
    # Load the CSV file
    df = pd.read_csv(input_file)

    # Extract the base name of the input file (without extension)
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    
    # Create the output file path with .bed extension in the specified output directory
    output_file = os.path.join(output_dir, f"{base_name}.bed")

    # Select the columns needed for the BED format
    bed_df = df[['Sequence_ID', 'Start', 'End']]

    # Save the DataFrame to the output file in BED format (tab-separated, no header, no index)
    bed_df.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Filter CSV and convert to BED format.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input CSV file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output directory.')

    # Parse the arguments
    args = parser.parse_args()

    # Run the main function with the provided arguments
    main(args.input, args.output)
