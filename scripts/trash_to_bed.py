import pandas as pd
import argparse

def main(input_file, output_file):
    # Load the CSV file
    df = pd.read_csv(input_file)

    # Filter the dataframe for rows where 'class' exists and 'width' is more than 100,000
    filtered_df = df.dropna(subset=['class'])

    # Assuming 'name' is the chromosome, 'start' and 'end' are positions,
    # and 'class' can be used as the name of the feature in the BED file.
    bed_df = filtered_df[['name', 'start', 'end', 'class']]

    # Save the filtered dataframe to a BED file
    bed_df.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Filter CSV and convert to BED format.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input CSV file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output BED file.')

    # Parse the arguments
    args = parser.parse_args()

    # Run the main function with the provided arguments
    main(args.input, args.output)
