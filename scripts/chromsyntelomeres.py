import pandas as pd
import argparse

# Function to process the telomere data
def process_telomeres(csv_path, fai_path, out_tsv):
    # Read the input CSV and FAI files
    telodf = pd.read_csv(csv_path)
    faidf = pd.read_csv(fai_path, sep='\t', usecols=[0, 1], names=['SeqName', 'SeqLen'], header=None)

    # Add new columns initialized to False
    faidf[['Tel5', 'Tel3']] = False
    faidf[['Tel5Len', 'Tel3Len', 'Trim5', 'Trim3', 'TelPerc']] = 0

    # Filter relevant columns from telodf
    telodf = telodf[['Sequence_ID', 'Relative Start']]

    # Process each row in telodf
    for index, row in telodf.iterrows():
        seq_name = row["Sequence_ID"]
        rel_start = row["Relative Start"]

        # Update faidf based on relative start values
        if 0 <= rel_start <= 0.1:
            faidf.loc[faidf["SeqName"] == seq_name, "Tel5"] = True
        elif 0.9 <= rel_start <= 1:
            faidf.loc[faidf["SeqName"] == seq_name, "Tel3"] = True

    # Save the updated DataFrame as a TSV file with headers
    faidf.to_csv(out_tsv, sep='\t', index=False)
    print(f"Processed file saved to {out_tsv}")

# Main function
if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process telomere data and update fasta index data.")
    parser.add_argument("-csv", required=True, help="Path to the input telomere CSV file.")
    parser.add_argument("-fai", required=True, help="Path to the input fasta index (FAI) file.")
    parser.add_argument("-outtsv", required=True, help="Path to save the output TSV file.")

    # Parse the arguments
    args = parser.parse_args()

    # Call the processing function
    process_telomeres(args.csv, args.fai, args.outtsv)