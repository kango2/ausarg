import sys
import os
import pandas as pd

def combine_csv_files(folder_path, output_file):
    # List to store dataframes from each CSV file
    dataframes = []

    # Iterate through CSV files in the folder
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.csv'):
            file_path = os.path.join(folder_path, file_name)
            df = pd.read_csv(file_path)
            dataframes.append(df)

    # Concatenate all dataframes into a single dataframe
    combined_df = pd.concat(dataframes)

    # Group by Read_Length and QV, then sum the Read_Numbers
    grouped_df = combined_df.groupby(['Read_Length', 'QV'], as_index=False)['Read_Numbers'].sum()

    # Save the final dataframe to the specified output CSV file
    grouped_df.to_csv(output_file, index=False)

    print("Combined CSV file saved:", output_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py folder_path output_file")
        sys.exit(1)

    folder_path = sys.argv[1]
    output_file = sys.argv[2]

    combine_csv_files(folder_path, output_file)

