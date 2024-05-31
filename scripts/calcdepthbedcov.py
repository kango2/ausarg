import pandas as pd
import argparse
import os
import concurrent.futures

def calculate_mean_median(scaffold, group):
    mean = round(group['read_depth'].mean(), 2)
    median = round(group['read_depth'].median(), 2)
    length = group['end'].max() - group['start'].min()
    return scaffold, mean, median, length

def calculate_read_depths(input_file, output_file):
    # Read the BED file into a DataFrame
    df = pd.read_csv(input_file, sep="\s+", header=None, names=['scaffold', 'start', 'end', 'read_depth'])

    # Group by scaffold
    grouped = df.groupby('scaffold')

    # Use multithreading to calculate mean and median
    results = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(calculate_mean_median, scaffold, group) for scaffold, group in grouped]
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

    # Convert results to DataFrame, sort by length, and save to file
    result_df = pd.DataFrame(results, columns=['scaffold', 'mean', 'median', 'length'])
    result_df = result_df.sort_values(by='length', ascending=False)
    result_df.to_csv(output_file, sep="\t", index=False, columns=['scaffold', 'mean', 'median'])

def main():
    parser = argparse.ArgumentParser(description='Calculate mean and median read depth per scaffold.')
    parser.add_argument('-i', '--input', required=True, help='Input BED file')
    parser.add_argument('-o', '--output', help='Output file (optional)')

    args = parser.parse_args()

    # Determine output file name
    if args.output:
        output_file = args.output
    else:
        base_name = os.path.basename(args.input)
        output_file = os.path.splitext(base_name)[0] + '_perscaf.tsv'

    calculate_read_depths(args.input, output_file)
    print(f'Results saved to {output_file}')

if __name__ == "__main__":
    main()
