import pandas as pd
import argparse
from collections import defaultdict

# Function to read AGP file and create a scaffold size dictionary
def read_agp_file(agp_file):
    scaffold_sizes = defaultdict(int)
    with open(agp_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            scaffold, scaffold_start, scaffold_end, *_ = parts
            scaffold_sizes[scaffold] = max(scaffold_sizes[scaffold], int(scaffold_end))
    return scaffold_sizes

# Function to create a dataframe with scaffold segments
def create_scaffold_dataframe(scaffold_sizes):
    data = []
    for scaffold, size in scaffold_sizes.items():
        for start in range(1, size, 1000):
            end = min(start + 999, size)
            data.append([scaffold, start, end, 'N/A'])
    return pd.DataFrame(data, columns=['Scaffold', 'Start', 'End', 'AvgDepth'])

# Function to parse the original depth file
def parse_depth_file(depth_file):
    depth_data = pd.read_csv(depth_file, header=None, names=['Scaffold', 'Start', 'End', 'AvgDepth'])
    return depth_data

# Function to adjust depth data based on AGP
def adjust_depth_data(scaffold_df, depth_data):
    for index, row in scaffold_df.iterrows():
        scaffold, start, end = row['Scaffold'], row['Start'], row['End']
        overlapping_depths = depth_data[(depth_data['Scaffold'] == scaffold) & (depth_data['End'] >= start) & (depth_data['Start'] <= end)]

        total_depth = 0
        total_length = 0
        for _, depth_row in overlapping_depths.iterrows():
            overlap_start = max(start, depth_row['Start'])
            overlap_end = min(end, depth_row['End'])
            overlap_length = overlap_end - overlap_start + 1
            total_depth += depth_row['AvgDepth'] * overlap_length
            total_length += overlap_length

        scaffold_df.at[index, 'AvgDepth'] = total_depth / total_length if total_length > 0 else 'N/A'

    return scaffold_df

# Main function to run the script
def main(agp_file, depth_file, output_file):
    scaffold_sizes = read_agp_file(agp_file)
    scaffold_df = create_scaffold_dataframe(scaffold_sizes)
    depth_data = parse_depth_file(depth_file)
    adjusted_df = adjust_depth_data(scaffold_df, depth_data)
    adjusted_df.sort_values(by=['Scaffold', 'Start', 'End'], inplace=True)
    adjusted_df.to_csv(output_file, index=False, header=False)
    print(f"Adjusted depth file saved to {output_file}")

# Command line argument parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adjust depth file based on AGP file.")
    parser.add_argument("-a", "--agp", required=True, help="Input AGP file")
    parser.add_argument("-d", "--depth", required=True, help="Input depth file")
    parser.add_argument("-o", "--output", required=True, help="Output file for adjusted depth")
    args = parser.parse_args()

    main(args.agp, args.depth, args.output)
