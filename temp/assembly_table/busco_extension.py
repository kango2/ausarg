import json
import csv
import argparse

# Create the argument parser
parser = argparse.ArgumentParser(description='Update CSV file with data from JSON file')
parser.add_argument('--json_file', type=str, help='Path to the JSON file')
parser.add_argument('--csv_file', type=str, help='Path to the CSV file')

# Parse the command-line arguments
args = parser.parse_args()

# Load the JSON data from the file
with open(args.json_file, 'r') as f:
    data = json.load(f)

results = data["results"]

# Add the "BUSCO_" prefix to the keys and skip the first and last columns
keys_to_use = list(results.keys())[1:-1]
busco_data = {"BUSCO_" + key: results[key] for key in keys_to_use}

# Read the existing CSV file and store the header and rows
header = []
rows = []

with open(args.csv_file, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)

    for row in reader:
        rows.append(row)

# Add the BUSCO data to the header and rows
header += list(busco_data.keys())
rows = [row + list(busco_data.values()) for row in rows]

# Write the updated data to the existing CSV file
with open(args.csv_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)

    for row in rows:
        writer.writerow(row)
