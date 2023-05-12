import json
import csv
import os

# Load the JSON data from a file
json_file = '/g/data/xl04/ka6418/species/Tiliqua_Rugosa/rTilRug0.5/assembly/evaluation/busco/run_vertebrata_odb10/short_summary.json'

with open(json_file, 'r') as f:
    data = json.load(f)

results = data["results"]

# Add the "BUSCO_" prefix to the keys and skip the first and last columns
keys_to_use = list(results.keys())[1:-1]
busco_data = {"BUSCO_" + key: results[key] for key in keys_to_use}

# Define your existing CSV file path
csv_file = '/g/data/xl04/ka6418/ausarg/experimental_scripts/assembly_table/rTilRug0.5.asm.hp1_stats.csv'

# Read the existing CSV file and store the header and rows
header = []
rows = []

with open(csv_file, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)

    for row in reader:
        rows.append(row)

# Add the BUSCO data to the header and rows
header += list(busco_data.keys())
rows = [row + list(busco_data.values()) for row in rows]

# Write the updated data to the existing CSV file
with open(csv_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)

    for row in rows:
        writer.writerow(row)
