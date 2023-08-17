#!/bin/bash

temp_dir="$1"
final_csv_output="$2"

csvs="$temp_dir/csvs"

python3 /g/data/xl04/ka6418/ausarg/temp/qv_frequency_concurrent/csv_combine.py "$csvs" "$final_csv_output"

# Remove the temporary directory
#rm -rf "$temp_dir"