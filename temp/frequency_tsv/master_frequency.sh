#!/bin/bash

# Define a function to submit a job using qsub
submit_job() {
  local input_file="$1"
  local id="$2"
  local sequencing_technology="$3"
  local output_dir="$4"
  local job_name="job_$(basename "$input_file" .fq.gz)"
  local output_log="$LOG_DIR/${job_name}.out"
  local error_log="$LOG_DIR/${job_name}.err"

  qsub \
    -N "$job_name" \
    -o "$output_log" \
    -e "$error_log" \
    -v INPUT_FILE="$input_file",SEQ="$sequencing_technology",ODIR="$output_dir",ID="$id" \
    /g/data/xl04/ka6418/ausarg/temp/frequency_tsv/frequency_script.sh
}

# Set the directory where log files will be stored
LOG_DIR="/g/data/xl04/ka6418/logs"

# Create the log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Check if any input files were provided
if [ $# -eq 0 ]; then
  echo "Usage: $0 -f <input_file1,id1,sequencing_technology1,output_dir1> [-f <input_file2,id2,sequencing_technology2,output_dir2>] ..."
  exit 1
fi

# Initialize an array to store input groups
input_groups=()

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -f)
      shift
      input_groups+=("$1")
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Iterate through input groups and submit a job for each group
for input_group in "${input_groups[@]}"; do
  IFS=',' read -r input_file id sequencing_technology output_dir <<< "$input_group"
  
  # Call the function to submit the job for the current input group
  submit_job "$input_file" "$id" "$sequencing_technology" "$output_dir"
done

echo "Job submission completed."
