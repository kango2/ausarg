#!/bin/bash

# Define a function to submit a job using qsub
submit_job() {
  local input_file="$1"
  local job_name="job_$(basename "$input_file" .fq.gz)"
  local output_log="$LOG_DIR/${job_name}.out"
  local error_log="$LOG_DIR/${job_name}.err"

  qsub \
    -N "$job_name" \
    -o "$output_log" \
    -e "$error_log" \
    -v INPUT_FILE="$input_file" \
    /g/data/xl04/ka6418/jellyfish/dump_to_histo_slave.sh
}

# Set the directory where log files will be stored
LOG_DIR="/g/data/xl04/ka6418/logs"

# Create the log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Check if any input files were provided
if [ $# -eq 0 ]; then
  echo "Usage: $0 <input_file1> [<input_file2> ...]"
  exit 1
fi

# Iterate through each provided input file and submit a job
for input_file in "$@"; do
  # Call the function to submit the job for the current input file
  submit_job "$input_file"
done

echo "Job submission completed."
