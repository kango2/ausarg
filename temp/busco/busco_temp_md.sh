#!/bin/bash
module load singularity

# Set default values for variables
FASTA_FILE=""
BUSCO_IMG=""
LINEAGE_DIR=""
OUTPUT_DIR=""
PBS_NCPUS=""
MODE="genome"  # Default mode

# Parse command-line arguments
while getopts ":f:i:l:o:c:m:" opt; do
  case $opt in
    f)
      FASTA_FILE=$OPTARG
      ;;
    i)
      BUSCO_IMG=$OPTARG
      ;;
    l)
      LINEAGE_DIR=$OPTARG
      ;;
    o)
      OUTPUT_DIR=$OPTARG
      ;;
    c)
      PBS_NCPUS=$OPTARG
      ;;
    m)
      MODE=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Check if required variables are provided
if [ -z "$FASTA_FILE" ] || [ -z "$BUSCO_IMG" ] || [ -z "$LINEAGE_DIR" ]; then
  echo "Error: Missing required input. Usage: $0 -f <FASTA_FILE> -i <BUSCO_IMG> -l <LINEAGE_DIR> [-o <OUTPUT_DIR>] [-c <PBS_NCPUS>] [-m <MODE>]"
  exit 1
fi

# Set default values if optional parameters are not provided
if [ -z "$OUTPUT_DIR" ]; then
  OUTPUT_DIR="$(dirname "$FASTA_FILE")/evaluation/busco"
fi

if [ -z "$PBS_NCPUS" ]; then
  PBS_NCPUS=""
fi

# Execute the BUSCO command with PBS directives
singularity exec $BUSCO_IMG busco \
  -i $FASTA_FILE \
  -o $OUTPUT_DIR \
  -l $LINEAGE_DIR \
  -m $MODE \
  $PBS_NCPUS \
  --offline

