module load singularity

# Set default values for variables
FASTA_FILE=""
BUSCO_IMG=""
LINEAGE_DIR=""
OUTPUT_DIR=""

# Parse command-line arguments
while getopts ":f:i:l:" opt; do
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
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Check if required variables are provided
if [ -z "$FASTA_FILE" ] || [ -z "$BUSCO_IMG" ] || [ -z "$LINEAGE_DIR" ]; then
  echo "Error: Missing required input. Usage: $0 -f <FASTA_FILE> -i <BUSCO_IMG> -l <LINEAGE_DIR>"
  exit 1
fi

# Generate the output directory path automatically
OUTPUT_DIR="$(dirname "$FASTA_FILE")/evaluation/busco"


#Execute the BUSCO command with PBS directives
singularity exec $BUSCO_IMG busco \
  -i $FASTA_FILE \
  -o $OUTPUT_DIR \
  -l $LINEAGE_DIR \
  -m genome \
  -c $PBS_NCPUS \
  --offline