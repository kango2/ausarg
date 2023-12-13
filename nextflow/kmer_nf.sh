#!/bin/bash
set -ex

# Load required modules
module load jellyfish/2.3.0 utils/0.0

# Initialize default values
inputfiles=""
sampleID=""
OUTDIR=""
klength=""
technology=""

# Parse command line flags
while getopts "i:s:o:l:t:" opt; do
  case $opt in
    i) inputfiles=$OPTARG ;;
    s) sampleID=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    l) klength=$OPTARG ;;
    t) technology=$OPTARG ;;
    \?) echo "Invalid option -$OPTARG" >&2
        exit 1
        ;;
  esac
done

# Check if all parameters are set
if [ -z "$inputfiles" ] || [ -z "$sampleID" ] || [ -z "$OUTDIR" ] || [ -z "$klength" ]; then
    echo "Usage: $0 -i inputfiles -s sampleID -o OUTDIR -l klength"
    exit 1
fi

sampleID="${sampleID}_${technology}_${klength}"

IFS=":" read -ra filelist <<< "$inputfiles"
filecount=${#filelist[@]}

for file in "${filelist[@]}"; do
	#tried the following to increase efficiency but monitoring suggests that pigz uses <100% cpus indicating one thread only
	#TODO: use 3 or 4 threads and see if the CPU usage increases
  echo "pigz -c -d -p 2 ${file}"  >> ${PBS_JOBFS}/filegenerator.cmds
done

##jellyfish --text option will output text strings directly. So no need to create counts database and then dump
##dump files are about 80-90GB when quality filters are used for 1KGP data
##sort uses 150GB of virtual memory, temporary files in text format are written to disk that are about 19GB in size each
##TODO: dump files can be compressed after sort to save disk space
##specifying -s 32G reaches up to 142.5GB of RAM usage during the count step
##TODO: incorporate other jellyfish parameters for quality filter if required
##TODO: perhaps combine all .done files into one file sensibly

runcmd.sh -c "jellyfish count --text -g ${PBS_JOBFS}/filegenerator.cmds -G ${filecount} --threads=$((PBS_NCPUS - filecount * 2)) -m ${klength} -o ${PBS_JOBFS}/$sampleID.dump -C -s 32G" -t ${sampleID}.kmercount -d ${OUTDIR}/$sampleID.kmercount.done -f false
runcmd.sh -c "jellyfish histo -t ${PBS_NCPUS} -o ${OUTDIR}/$sampleID.histo ${PBS_JOBFS}/$sampleID.dump" -t ${sampleID}.kmerhisto -d ${OUTDIR}/$sampleID.kmerhisto.done -f false
runcmd.sh -c "sort -k1,1 --parallel=${PBS_NCPUS} --buffer-size=80% -T ${PBS_JOBFS} -o ${PBS_JOBFS}/$sampleID.sorted.dump ${PBS_JOBFS}/$sampleID.dump" -t ${sampleID}.kmersort -d ${OUTDIR}/$sampleID.kmersort.done -f false
runcmd.sh -c "pigz -q -p ${PBS_NCPUS} ${PBS_JOBFS}/$sampleID.sorted.dump" -t ${sampleID}.kmerzip -d ${OUTDIR}/$sampleID.kmerzip.done -f false
runcmd.sh -c "rsync -a ${PBS_JOBFS}/$sampleID.sorted.dump.gz ${OUTDIR}/" -t ${sampleID}.kmerxfer -d ${OUTDIR}/$sampleID.kmerxfer.done -f false                              