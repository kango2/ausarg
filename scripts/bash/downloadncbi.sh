#!/bin/bash
#PBS -N downloadncbi
#PBS -P xl04
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -l mem=2GB
#PBS -l ncpus=1
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -j oe

set -ex 

module load ncbi-datasets-cli/16.6.0 htslib utils

## Usage -v 
# outdir : output directory
# taxon : taxon name 

outfile_prefix="$taxon"

# Replace spaces in prefix with underscores
outfile_prefix=${outfile_prefix// /_}

random=${RANDOM}

# Generate JSON summary
datasets summary genome taxon "$taxon" --assembly-level chromosome --as-json-lines > ${outdir}/${outfile_prefix}.$(date +"%Y%m%d").${random}.json

# Convert JSON to TSV
dataformat tsv genome --fields accession,organism-common-name,organism-name,organism-tax-id,assminfo-status,assminfo-refseq-category \
  --inputfile ${outdir}/${outfile_prefix}.$(date +"%Y%m%d").${random}.json | grep current | grep 'reference genome' | cut -f1,3,4 > ${outdir}/${outfile_prefix}.$(date +"%Y%m%d").${random}.tsv 

rm ${outdir}/${outfile_prefix}.$(date +"%Y%m%d").${random}.json

# Process each entry in the TSV
while IFS=$'\t' read -r speciesid speciesname txid; do

    # Check if the directory exists
    if [ -d "${outdir}/${speciesid}.${speciesname// /_}.txid${txid}" ]; then
        echo "Directory ${outdir}/${speciesid}.${speciesname// /_}.txid${txid} already exists. Skipping..."
        continue
    fi 

    # Create directory and download genome
    mkdir -p "${outdir}/${speciesid}.${speciesname// /_}.txid${txid}"
    filename="${outdir}/${speciesid}.${speciesname// /_}.txid${txid}/${speciesid}.zip"
    datasets download genome accession ${speciesid} --include genome --filename ${filename}

    # Unzip and process the genome data
    cd "${outdir}/${speciesid}.${speciesname// /_}.txid${txid}"
    unzip ${filename}
    mv ncbi_dataset/data/assembly_data_report.jsonl "${speciesid}_data_report.jsonl"
    mv ncbi_dataset/data/${speciesid}/*.fna .

    rm -rf ncbi_dataset

done < ${outdir}/${outfile_prefix}.$(date +"%Y%m%d").${random}.tsv
