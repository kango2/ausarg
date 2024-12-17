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

##Usage -v 
#outdir : output directory

if [[ -z "$outdir" ]]; then
    echo "Error: outdir variable is missing. Terminating script."
    exit 1
fi

datasets summary genome taxon Squamata --assembly-level chromosome --as-json-lines > ${outdir}/squamata.$(date +"%Y%m%d").json

dataformat tsv genome --fields accession,organism-common-name,organism-name,organism-tax-id,assminfo-status,assminfo-refseq-category \
  --inputfile ${outdir}/squamata.$(date +"%Y%m%d").json | grep current | grep 'reference genome' | cut -f1,3,4 > ${outdir}/squamata.$(date +"%Y%m%d").tsv 

rm ${outdir}/squamata.$(date +"%Y%m%d").json

while IFS=$'\t' read -r speciesid speciesname txid; do

#speciesid to asmaccession
#zip the genome
#runcmd.sh from util 
#add checkpoint

    if [ -d "${outdir}/${speciesid}.${speciesname// /_}.txid${txid}" ]; then
        echo "Directory ${dirpath} already exists. Skipping..."
        continue
    fi 

    mkdir -p "${outdir}/${speciesid}.${speciesname// /_}.txid${txid}"
    filename="${outdir}/${speciesid}.${speciesname// /_}.txid${txid}/${speciesid}.zip"
    datasets download genome accession ${speciesid} --include genome --filename ${filename}
    cd "${outdir}/${speciesid}.${speciesname// /_}.txid${txid}"
    unzip ${filename}
    mv ${outdir}/${speciesid}.${speciesname// /_}.txid${txid}/ncbi_dataset/data/assembly_data_report.jsonl "${outdir}/${speciesid}.${speciesname// /_}.txid${txid}/${speciesid}_data_report.jsonl"
    mv ${outdir}/${speciesid}.${speciesname// /_}.txid${txid}/ncbi_dataset/data/${speciesid}/*.fna "${outdir}/${speciesid}.${speciesname// /_}.txid${txid}"
    rm -rf ncbi_dataset

done < ${outdir}/squamata.$(date +"%Y%m%d").tsv 
