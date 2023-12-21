#!/bin/bash -ue
joined_files=""

for file in HIC_R1.fastq.gz HIC_R2.fastq.gz; do
    if [ -z "$joined_files" ]; then
        joined_files="$file"
    else
        joined_files="${joined_files}:${file}"
    fi
done

/g/data/xl04/ka6418/github/ausarg/nextflow/kmer_nf.sh -i ${joined_files} -s BASDU -o /g/data/xl04/ka6418/github/ausarg/nextflow/outtest/rawdata/kmers -l 21 -t Illumina
