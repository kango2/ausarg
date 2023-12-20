#!/bin/bash -ue
echo longread_ont.fq.gz >> /g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/test.txt



joined_files=""

for file in longread_ont.fq.gz; do
    if [ -z "$joined_files" ]; then
        joined_files="$file"
    else
        joined_files="${joined_files}:${file}"
    fi
done

echo /g/data/xl04/ka6418/github/ausarg/nextflow/kmer_nf.sh -i ${joined_files} -s BASDU -o /g/data/xl04/ka6418/github/ausarg/nextflow/outtest -l 17 -t ONT >> /g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/test.txt
