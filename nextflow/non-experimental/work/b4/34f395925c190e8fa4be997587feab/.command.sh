#!/bin/bash -ue
joined_files=""

for file in longread_pb.fq.gz longread_pb_two.fq.gz; do
    if [ -z "$joined_files" ]; then
        joined_files="$file"
    else
        joined_files="${joined_files}:${file}"
    fi
done

/g/data/xl04/ka6418/github/ausarg/nextflow/kmer_nf.sh -i ${joined_files} -s BASDU -o /g/data/xl04/ka6418/github/ausarg/nextflow/outtest -l 17 -t PacBio