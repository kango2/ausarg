#!/bin/bash -ue
for file in longread_pb.fq.gz longread_pb_two.fq.gz; do
    echo $file >> /g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/test.txt
done
