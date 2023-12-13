#!/bin/bash -ue
for file in [/g/data/xl04/ka6418/nextflow_testing/testdata/longread_pb.fq.gz, /g/data/xl04/ka6418/nextflow_testing/testdata/longread_pb_two.fq.gz]; do
    echo $file >> /g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/test.txt
done
