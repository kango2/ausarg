#!/bin/bash -ue
pacbio_joined=$(echo ${longread_pb.fq.gz longread_pb_two.fq.gz.join(':')})
echo $pacbio_joined > /g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/test.txt
