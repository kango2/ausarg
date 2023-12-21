#!/bin/bash -ue
export input=/g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta
export output=pipelinetest
export permatch=90
export copies=100

bash /g/data/xl04/ka6418/github/ausarg/scripts/find_telomeres.sh
