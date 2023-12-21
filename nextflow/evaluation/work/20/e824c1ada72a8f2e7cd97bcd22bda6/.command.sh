#!/bin/bash -ue
export inputfasta=/g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta
export outputdir=/g/data/xl04/ka6418/nextflow_testing/pipelinetest

bash /g/data/xl04/ka6418/github/ausarg/scripts/centromeres.sh
