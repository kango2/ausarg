#!/bin/bash -ue
module load pythonlib 
python3 /g/data/xl04/ka6418/github/ausarg/scripts/asm_to_sequencetable.py -fasta /g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta -outputdir /g/data/xl04/ka6418/nextflow_testing/pipelinetest
