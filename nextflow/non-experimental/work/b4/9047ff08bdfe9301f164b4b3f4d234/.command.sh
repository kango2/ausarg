#!/bin/bash -ue
module load yahs samtools
samtools faidx /g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta
label=(basename /g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta .fasta)
yahs -e GATC,GANTC,CTNAG,TTAA -o /g/data/xl04/ka6418/github/ausarg/nextflow/outtest/yahs/$label /g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta /g/data/xl04/ka6418/github/ausarg/nextflow/outtest/arima/*_ArimaHiC.bam
mv /g/data/xl04/ka6418/github/ausarg/nextflow/outtest/yahs/${label}_scaffolds_final.fa /g/data/xl04/ka6418/github/ausarg/nextflow/outtest/yahs/${label}_YAHS.fasta
# -r 10000,20000,50000,100000,200000,500000,1000000
