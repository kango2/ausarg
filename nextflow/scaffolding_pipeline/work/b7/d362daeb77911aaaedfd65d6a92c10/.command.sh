#!/bin/bash -ue
module load yahs

yahs -e GATC,GANTC,CTNAG,TTAA -o /g/data/xl04/ka6418/nextflow_testing/dedup/dedupp /g/data/xl04/ka6418/nextflow_testing/testdata/asm.contigs.fasta /g/data/xl04/ka6418/nextflow_testing/arimatest/asm.contigs_ArimaHiC.bam
# -r 10000,20000,50000,100000,200000,500000,1000000
