#!/bin/bash -ue
trim_galore  -o shortread_qc --cores ${PBS_NCPUS} --paired test_illum_R1.fastq.gz test_illum_R2.fastq.gz
