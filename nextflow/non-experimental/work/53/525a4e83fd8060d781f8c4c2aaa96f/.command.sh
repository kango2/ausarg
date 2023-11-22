#!/bin/bash -ue
echo /g/data/xl04/ka6418/bassiana/hifiasm_bassiana/hifiasm/hifiasm -t ${PBS_NCPUS} -o "/g/data/xl04/ka6418/github/ausarg/nextflow/outtest/assembly/BASDU" --ul testont.fq.gz --h1 test_illum_R1.fastq.gz test_illumtwo_R1.fastq.gz --h2 test_illum_R2.fastq.gz test_illumtwo_R2.fastq.gz  testhifi.fq.gz testhifitwo.fq.gz
