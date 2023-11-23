#!/bin/bash -ue
/g/data/xl04/ka6418/bassiana/hifiasm_bassiana/hifiasm/hifiasm -t ${PBS_NCPUS} -o "/g/data/xl04/ka6418/github/ausarg/nextflow/outtest/assembly/BASDU"  longread_pb.fq.gz 
# --ul longread_ont.fq.gz --h1 HIC_R1.fastq.gz --h2 HIC_R2.fastq.gz 

awk '/^S/{print ">"$2;print $3}' /g/data/xl04/ka6418/github/ausarg/nextflow/outtest/assembly/BASDU*.p_ctg.gfa > /g/data/xl04/ka6418/github/ausarg/nextflow/outtest/assembly/BASDU_HifiASM.fasta