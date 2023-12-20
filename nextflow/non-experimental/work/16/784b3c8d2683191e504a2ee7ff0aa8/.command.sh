#!/bin/bash -ue
/g/data/xl04/ka6418/bassiana/hifiasm_bassiana/hifiasm/hifiasm -t ${PBS_NCPUS} -o "/g/data/xl04/ka6418/nextflow_testing/pipelinetest/assembly/BASDU"  longread_pb.fq.gz longread_pb_two.fq.gz 
# --ul longread_ont.fq.gz --h1 HIC_R1.fastq.gz PE_1.fq.gz --h2 HIC_R2.fastq.gz PE_2.fq.gz 

awk '/^S/{print ">"$2;print $3}' /g/data/xl04/ka6418/nextflow_testing/pipelinetest/assembly/BASDU*.p_ctg.gfa > /g/data/xl04/ka6418/nextflow_testing/pipelinetest/assembly/BASDU_HifiASM.fasta
