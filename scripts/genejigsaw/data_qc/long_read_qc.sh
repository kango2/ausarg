#!/bin/bash

source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate nanoplot_env

pacbio_fastq=/g/data/xl04/bpadata/Bassiana_duperreyi/raw/pacbio/350730_AusARG_AGRF_PacBio_DA060219.fq.gz:/g/data/xl04/bpadata/Bassiana_duperreyi/raw/pacbio/350730_AusARG_AGRF_PacBio_DA060220.fq.gz
ont_fastq=/g/data/xl04/ausarg_data/Bassiana_duperreyi/raw/ont/350748_AusARG_RamaciottiGarvan_ONTPromethION_PAF33103.fq.gz:/g/data/xl04/ausarg_data/Bassiana_duperreyi/raw/ont/350765_AusARG_RamaciottiGarvan_ONTPromethION_PAF30832.fq.gz:/g/data/xl04/ausarg_data/Bassiana_duperreyi/raw/ont/350766_AusARG_RamaciottiGarvan_ONTPromethION_PAG07318.fq.gz:/g/data/xl04/ausarg_data/Bassiana_duperreyi/raw/ont/350767_AusARG_RamaciottiGarvan_ONTPromethION_PAG18256.fq.gz:/g/data/xl04/ausarg_data/Bassiana_duperreyi/raw/ont/350782_AusARG_RamaciottiGarvan_ONTPromethION_PAG18312.fq.gz
job_name=qc
outdir=/g/data/xl04/ka6418/bassiana/raw_data_evaluation
storage=gdata/xl04+gdata/if89
logs=/g/data/xl04/ka6418/bassiana/raw_data_evaluation
project=xl04

mkdir -p ${outdir}/pacbio
mkdir -p ${outdir}/ont

# Convert colon-separated list to array for pacbio_fastq
IFS=':' read -ra pacbio_files <<< "$pacbio_fastq"
for file in "${pacbio_files[@]}"; do
    qsub -P ${project} -N ${job_name} -l storage=${storage} -v outdir=${outdir}/pacbio,fastq=$file /g/data/xl04/ka6418/github/ausarg/scripts/nanoplot.sh
done

# Convert colon-separated list to array for ont_fastq
IFS=':' read -ra ont_files <<< "$ont_fastq"
for file in "${ont_files[@]}"; do
    qsub -P ${project} -N ${job_name} -l storage=${storage} -v outdir=${outdir}/ont,fastq=$file /g/data/xl04/ka6418/github/ausarg/scripts/nanoplot.sh
done
