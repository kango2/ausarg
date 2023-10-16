#!/bin/bash
source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate nanoplot_env

pacbio_fastq=
ont_fastq=
illumina_fastq=
illumina_cleaneddir=
job_name=
outdir=
storage=
logs=
project=

mkdir -p ${outdir}/pacbio
mkdir -p ${outdir}/ont
mkdir -p ${outdir}/illumina
mkdir -p ${outdir}/illumina/beforetrim
mkdir -p ${outdir}/illumina/aftertrim
mkdir -p ${illumina_cleaneddir}

for every file in ${pacbio_fastq}:
        qsub -v outdir=${outdir}/pacbio,fastq=file nanoplot.sh

for every file in ${ont_fastq}:
        qsub -v outdir=${outdir}/ont,fastq=file nanoplot.sh

for every file in ${illumina_fastq}:
        qsub -v outdir=${outdir}/illumina/beforetrim,fastq=file fastqc.sh
    
for every file in ${illumina_fastq}:
        qsub -v outdir_trim=${illumina_cleaneddir},outdir_qc=${outdir}/illumina/aftertrim fastq=file trimgalore.sh

