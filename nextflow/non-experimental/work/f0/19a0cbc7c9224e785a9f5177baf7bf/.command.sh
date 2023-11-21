#!/bin/bash -ue
source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate genejigsaw
trim_galore -o 'shortread_qc' --cores ${PBS_NCPUS} --paired 'test_illum_R1.fastq.gz' 'test_illum_R2.fastq.gz'
