#!/bin/bash -ue
source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate genejigsaw
trim_galore -o 'shortread_qc' --cores ${PBS_NCPUS} --paired 'HIC_R1.fastq.gz' 'HIC_R2.fastq.gz'
