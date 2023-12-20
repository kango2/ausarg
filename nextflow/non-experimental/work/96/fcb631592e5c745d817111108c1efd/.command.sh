#!/bin/bash -ue
source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate genejigsaw
trim_galore -o 'qc' --cores ${PBS_NCPUS} --paired 'PE_1.fq.gz' 'PE_2.fq.gz'
