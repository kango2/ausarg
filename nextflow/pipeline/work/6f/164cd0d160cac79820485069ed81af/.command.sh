#!/bin/bash -ue
module load biopython/1.79

python3 /g/data/xl04/ka6418/github/ausarg/nextflow/long-read-qv/long-read-qv.py --i testhifi.fq.gz --o pipeline --s BASDU --p PACBIO_SMRT --f PAF987
