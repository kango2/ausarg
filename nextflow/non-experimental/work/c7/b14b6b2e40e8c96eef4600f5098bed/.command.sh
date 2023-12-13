#!/bin/bash -ue
echo /g/data/xl04/ka6418/nextflow_testing/testdata/fastq/HIC_R1.fastq.gz:/g/data/xl04/ka6418/nextflow_testing/testdata/fastq/HIC_R2.fastq.gz:/g/data/xl04/ka6418/nextflow_testing/testdata/PE_1.fq.gz:/g/data/xl04/ka6418/nextflow_testing/testdata/PE_1.fq.gz >> /g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/test.txt
if ["ILLUMINA" = "ILLUMINA"]; then 

    echo /g/data/xl04/ka6418/github/ausarg/nextflow/kmer_nf.sh -i /g/data/xl04/ka6418/nextflow_testing/testdata/fastq/HIC_R1.fastq.gz:/g/data/xl04/ka6418/nextflow_testing/testdata/fastq/HIC_R2.fastq.gz:/g/data/xl04/ka6418/nextflow_testing/testdata/PE_1.fq.gz:/g/data/xl04/ka6418/nextflow_testing/testdata/PE_1.fq.gz -s BASDU -o /g/data/xl04/ka6418/github/ausarg/nextflow/outtest -l 17 -t ILLUMINA >> /g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/test.txt

else 

    joined_files=""

    for file in /g/data/xl04/ka6418/nextflow_testing/testdata/fastq/HIC_R1.fastq.gz:/g/data/xl04/ka6418/nextflow_testing/testdata/fastq/HIC_R2.fastq.gz:/g/data/xl04/ka6418/nextflow_testing/testdata/PE_1.fq.gz:/g/data/xl04/ka6418/nextflow_testing/testdata/PE_1.fq.gz; do
        if [ -z "$joined_files" ]; then
            joined_files="$file"
        else
            joined_files="${joined_files}:${file}"
        fi
    done

    echo /g/data/xl04/ka6418/github/ausarg/nextflow/kmer_nf.sh -i ${joined_files} -s BASDU -o /g/data/xl04/ka6418/github/ausarg/nextflow/outtest -l 17 -t ILLUMINA >> /g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/test.txt

fi
