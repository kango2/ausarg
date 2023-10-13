#!/bin/bash
module load biopython

# Usage: 
# ./asm_qc.sh -f <fasta_file> -i <illumina_data> -p <pacbio_data> -o <ont_data> -d <output_directory> -l <log_directory> -P <project> -s <storage>

usage() {
    echo "Usage: $0 [-f FASTA_FILE] [-i ILLUMINA_DATA] [-p PACBIO_DATA] [-o ONT_DATA] [-d OUTPUT_DIRECTORY] [-l LOG_DIRECTORY] [-P PROJECT] [-s STORAGE]"
    echo 
    echo "Required arguments:"
    echo "  -f    Path to FASTA file"
    echo "  -i    Path to Illumina data. For paired-end data, separate the forward and reverse reads with a ';'. Different samples should be separated with a ':'. Example: sample1_R1.fastq;sample1_R2.fastq:sample2_R1.fastq;sample2_R2.fastq"
    echo "  -p    Path to PacBio data. Multiple files should be separated with a ':'. Example: sample1_pacbio.fastq:sample2_pacbio.fastq"
    echo "  -o    Path to ONT data. Multiple files should be separated with a ':'. Example: sample1_ont.fastq:sample2_ont.fastq"
    echo "  -d    Path to the output directory"
    echo "  -l    Path to the log directory"
    echo "  -P    Project name"
    echo "  -s    Storage option for qsub. Eg gdata/xl04+gdata/if89"
    echo
    echo "Description:"
    echo "This script performs various assembly quality control and data processing operations. It computes GC content, finds telomeres and centromeres, generates read depths for different platforms, and converts an assembly to a sequence table."
    exit 1
}

while getopts ":f:i:p:o:d:l:P:s:" opt; do
    case ${opt} in
        f ) fasta=$OPTARG ;;
        i ) illumina=$OPTARG ;;
        p ) pacbio=$OPTARG ;;
        o ) ont=$OPTARG ;;
        d ) outputdir=$OPTARG ;;
        l ) logsdir=$OPTARG ;;
        P ) project=$OPTARG ;;
        s ) storage=$OPTARG ;;
        \? ) usage ;;
    esac
done

if [ -z "${fasta}" ] || [ -z "${illumina}" ] || [ -z "${pacbio}" ] || [ -z "${ont}" ] || [ -z "${outputdir}" ] || [ -z "${logsdir}" ] || [ -z "${project}" ] || [ -z "${storage}" ]; then
    usage
fi

ref_base=$(basename "${fasta}" .fasta)
ref=${fasta}
gc=/g/data/xl04/ka6418/github/ausarg/scripts/gc_content.sh 
telomere=/g/data/xl04/ka6418/github/ausarg/scripts/find_telomeres.sh
centromere=/g/data/xl04/ka6418/github/ausarg/scripts/centromeres.sh
read_depth=/g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth_launcher_v2.sh
seq_table=/g/data/xl04/ka6418/github/ausarg/scripts/asm_to_sequencetable.py

echo qsub -P ${project} -l "storage=${storage}" -N "${ref_base}_GC" -j oe -o ${logsdir} -v input=${fasta},output=${outputdir} ${gc} 

echo qsub -P ${project} -l "storage=${storage}" -N "${ref_base}_Tel" -j oe -o ${logsdir} -v input=${fasta},output=${outputdir},permatch=90,copies=100 ${telomere}

echo qsub -P ${project} -l "storage=${storage}" -N "${ref_base}_Cen" -j oe -o ${logsdir} -v inputfasta=${fasta},outputdir=${outputdir} ${centromere}

mkdir -p ${outputdir}/"${ref_base}_depth"

qsub -q express -N "${ref_base}_illumina" -o ${logsdir} -P ${project} -l "storage=${storage}" -v platform=illumina,rawreads=${illumina},reference=${ref},output=${outputdir}/"${ref_base}_depth"/${ref_base}_illumina /g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth.sh
qsub -q express -N "${ref_base}_ont" -o ${logsdir} -P ${project} -l "storage=${storage}" -v platform=ont,rawreads=${ont},reference=${ref},output=${outputdir}/"${ref_base}_depth"/${ref_base}_ont /g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth.sh
qsub -q express -N "${ref_base}_pacbio" -o ${logsdir} -P ${project} -l "storage=${storage}" -v platform=pacbio,rawreads=${pacbio},reference=${ref},output=${outputdir}/"${ref_base}_depth"/${ref_base}_pacbio /g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth.sh

#python3 ${seq_table} -fasta ${fasta} -outputdir ${outputdir}