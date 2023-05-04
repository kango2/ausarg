Usage Instructions for "pacbio_ont_alignment.pbs"

The script supports sequence alignment for PacBio/ONT files with merge support for multiple files. You can supply the file names as a separate text file (query_fastq_list.txt included as an example in the directory)

An example qsub command to launch on Gadi : 
qsub -v REFERENCE_ASSEMBLY="/g/data/xl04/tb3184/assembly/rTilRug0.5/assembly.haplotype1.fasta",SEQUENCING_TECHNOLOGY="pacbio",DIRECTORY_PATH="/g/data/xl04/ka6418/sequence_alignment/bam_outputs/tiliqua_rugosa",QUERY_FASTQ_LIST_FILE="/g/data/xl04/ka6418/sequence_alignment/script/ausarg-script/query_fastq_list.txt" -N pacbio_ont_alignment.pbs pacbio_ont_alignment.pbs 