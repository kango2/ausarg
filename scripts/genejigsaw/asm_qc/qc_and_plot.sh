#!/bin/bash
module load biopython

fasta=${fasta}
illumina=${illumina}
pacbio=${pacbio}
ont=${ont}
outputdir=${outputdir}
logsdir=${logsdir}
project=${project}
storage=${storage}

gc=/g/data/xl04/ka6418/github/ausarg/scripts/gc_content.sh 
telomere=/g/data/xl04/ka6418/github/ausarg/scripts/find_telomeres.sh
centromere=/g/data/xl04/ka6418/github/ausarg/scripts/centromeres.sh
read_depth=/g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth_launcher_v2.sh
seq_table=/g/data/xl04/ka6418/github/ausarg/scripts/asm_to_sequencetable.py

echo qsub -P ${project} -l ${storage} -N "GC" -j oe -o ${logsdir} -v input=${fasta},output=${outputdir} ${gc} 

echo qsub -P ${project} -l ${storage} -N "Tel" -j oe -o ${logsdir} -v input=${fasta},output=${outputdir},permatch=90,copies=100 ${telomere}

echo qsub -P ${project} -l ${storage} -N "Cen" -j oe -o ${logsdir} -v inputfasta=${fasta},outputdir=${outputdir} ${centromere}

echo logsdir=${logsdir} illumina=${illumina} pacbio=${pacbio} ont=${ont} ref=${fasta} outputdir=${outputdir} ${read_depth}  \

ref_base=$(basename "${ref}" .fasta)
mkdir -p ${outputdir}/${ref_base}

qsub -N "illumina" -o ${logsdir} -P ${project} -l ${storage} -v platform=illumina,rawreads=${illumina},reference=${ref},output=${outputdir}/${ref_base}/${ref_base}_illumina /g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth.sh
qsub -N "ont" -o ${logsdir} -P ${project} -l ${storage} -v platform=ont,rawreads=${ont},reference=${ref},output=${outputdir}/${ref_base}/${ref_base}_ont /g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth.sh
qsub -N "pacbio" -o ${logsdir} -P ${project} -l ${storage} -v platform=pacbio,rawreads=${pacbio},reference=${ref},output=${outputdir}/${ref_base}/${ref_base}_pacbio /g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth.sh


echo python3 ${seq_table} -fasta ${fasta} -outputdir ${outputdir}






