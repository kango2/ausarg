#!/bin/bash

#TODO : functionality to launch separate jobs if needed only for one tech
#TODO : custom job names like ref_tech_readdepth

logsdir=${logsdir}
illumina=${illumina}
pacbio=${pacbio}
ont=${ont}
ref=${ref}
outputdir=${outputdir}
ref_base=$(basename ${ref}.fasta)
mkdir -p ${outputdir}/${ref_base}


echo qsub -o ${logsdir} -l storage=gdata/xl04+gdata/if89 -v platform=illumina,rawreads=${illumina},reference=${ref},output=${outputdir}/${ref_base}/${ref_base}_illumina /g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth.sh
echo qsub -o ${logsdir} -l storage=gdata/xl04+gdata/if89 -v platform=ont,rawreads=${ont},reference=${ref},output=${outputdir}/${ref_base}/${ref_base}_ont /g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth.sh
echo qsub -o ${logsdir} -l storage=gdata/xl04+gdata/if89 -v platform=pacbio,rawreads=${pacbio},reference=${ref},output=${outputdir}/${ref_base}/${ref_base}_pacbio /g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth.sh
