#!/bin/bash
#PBS -N deepconsensus
#PBS -P xl04
#PBS -q normal
#PBS -l storage=gdata/xl04+gdata/if89+gdata/te53
#PBS -l walltime=48:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
#PBS -l wd
#PBS -j oe

#sampleid to runid for NCBI nomenclature

module load nextflow

cd ${outdir}

rm -f "${outdir}/${sample}_report.html" "${outdir}/${sample}_trace.txt" "${outdir}/${sample}_timeline.html" "${outdir}/${sample}_graph.dot"

nextflow -log ${outdir}/${sample}.log \
 run /g/data/xl04/ka6418/github/ausarg/workflows/deepconsensus/deepconsensus.nf \
 -work-dir ${outdir}/${sample}_workdir \
 -config /g/data/xl04/ka6418/github/ausarg/workflows/deepconsensus/deepconsensus_prod.config \
 -profile NCI \
 --subreads ${subreads} \
 --output ${outdir} \
 --chunks ${chunks} \
 --sample ${sample} \
 --workdir ${sampledir} \
 -with-report "${outdir}/${sample}_report.html" \
 -with-trace "${outdir}/${sample}_trace.txt"  \
 -with-timeline "${outdir}/${sample}_timeline.html" \
 -with-dag "${outdir}/${sample}_graph.dot" \
 -resume