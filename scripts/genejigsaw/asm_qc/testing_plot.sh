#!/bin/bash
#PBS -N alignment
#PBS -P xl04
#PBS -q express
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l walltime=05:00:00
#PBS -l mem=32GB
#PBS -l ncpus=8
#PBS -l wd

module load Rlib

Rscript /g/data/xl04/ka6418/github/ausarg/scripts/genejigsaw/asm_qc/asm_qc_plot.R --seqtable "/g/data/xl04/ka6418/chromosome_graph/rTilRug_HiC_pctg_seqtable.csv" --telomeres "/g/data/xl04/ka6418/chromosome_graph/rTilRug_HiC_pctg_Telomeres.csv" --repetitive "/g/data/xl04/ka6418/chromosome_graph/TRASH_output/Summary.of.repetitive.regions.rTilRug_HiC_pctg.fasta.csv" --ilmnrd "/g/data/xl04/ka6418/chromosome_graph/rTilRug_HiC_pctg_Illumina_sorted.bam.binned.fixed.depth.csv" --ontrd "/g/data/xl04/ka6418/chromosome_graph/rTilRug_HiC_pctg_ONT_sorted.bam.binned.fixed.depth.csv" --pacbio "/g/data/xl04/ka6418/chromosome_graph/rTilRug_HiC_pctg_PacBio_sorted.bam.binned.fixed.depth.csv" --gc "/g/data/xl04/ka6418/chromosome_graph/rTilRug_HiC_pctg.fasta_GC.csv"
