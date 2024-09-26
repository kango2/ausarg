#!/bin/bash
#PBS -l ncpus=1,mem=80GB,walltime=00:30:00,storage=gdata/if89+gdata/xl04
#PBS -N PlotSexChrReadDepth
#PBS -j oe
#PBS -M z5205618@ad.unsw.edu.au
#PBS -m ae

module use /g/data/if89/apps/modulefiles
module load Rlib/4.3.1

# This plotting script should be run after you've generated the inputs using the subscripts
# Edit path to input files here after the export=
export chrLengths="/path/to/BASDU.tsv"
export chrFilter="/path/to/BASDU.filter"
export RD_male="/path/to/male_fixed.bedgraph"
export RD_female="/path/to/female_fixed.bedgraph"
export blast="/path/to/Ykmerblast.tsv"
export density="/path/to/repeat_density.bed"
export output="/path/to/output.pdf"

# You can manually generate the chrLengths file using this command
### module load seqtk/1.4
### seqtk comp /path/to/genome.fasta | cut -f1,2 | awk '{print $1"\t0\t"$2"\tA"NR"\tgneg"}' | sed $'1ichr\tstart\tend\tname\tgieStain' > /path/to/output/BASDU.tsv

# chrFilter should be a .lst file with the fasta headers of the chromosomes you want to plot
# RD_male should contain the fixed bamcoverage output from a male individual mapped to genome
# RD_female should contain the fixed bamcoverage output from a female individual mapped to genome
# blast should contain the output of the Y enriched kmers blased against the genome
# density should contain the repeat density obtained from repeatmasker gff
# output should be the path to the output pdf

Rscript Plotting.R \
${chrLengths} \
${chrFilter} \
${RD_male} \
${RD_female} \
${blast} \
${density} \
${output}
