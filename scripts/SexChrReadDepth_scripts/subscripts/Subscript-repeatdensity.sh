# scripts to generate the repeat density input, can run manually as they are not computationally intensive
module use /g/data/if89/apps/modulefiles
module load bedtools/2.31.0 seqtk/1.4

seqtk comp /path/to/genome.fasta | cut -f1,2 > BASDU.tsv
# one for repeat density in 1MB windows and one for 20KB windows
bedtools makewindows -g BASDU.tsv -w 1000000 > BASDU_1MBwindows.bed
bedtools makewindows -g BASDU.tsv -w 20000 > BASDU_20KBwindows.bed
bedtools sort -i /path/to/repeatmasker/output.gff | bedtools merge > BASDU.bed3
bedtools coverage -a BASDU_1MBwindows.bed -b BASDU.bed3 > BASDU_1MBcoverage.bed
bedtools coverage -a BASDU_20KBwindows.bed -b BASDU.bed3 > BASDU_20KBcoverage.bed
# coverage (in % of bases) in this case is the number of bases in the 1MB window that are covered by the repeatmasker gff features, which are repeats, so it will range from 0-1, nothing masked or 100$ masked
