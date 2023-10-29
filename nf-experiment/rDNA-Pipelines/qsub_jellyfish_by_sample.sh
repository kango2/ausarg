# ONT samples
for sampleID in $(cat matching-names.txt); do for klength in 9 11 13 17; do OUTDIR="/g/data/te53/nj8315/jellyfish_kmerlength${klength}"; mkdir -p $OUTDIR; qsub -l storage=gdata/if89+gdata/te53 -o /g/data/te53/nj8315/logs -P te53 -v inputfiles=/g/data/te53/rdna/converted_fastq/${sampleID}.fastq.gz,OUTDIR=${OUTDIR},sampleID=${sampleID},klength=${klength} /g/data/te53/nj8315/jellyfish.sh; done; done

# Illumina pacbio samples
for sampleID in $(cat list.txt); do for klength in 9 11 13 17; do OUTDIR="/g/data/te53/nj8315/illumina_k${klength}"; mkdir -p $OUTDIR; qsub -l storage=gdata/if89+gdata/te53 -o /g/data/te53/nj8315/logs -P te53 -v inputfiles=/g/data/te53/rdna/converted-illumin-fastq/${sampleID}.fastq.gz,OUTDIR=${OUTDIR},sampleID=${sampleID},klength=${klength} /g/data/te53/nj8315/jellyfish.sh; done; done
