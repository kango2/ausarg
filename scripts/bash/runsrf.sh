#!/bin/bash
#PBS -l ncpus=48,mem=190GB,walltime=01:00:00,storage=gdata/if89+gdata/te53+gdata/xl04,jobfs=400GB
#PBS -q normal
#PBS -P xl04
#PBS -N srf
#PBS -j oe
#PBS -o /g/data/xl04/genomeprojects/Pogona_vitticeps/logs/

set -ex

#example usage: qsub -v inputfasta=asm.fasta,OUTPUTDIR=/path/to/output/directory/,window2refchains=window2ref.chains runtrf.sh
#todo: may not be issues in the longer run once workflow is settled
#1. output is only taking care of .fasta.gz extension, need to work with variety of extensions
#2. assumes bgzip compressed file for samtools faidx

module load KMC/3.2.4 k8/1.0 minimap2/2.28 parallel/20191022 samtools/1.19.2 pythonlib/3.9.2 TRF-mod/4.10.0 utils/0.0

##need to put this in if89 or other central location later
export SRFBIN=/g/data/te53/hrp561/tmp/srf/
#export MASKFA=/g/data/te53/hrp561/refgen/python/maskFasta.py

mkdir -p ${OUTPUTDIR}
cd ${PBS_JOBFS}

# Check if the process has already been completed
# if [ -e ${OUTPUTDIR}/${sampleid}.srftrf.percontig.done ]
# then
#     echo "Skipping ${sampleid} as analysis is already done"
#     exit 0
# fi

# # Run the srf analysis
# ## run srf on each contig. This will avoid collapsing alpha-HORs that may be "identical" between contigs/chromosomes
# ## a pooled version for all contigs can also be done if needed to identify alpha-HORs that are shared between individuals and contigs/chromosomes. this is a later problem
# ## create index if required

# if [ ! -e ${inputfasta}.fai ]
# then
#     samtools faidx ${inputfasta}
# fi

# # Loop through each contig in the input fasta index file
# for contig in $(cut -f1 ${inputfasta}.fai)
# do
#     # Check if the analysis for this contig is already done
#     if [ -s ${sampleid}.${contig}.srf.done ]
#     then
#         echo "Skipping ${sampleid}.${contig} as analysis is already done"
#         continue
#     fi
    
#     # Extract the contig sequence from the input fasta file
#     samtools faidx ${inputfasta} ${contig} >${sampleid}.${contig}.fa
    
#     # Run KMC to count k-mers in the contig sequence
#     kmc -fm -k${klen} -t${PBS_NCPUS} -ci20 -cs100000 ${sampleid}.${contig}.fa ${sampleid}.${contig}.kmc ./ &>${sampleid}.${contig}.kmc.log
#     kmc_dump ${sampleid}.${contig}.kmc ${sampleid}.${contig}.counts.txt &>>${sampleid}.${contig}.kmc.log
    
#     # Check if k-mer counting was successful
#     if [ -s ${sampleid}.${contig}.counts.txt ]
#     then
#         # Run SRF analysis on the contig
#         ${SRFBIN}/srf -p ${sampleid}.${contig} ${sampleid}.${contig}.counts.txt >${sampleid}.${contig}.srf.fa
        
#         # Align the SRF sequences back to the contig
#         minimap2 -t ${PBS_NCPUS} -c -N1000000 -f1000 -r100,100 <(${SRFBIN}/srfutils.js enlong ${sampleid}.${contig}.srf.fa) ${sampleid}.${contig}.fa >${sampleid}.${contig}.asm2srf.paf
        
#         # Convert the alignment to BED format
#         ${SRFBIN}/srfutils.js paf2bed ${sampleid}.${contig}.asm2srf.paf > ${sampleid}.${contig}.asm2srf.bed
        
#         # Generate abundance information from the BED file
#         ${SRFBIN}/srfutils.js bed2abun ${sampleid}.${contig}.asm2srf.bed > ${sampleid}.${contig}.asm2srf.abun
        
#         # Mask the contig sequence based on the BED file
#         maskFasta.py --fasta_file ${sampleid}.${contig}.fa --bed_file ${sampleid}.${contig}.asm2srf.bed --output_file ${sampleid}.${contig}.masked.fa
#     else
#         continue
#     fi
#     # Mark the contig as done
#     touch ${sampleid}.${contig}.srf.done
# done

# cat ${sampleid}.*.asm2srf.abun >${OUTPUTDIR}/${sampleid}.asm2srf.percontig.abun
# cat ${sampleid}.*.asm2srf.bed >${OUTPUTDIR}/${sampleid}.asm2srf.percontig.bed
# cat ${sampleid}.*.asm2srf.paf >${OUTPUTDIR}/${sampleid}.asm2srf.percontig.paf
# cat ${sampleid}.*.srf.fa >${OUTPUTDIR}/${sampleid}.srf.percontig.fa
# cat ${sampleid}.*.kmc.log >${OUTPUTDIR}/${sampleid}.kmc.percontig.log

# Check if the analysis for this contig is already done
if [ -s ${sampleid}.srf.genome.done ]
then
    echo "Skipping ${sampleid} as analysis is already done"
    exit 0
fi

##Run srf for the whole genome
# Run KMC to count k-mers in the genome
kmc -fm -k${klen} -t${PBS_NCPUS} -ci20 -cs100000 ${inputfasta} ${sampleid}.kmc.genome ./ &>${sampleid}.kmc.genome.log
kmc_dump ${sampleid}.kmc.genome ${sampleid}.counts.genome.txt &>>${sampleid}.kmc.genome.log

# Run SRF analysis on the genome
${SRFBIN}/srf -p ${sampleid} ${sampleid}.counts.genome.txt >${sampleid}.srf.genome.fa

# Align the SRF sequences back to the genome
minimap2 -t ${PBS_NCPUS} -c -N1000000 -f1000 -r100,100 <(${SRFBIN}/srfutils.js enlong ${sampleid}.srf.genome.fa) ${inputfasta} >${sampleid}.asm2srf.genome.paf

# Convert the alignment to BED format
${SRFBIN}/srfutils.js paf2bed ${sampleid}.asm2srf.genome.paf > ${sampleid}.asm2srf.genome.bed

# Generate abundance information from the BED file
${SRFBIN}/srfutils.js bed2abun ${sampleid}.asm2srf.genome.bed > ${sampleid}.asm2srf.genome.abun

# Copy files to ${OUTPUTDIR}
cp ${sampleid}.srf.genome.fa ${OUTPUTDIR}/
cp ${sampleid}.asm2srf.genome.paf ${OUTPUTDIR}/
cp ${sampleid}.asm2srf.genome.bed ${OUTPUTDIR}/
cp ${sampleid}.asm2srf.genome.abun ${OUTPUTDIR}/
cp ${sampleid}.kmc.genome.log ${OUTPUTDIR}/

# Mark the process as done
touch ${OUTPUTDIR}/${sampleid}.srf.genome.done

## following is optional since PBS_JOBFS will be cleaned up by PBS
## this is for insurance purposes where PBS_JOBFS is not used
for contig in $(cut -f1 ${inputfasta}.fai)
do
    rm -f ${sampleid}.${contig}.fa
    rm -f ${sampleid}.${contig}.kmc
    rm -f ${sampleid}.${contig}.counts.txt
    rm -f ${sampleid}.${contig}.srf.fa
    rm -f ${sampleid}.${contig}.asm2srf.paf
    rm -f ${sampleid}.${contig}.asm2srf.bed
    rm -f ${sampleid}.${contig}.asm2srf.abun
    rm -f ${sampleid}.${contig}.masked.fa
    rm -f ${sampleid}.${contig}.kmc.log
    rm -f ${sampleid}.${contig}.kmc.kmc_suf
    rm -f ${sampleid}.${contig}.kmc.kmc_pre
#    rm -f ${sampleid}.${contig}.masked.trf
done


