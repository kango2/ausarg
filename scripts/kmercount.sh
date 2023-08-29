#!/bin/bash

#PBS -l ncpus=48
#PBS -l walltime=48:00:00
#PBS -l mem=190GB
#PBS -N kmercount
#PBS -j oe
#PBS -l jobfs=400GB

#Authors : Hardip Patel, Kirat Alreja, Arthur Georges

set -ex

usage() {
	echo "Usage: qsub -l storage=gdata/if89+gdata/projectcode -o /path/to/stdouterr -P projectcode -v inputfiles=/path/to/input1.fq.gz:/path/to/input2.fq.gz,OUTDIR=/path/to/output/directory,sampleID=myfavsample,klength=17 ./kmercount.sh" >&2
	echo
	exit 1
}

#check if input files are provided
[ -z "${inputfiles}" ] && usage
#output directory
[ -z "${OUTDIR}" ] && usage
#sample id
[ -z "${sampleID}" ] && usage
#kmer length
[ -z "${klength}" ] && usage

module load jellyfish/2.3.0 utils/0.0

IFS=":" read -ra filelist <<< "$inputfiles"

filecount=${#filelist[@]}

for file in "${filelist[@]}"; do
	#tried the following to increase efficiency but monitoring suggests that pigz uses <100% cpus indicating one thread only
	#TODO: use 3 or 4 threads and see if the CPU usage increases
  echo "pigz -c -d -p 2 ${file}"  >> ${PBS_JOBFS}/filegenerator.cmds
done

##jellyfish --text option will output text strings directly. So no need to create counts database and then dump
##dump files are about 80-90GB when quality filters are used for 1KGP data
##sort uses 150GB of virtual memory, temporary files in text format are written to disk that are about 19GB in size each
##TODO: dump files can be compressed after sort to save disk space
##specifying -s 32G reaches up to 142.5GB of RAM usage during the count step
##TODO: incorporate other jellyfish parameters for quality filter if required
##TODO: perhaps combine all .done files into one file sensibly

runcmd.sh -c "jellyfish count --text -g ${PBS_JOBFS}/filegenerator.cmds -G ${filecount} --threads=$((PBS_NCPUS - filecount * 2)) -m ${klength} -o ${PBS_JOBFS}/$sampleID.dump -C -s 32G" -t ${sampleID}.kmercount -d ${OUTDIR}/$sampleID.kmercount.done -f false
runcmd.sh -c "jellyfish histo -t ${PBS_NCPUS} -o ${OUTDIR}/$sampleID.histo ${PBS_JOBFS}/$sampleID.dump" -t ${sampleID}.kmerhisto -d ${OUTDIR}/$sampleID.kmerhisto.done -f false
runcmd.sh -c "sort -k1,1 --parallel=${PBS_NCPUS} --buffer-size=80% -T ${PBS_JOBFS} -o ${PBS_JOBFS}/$sampleID.sorted.dump ${PBS_JOBFS}/$sampleID.dump" -t ${sampleID}.kmersort -d ${OUTDIR}/$sampleID.kmersort.done -f false
runcmd.sh -c "pigz -q -p ${PBS_NCPUS} ${PBS_JOBFS}/$sampleID.sorted.dump" -t ${sampleID}.kmerzip -d ${OUTDIR}/$sampleID.kmerzip.done -f false
runcmd.sh -c "rsync -a ${PBS_JOBFS}/$sampleID.sorted.dump.gz ${OUTDIR}/" -t ${sampleID}.kmerxfer -d ${OUTDIR}/$sampleID.kmerxfer.done -f false                              