#!/bin/bash
#PBS -P xl04
#PBS -N EEL
#PBS -q gpuvolta
#PBS -l ncpus=48
#PBS -l ngpus=4
#PBS -l mem=384GB
#PBS -l walltime=48:00:00
#PBS -l wd
#PBS -l storage=gdata/if89+scratch/xl04+gdata/xl04

set -ex

###################################################################

# Change this to the model you want to use
#MODEL=dna_r10.4.1_e8.2_400bps_5khz_sup.cfg
# MODEL=dna_r10.4.1_e8.2_400bps_sup.cfg
# MODEL=dna_r10.4.1_e8.2_400bps_hac_prom.cfg
#MODEL=dna_r9.4.1_450bps_sup.cfg
#MODEL=dna_r9.4.1_450bps_hac_prom.cfg

###################################################################

# Make sure to change:
# 1. wv19 to your own project
# 2. the name of the Guppy model
# 3. optionally, if you want to use the A100 GPU queue instead of the V100 queue, change "gpuvolta" to "dgxa100" and change the number of CPUs to 64 (dgxa100 requires at least 16 CPUs per GPU)

# to run:
# qsub -v MERGED_SLOW5=/path/to/reads.blow5,BASECALL_OUT=/path/to/out/dir ./buttery-eel.pbs.sh

###################################################################

#directory where basecalls should be written to
[ -z "${BASECALL_OUT}" ] && usage
# merged BLOW5
[ -z "${MERGED_SLOW5}" ] && usage

module load /g/data/if89/apps/modulefiles/buttery-eel/0.5.1+dorado7.4.12 htslib

###################################################################

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

#https://unix.stackexchange.com/questions/55913/whats-the-easiest-way-to-find-an-unused-local-port
PORT=5000
get_free_port() {
	for port in $(seq 5000 65000); do
		echo "trying port $port" >&2
		PORT=$port
		ss -lpna | grep -q ":$port " || break
	done
}

get_free_port
test -z "${PORT}" && die "Could not find a free port"
echo "Using port ${PORT}"

ONT_DORADO_PATH=$(which dorado_basecall_server | sed "s/dorado\_basecall\_server$//")/
${ONT_DORADO_PATH}/dorado_basecall_server --version || die "Could not find dorado_basecall_server"

test -e ${MERGED_SLOW5} || die "${MERGED_SLOW5} not found. Exiting."

cd ${BASECALL_OUT} || die "${MERGED_SLOW5} not found. Exiting."

buttery-eel -i ${MERGED_SLOW5} -o "${BASECALL_OUT}/$(basename ${MERGED_SLOW5} .blow5).${MODE}.fastq" -g ${ONT_DORADO_PATH} --port ${PORT} --use_tcp --config ${MODEL} -x cuda:all --slow5_threads 10 --slow5_batchsize 4000 --procs 20 --trim_adapters --qscore 7 || die "basecalling failed"

md5sum "$(basename ${MERGED_SLOW5} .blow5).${MODE}.pass.fastq" > $(basename ${MERGED_SLOW5} .blow5).${MODE}.fastq.md5
md5sum "$(basename ${MERGED_SLOW5} .blow5).${MODE}.fail.fastq" >> $(basename ${MERGED_SLOW5} .blow5 ).${MODE}.fastq.md5

bgzip --index -@ ${PBS_NCPUS} "${BASECALL_OUT}/$(basename ${MERGED_SLOW5} .blow5).${MODE}.pass.fastq"
bgzip --index -@ ${PBS_NCPUS} "${BASECALL_OUT}/$(basename ${MERGED_SLOW5} .blow5).${MODE}.fail.fastq"

md5sum "$(basename ${MERGED_SLOW5} .blow5).${MODE}.pass.fastq.gz" >> $(basename ${MERGED_SLOW5} .blow5).${MODE}.fastq.md5
md5sum "$(basename ${MERGED_SLOW5} .blow5).${MODE}.pass.fastq.gz.gzi" >> $(basename ${MERGED_SLOW5} .blow5).${MODE}.fastq.md5
md5sum "$(basename ${MERGED_SLOW5} .blow5).${MODE}.fail.fastq.gz" >> $(basename ${MERGED_SLOW5} .blow5).${MODE}.fastq.md5
md5sum "$(basename ${MERGED_SLOW5} .blow5).${MODE}.fail.fastq.gz.gzi" >> $(basename ${MERGED_SLOW5}.blow5).${MODE}.fastq.md5

echo "basecalling success"