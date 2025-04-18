#!/bin/bash
#PBS -P xl04
#PBS -N EEL
#PBS -q gpuvolta
#PBS -l ncpus=48
#PBS -l ngpus=4
#PBS -l mem=384GB
#PBS -l walltime=0:20:00
#PBS -l wd
#PBS -l storage=gdata/if89+scratch/xl04+gdata/xl04

###################################################################

# Change this to the model you want to use
#MODEL=dna_r10.4.1_e8.2_400bps_5khz_sup.cfg
# MODEL=dna_r10.4.1_e8.2_400bps_sup.cfg
# MODEL=dna_r10.4.1_e8.2_400bps_hac_prom.cfg
#MODEL=dna_r9.4.1_450bps_sup.cfg
MODEL=dna_r9.4.1_450bps_hac_prom.cfg

###################################################################

# Make sure to change:
# 1. wv19 to your own project
# 2. the name of the Guppy model
# 3. optionally, if you want to use the A100 GPU queue instead of the V100 queue, change "gpuvolta" to "dgxa100" and change the number of CPUs to 64 (dgxa100 requires at least 16 CPUs per GPU)

# to run:
# qsub -v MERGED_SLOW5=/path/to/reads.blow5,BASECALL_OUT=/path/to/out/dir ./buttery-eel.pbs.sh

###################################################################

usage() {
	echo "Usage: qsub -v MERGED_SLOW5=/g/data/wv19/public/hg2_prom_lsk114_subsubsample/reads.blow5,BASECALL_OUT=/scratch/wv19/hg1112/tmp/hg2_prom_lsk114_subsubsample/ ./buttery-eel.pbs.sh" >&2
	echo
	exit 1
}

#directory where basecalls should be written to
[ -z "${BASECALL_OUT}" ] && usage
# merged BLOW5
[ -z "${MERGED_SLOW5}" ] && usage

module load /g/data/if89/apps/modulefiles/buttery-eel/0.4.2+dorado7.2.13

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

test -d ${BASECALL_OUT} && die "Output directory ${BASECALL_OUT} already exists. Please delete it first or give an alternate location. Exiting."

test -e ${MERGED_SLOW5} || die "${MERGED_SLOW5} not found. Exiting."

mkdir ${BASECALL_OUT} || die "Creating directory ${BASECALL_OUT} failed. Exiting."
cd ${BASECALL_OUT} || die "${MERGED_SLOW5} not found. Exiting."

/usr/bin/time -v  buttery-eel -i ${MERGED_SLOW5} -o ${BASECALL_OUT}/reads.fastq --guppy_bin ${ONT_DORADO_PATH} --port ${PORT} --use_tcp --config ${MODEL} -x cuda:all --guppy_batchsize 20000 --max_queued_reads 20000 --slow5_threads 10 --slow5_batchsize 100 --procs 20 --detect_mid_strand_adapter --trim_adapters --detect_adapter --do_read_splitting --qscore 7 || die "basecalling failed"

echo "basecalling success"