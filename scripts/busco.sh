#PBS -N BUSCO
#PBS -P xl04
#PBS -q normalsr
#PBS -l walltime=3:00:00
#PBS -l mem=512GB
#PBS -l ncpus=104
#PBS -j oe
#PBS -l storage=gdata/xl04+gdata/if89+gdata/te53
#PBS -l wd

# Inputs:
#   - fasta: The path to the input fasta file (can be .fa, .fasta, .fa.gz, or .fasta.gz)
#   - outdir: The directory where the BUSCO output will be stored

usage() {
	echo "Usage: qsub -v fasta=/path/to/genome.fasta,outdir=/path/to/store/busco/results ./busco.sh" >&2
	echo
	exit 1
}

[ -z "${fasta}" ] && usage
[ -z "${outdir}" ] && usage


die() {
	echo "$1" >&2
	echo
	exit 1
}

set -ex 
module load singularity

lineage=/g/data/if89/datalib/busco
img=/g/data/if89/singularityimg/busco-5.4.7.sif

prefix=$(basename "${fasta}" | sed -E 's/\.(fa|fasta)(\.gz)?$//')

echo "Running BUSCO on ${fasta}, storing output in ${outdir}"

singularity exec ${img} busco --out_path ${outdir} \
 -o run_$prefix --offline -i ${fasta} -l sauropsida_odb10 \
 --download_path ${lineage} --cpu ${PBS_NCPUS} -m genome --tar -f || die "BUSCO failed"