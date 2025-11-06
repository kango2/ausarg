#PBS -N BUSCO
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=3:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l jobfs=400GB
#PBS -j oe
#PBS -l storage=gdata/xl04+gdata/if89+gdata/te53
#PBS -l wd

# Inputs:
#   - fasta: The path to the input fasta file (can be .fa, .fasta, .fa.gz, or .fasta.gz)
#   - outdir: The directory where the BUSCO output will be stored
#   - prefix: prefix for busco result (folder)

die() {
	echo "$1" >&2
	echo
	exit 1
}

set -ex 
module load singularity
module load ezlabgva/busco/v5.8.2_cv1

lineage=/g/data/if89/datalib/busco

echo "Running BUSCO on ${fasta}, storing output in ${outdir}"



# Check if the FASTA is gzipped
if [[ "${fasta}" == *.gz ]]; then
    echo "FASTA file is gzipped. Decompressing to jobfs with original base name."
	fasta_base=$(basename "${fasta}" .gz)
    fasta_to_use="${PBS_JOBFS}/${fasta_base}"
    pigz -dc "${fasta}" > "${fasta_to_use}" || die "Failed to decompress FASTA file with pigz."
else
    fasta_to_use=${fasta}
fi


busco --mode proteins --out_path ${outdir} -o ${folder} --offline -i ${fasta_to_use} -l sauropsida_odb10 --download_path ${lineage} --cpu ${PBS_NCPUS} --tar -f || die "BUSCO failed"