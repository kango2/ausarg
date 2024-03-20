#!/bin/bash
#PBS -N paftosyntent
#PBS -P xl04
#PBS -q express
#PBS -l walltime=1:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

#INPUTS 

#PAF alignment ${paf}
#RefName ${refname}
#TargetName ${targetname}
#RefFASTA ${ref}
#TargetFASTA ${target}
#minlen of chromosomes to include ${minlen}
#OutDir ${outdit}

#OUTPUT 

#chromsyn.pdf 

module load samtools Rlib

cd ${outdir}

/g/data/xl04/ka6418/temp/buttonpaf/pafconvert.sh ${paf} ${refname} ${targetname} "${outdir}/regdata.tsv"

samtools faidx ${ref}
samtools faidx ${target}


FAIDX_FILE="${ref}.fai"
OUTPUT_FILE="${outdir}/${refname}.tdt"

# Write the header to the output file
echo -e "SeqName\tSeqLen\tTel5\tTel3\tTel5Len\tTel3Len\tTrim5\tTrim3\tTelPerc" > "$OUTPUT_FILE"

# Process the FAIDX file
while IFS=$'\t' read -r seqName seqLen rest; do
    # Ignore lines that do not conform to expected format
    if [[ ! $seqName || ! $seqLen ]]; then
        continue
    fi
    # Print the required format with correct values
    echo -e "${seqName}\t${seqLen}\tFalse\tFalse\t0\t0\t0\t0\t0.0" >> "$OUTPUT_FILE"
done < "$FAIDX_FILE"

echo "${refname} $OUTPUT_FILE" > "${outdir}/sequences.fofn"

FAIDX_FILE="${target}.fai"
OUTPUT_FILE="${outdir}/${targetname}.tdt"

# Write the header to the output file
echo -e "SeqName\tSeqLen\tTel5\tTel3\tTel5Len\tTel3Len\tTrim5\tTrim3\tTelPerc" > "$OUTPUT_FILE"

# Process the FAIDX file
while IFS=$'\t' read -r seqName seqLen rest; do
    # Ignore lines that do not conform to expected format
    if [[ ! $seqName || ! $seqLen ]]; then
        continue
    fi
    # Print the required format with correct values
    echo -e "${seqName}\t${seqLen}\tFalse\tFalse\t0\t0\t0\t0\t0.0" >> "$OUTPUT_FILE"
done < "$FAIDX_FILE"

echo "${targetname} $OUTPUT_FILE" >> "${outdir}/sequences.fofn"

Rscript /g/data/xl04/ka6418/temp/chromsyn_testrun/github/chromsyn/chromsyn.R sequences="${outdir}/sequences.fofn" regdata="${outdir}/regdata.tsv" minlen=${minlen} regmirror=True












