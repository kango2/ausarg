#!/bin/bash
#PBS -N alignment
#PBS -P xl04
#PBS -o /g/data/xl04/ka6418/bassiana/all_assemblies/depth
#PBS -e /g/data/xl04/ka6418/bassiana/all_assemblies/depth
#PBS -q normal
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l walltime=10:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l wd

module load minimap2 samtools

platform=$platform
rawreads=$rawreads
reference=$reference
output=$output

# Step 1: Index the reference genome (minimap2 doesn't require explicit indexing)

if [ "$platform" == "illumina" ]; then
    # Splitting paired-end read sets by semicolon
    IFS=';' read -ra PAIRS <<< "$rawreads"
    for pair in "${PAIRS[@]}"; do
        IFS=':' read -ra ADDR <<< "$pair"
        R1="${ADDR[0]}"
        R2="${ADDR[1]}"
        minimap2 -t ${PBS_NCPUS} -a "$reference" "$R1" "$R2" >> "${output}.sam"
    done
elif [ "$platform" == "ont" ]; then
    IFS=':' read -ra ADDR <<< "$rawreads"
    for read_file in "${ADDR[@]}"; do
        minimap2 -t ${PBS_NCPUS} -ax map-ont "$reference" "$read_file" >> "${output}.sam"
    done
elif [ "$platform" == "pacbio" ]; then
    IFS=':' read -ra ADDR <<< "$rawreads"
    for read_file in "${ADDR[@]}"; do
        minimap2 -t ${PBS_NCPUS} -ax map-pb "$reference" "$read_file" >> "${output}.sam"
    done
else
    echo "Invalid platform. Choose from 'illumina', 'ont', or 'pacbio'."
    exit 1
fi

# Step 2: Convert SAM to BAM
samtools view -bS "${output}.sam" > "${output}.bam"

# Step 3: Sort and index BAM file
samtools sort -@ ${PBS_NCPUS} "${output}.bam" -o "${output}_sorted.bam"
samtools index "${output}_sorted.bam"

# Clean up intermediate files
rm "${output}.sam" "${output}.bam"

samtools depth -a "${output}_sorted.bam" | \
awk '{
    total+=$3;
    if (NR % 1000 == 0) {
        print $1 "," NR-999 "," $2 "," total/1000;
        total=0;
    }
}
END {
    if (NR % 1000 != 0) {
        print $1 "," NR-(NR%1000)+1 "," $2 "," total/(NR%1000);
    }
}' | \
awk -F, 'BEGIN {OFS=FS; prev=""; count=1} NR==1 {print $0; next} $1 != prev {count=1; prev=$1} {print $1, count, count+999, $4; count+=1000}' > "${output}_sorted.bam.binned.depth.csv"
