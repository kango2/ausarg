#!/bin/bash
#PBS -N alignment
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

module load minimap2 samtools

rawreads="$raw_reads"
reference="$reference_fasta"
output="$outputname"
technology="$tech"  # PacBio, ONT, or Illumina

# Step 1: Align reads to the reference genome based on technology
IFS=':' read -ra ADDR <<< "$rawreads"
for ((i=0; i<${#ADDR[@]}; i++)); do
    read_file="${ADDR[$i]}"
    case $technology in
        PacBio)
            minimap2 -t 48 -ax map-pb "$reference" "$read_file" >> "$outputname.sam"
            ;;
        ONT)
            minimap2 -t 48 -ax map-ont "$reference" "$read_file" >> "$outputname.sam"
            ;;
        Illumina)
            # Assuming paired-end reads: process read_file and the next file in the ADDR array as its R2
            if [[ $((i+1)) -lt ${#ADDR[@]} ]]; then
                read_file2="${ADDR[$((i+1))]}"
                minimap2 -t 48 -ax sr "$reference" "$read_file" "$read_file2" >> "$outputname.sam"
                # Increment the counter to skip the next file since it's already processed as R2
                ((i++))
            else
                echo "Mismatch in paired-end files. Exiting."
                exit 1
            fi
            ;;
        *)
            echo "Unknown sequencing technology. Exiting."
            exit 1
            ;;
    esac
done

# Step 2: Convert SAM to BAM
samtools view -bS "$outputname.sam" > "$outputname.bam"

# Step 3: Sort and index BAM file
samtools sort "$outputname.bam" -o "${outputname}_sorted.bam"
samtools index "${outputname}_sorted.bam"

# Clean up intermediate files
rm "$outputname.sam" "$outputname.bam"

# Step 4: Generate binned depth over 1000 base intervals in CSV format
samtools depth -a "${outputname}_sorted.bam" | awk 'BEGIN {OFS=","; print "Chromosome","Start","End","AverageDepth"} {
    total+=$3;
    if (NR % 1000 == 0) {
        print $1, NR-999, $2, total/1000;
        total=0;
    }
}
END {
    if (NR % 1000 != 0) {
        print $1, NR-(NR%1000)+1, $2, total/(NR%1000);
    }
}' > "${outputname}_sorted.bam.binned.depth.csv"