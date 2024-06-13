import sys
import os
import argparse
import pysam
import pandas as pd
import matplotlib.pyplot as plt

def read_fasta_file(fasta_file):
    sequence_dict = {}
    ref_lengths = {}
    with pysam.FastaFile(fasta_file) as fasta:
        for seq in fasta.references:
            sequence_dict[seq] = fasta.fetch(seq)
            ref_lengths[seq] = fasta.get_reference_length(seq)
    return sequence_dict, ref_lengths

def read_alignment_file(filename, ref_lengths):
    alnInfo = []  # List to store data for DataFrame
    unique_sequences = {} # identify unique read sequences that match alignment criteria

    with pysam.AlignmentFile(filename, "r") as file:
        for alignment in file:
            if not alignment.is_unmapped and alignment.query_alignment_length >= 500:
                alignment_score = alignment.get_tag('AS') if alignment.has_tag('AS') else 0
                if alignment_score >= 300:
                    ref_length = ref_lengths.get(alignment.reference_name, 0)
                    len_normalised_aln_score = alignment_score / ref_length if ref_length else 0
                    refid, refType = alignment.reference_name.split('_') if '_' in alignment.reference_name else (alignment.reference_name, '')
                    if alignment.query_name not in unique_sequences:
                        unique_sequences[alignment.query_name] = alignment.query_sequence

                    alnInfo.append({
                        "query_name": alignment.query_name,
                        "reference_name": alignment.reference_name,
                        "ref_length": ref_length,
                        "alignment_score": alignment_score,
                        "len_normalised_aln_score": len_normalised_aln_score,
                        "refid": refid,
                        "refType": refType
                    })
    return pd.DataFrame(alnInfo), unique_sequences

def process_bam_files(bam_files, ref_lengths, seqtech, output_prefix):
    combined_alnInfo = pd.DataFrame()

    for bam_file in bam_files:
        print(f"Reading {bam_file}")
        alnInfo, query_sequences = read_alignment_file(bam_file, ref_lengths)
        alnInfo['filename'] = bam_file
        combined_alnInfo = pd.concat([combined_alnInfo, alnInfo])

        with open(f'{output_prefix}_unique_sequences.fasta', 'a+') as file:
            for query_name, sequence in query_sequences.items():
                file.write(f">{query_name}\n{sequence}\n")

    combined_alnInfo['seqtech'] = seqtech
    return combined_alnInfo

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process BAM files and generate alignment data.")
    parser.add_argument('-p', '--pacbiobam', nargs='+', help='BAM files generated using PacBio data')
    parser.add_argument('-o', '--ontbam', nargs='+', help='BAM files generated using ONT data')
    parser.add_argument('-r', '--reference', required=True, help='Reference FASTA file')
    parser.add_argument('-d', '--output_dir', default='.', help='Output directory for the files')

    args = parser.parse_args()

    if not args.pacbiobam and not args.ontbam:
        parser.error('At least one of --pacbiobam or --ontbam must be provided.')

    if not os.path.isfile(args.reference):
        parser.error(f"Reference file '{args.reference}' does not exist.")

    if os.path.exists(args.output_dir):
        print(f"Output directory '{args.output_dir}' already exists. Please provide a non-existing directory.")
        sys.exit(1)
    else:
        os.makedirs(args.output_dir)

    return args

def generate_output(all_data, sequence_dict, output_dir):

    max_scores = all_data.groupby(['refType', 'seqtech'])['alignment_score'].max().reset_index().rename(columns={'alignment_score': 'max_score'})
    all_data = all_data.merge(max_scores, on=['refType', 'seqtech'])
    all_data['tech_normalised_score'] = all_data['alignment_score'] / all_data['max_score']
    all_data.to_csv(f'{output_dir}/reads2rdna_alnInfo.txt', index=False, sep ='\t')

    summed_scores = all_data.groupby(['refType', 'reference_name'])['tech_normalised_score'].sum().reset_index()
    best_matches = summed_scores.loc[summed_scores.groupby('refType')['tech_normalised_score'].idxmax()]

    with open(f'{output_dir}/best_rDNA_seqs.fa', 'w') as file:
        for _, row in best_matches.iterrows():
            ref_type = row['refType']
            ref_name = row['reference_name']
            if ref_name in sequence_dict:
                sequence = sequence_dict[ref_name]
                file.write(f">{ref_name}\n{sequence}\n")
            else:
                print(f"Sequence not found for {ref_name}")

def main():
    args = parse_arguments()
    sequence_dict, ref_lengths = read_fasta_file(args.reference)

    all_data = pd.DataFrame()
    if args.pacbiobam:
        pacbio_data = process_bam_files(args.pacbiobam, ref_lengths, "pacbio", f'{args.output_dir}/pacbio')
        all_data = pd.concat([all_data, pacbio_data])

    if args.ontbam:
        ont_data = process_bam_files(args.ontbam, ref_lengths, "ont", f'{args.output_dir}/ont')
        all_data = pd.concat([all_data, ont_data])

    generate_output(all_data, sequence_dict, args.output_dir)

if __name__ == "__main__":
    main()
