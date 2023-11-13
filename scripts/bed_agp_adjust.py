import sys

def read_agp_file(agp_file):
    agp_data = {}
    with open(agp_file, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if parts[4] == 'W':  # We only care about lines representing contigs
                scaffold = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                contig = parts[5]
                orientation = parts[8]
                agp_data[contig] = (scaffold, start, end, orientation)
    return agp_data

def adjust_bed_file(bed_file, agp_data, output_file):
    with open(bed_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('track') or line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            contig, start, end = parts[0], int(parts[1]), int(parts[2])
            if contig in agp_data:
                scaffold, scaffold_start, scaffold_end, orientation = agp_data[contig]
                if orientation == '+':
                    new_start = scaffold_start + start - 1
                    new_end = scaffold_start + end - 1
                else:
                    new_start = scaffold_end - end + 1
                    new_end = scaffold_end - start + 1
                outfile.write(f'{scaffold}\t{new_start}\t{new_end}\t' + "\t".join(parts[3:]) + '\n')

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input.bed> <input.agp>")
        sys.exit(1)

    bed_file = sys.argv[1]
    agp_file = sys.argv[2]
    output_file = "adjusted_annotations.bed"

    agp_data = read_agp_file(agp_file)
    adjust_bed_file(bed_file, agp_data, output_file)
    print(f"Adjusted annotations have been written to {output_file}")

if __name__ == "__main__":
    main()
