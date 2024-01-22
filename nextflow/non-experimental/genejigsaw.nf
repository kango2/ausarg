include { fromQuery } from 'plugin/nf-sqldb'
include { hifiasm_assembly } from '/g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/assembly_jigsaw.nf'
include { arima_mapping } from '/g/data/xl04/ka6418/github/ausarg/nextflow/scaffolding_pipeline/scaffold_jigsaw.nf'
include { yahs } from '/g/data/xl04/ka6418/github/ausarg/nextflow/scaffolding_pipeline/scaffold_jigsaw.nf'
include { generate_hicmap } from '/g/data/xl04/ka6418/github/ausarg/nextflow/scaffolding_pipeline/scaffold_jigsaw.nf'
include { setup_directory } from '/g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/setup_folders.nf'

process longread_qc {
    conda '/g/data/xl04/ka6418/miniconda/envs/genejigsaw'
    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=1,mem=8GB,storage=gdata/if89+gdata/xl04'

    input:
    tuple path (fastq_file),val (sample),val (flowcell),path (output),val (platform)

    output:
    tuple path ("$output/*_stats.csv"), path ("$output/*_quality_freq.csv") , path ("$output/*_length_freq.csv"), path ("$output/*_log.txt")


    script:

    """
#!/usr/bin/env python
from Bio import SeqIO
import subprocess
import sys
import os
import csv
import gzip
import argparse
import datetime

bin_size = 100

def calculate_n50_n90(read_lengths):
    sorted_lengths = sorted(read_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    cumulative_length = 0
    n50 = n90 = l50 = l90 = 0

    for i, length in enumerate(sorted_lengths):
        cumulative_length += length
        if not n50 and cumulative_length >= total_length * 0.5:
            n50 = length
            l50 = i + 1
        if not n90 and cumulative_length >= total_length * 0.9:
            n90 = length
            l90 = i + 1
            break

    return n50, n90, l50, l90

def log_progress(message, log_file, flowcell_id, input_file):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as log:
        log.write(f"{timestamp} - {message} - Flowcell: {flowcell_id} - File: {input_file}")


def process_fastq(input_fastq, output_path, flowcell_id, platform,sample, log_file):
    bins = {}
    length_sums = {}
    read_lengths = []
    total_bases = 0
    total_reads = 0
    total_ns = 0

    log_progress("Starting FASTQ file processing.", log_file,flowcell_id, input_fastq)

    with gzip.open(input_fastq, "rt") as f:
        for record in SeqIO.parse(f, "fastq-sanger"):
            sequence_length = len(record.seq)
            avg_qv = round(sum(record.letter_annotations["phred_quality"]) / sequence_length)
            total_bases += sequence_length
            total_reads += 1
            total_ns += record.seq.count("N")
            read_lengths.append(sequence_length)

            bin_number = (sequence_length) // bin_size * bin_size
            bin_key = (bin_number, avg_qv)

            if bin_key not in bins:
                bins[bin_key] = {"total_qv": 0, "count": 0}
            bins[bin_key]["total_qv"] += avg_qv
            bins[bin_key]["count"] += 1

            if bin_number not in length_sums:
                length_sums[bin_number] = 0
            length_sums[bin_number] += 1

    n50, n90, l50, l90 = calculate_n50_n90(read_lengths)
    average_read_length = total_bases / total_reads if total_reads > 0 else 0

    log_progress("FASTQ file processing completed. Writing CSV files...", log_file, flowcell_id, input_fastq)

    file_prefix = f"{sample}_{flowcell_id}_{platform}"
    quality_output_csv = os.path.join(output_path, f"{file_prefix}_quality_freq.csv")
    length_output_csv = os.path.join(output_path, f"{file_prefix}_length_freq.csv")
    stats_output_csv = os.path.join(output_path, f"{file_prefix}_stats.csv")

    # Write quality frequency data
    with open(quality_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["File_Path","Sample","Flowcell_ID", "Platform", "Read_Length", "QV", "Read_Numbers"])
        for bin_key, bin_data in bins.items():
            length_bin, qv_bin = bin_key
            frequency = bin_data["count"]
            csv_writer.writerow([input_fastq,sample, flowcell_id, platform, length_bin, qv_bin, frequency])

    # Write length frequency data
    with open(length_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["File_Path","Sample", "Flowcell_ID", "Platform", "Read_Length", "Summed_Read_Numbers"])
        for length_bin, count in length_sums.items():
            csv_writer.writerow([input_fastq,sample, flowcell_id, platform, length_bin, count])

    # Write stats data
    with open(stats_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["File_Path","Sample", "Flowcell_ID", "Platform", "Total_Bases", "Total_Reads", "Average_Read_Length", "N50", "N90", "L50", "L90", "Total_Ns"])
        csv_writer.writerow([input_fastq,sample, flowcell_id, platform, total_bases, total_reads, average_read_length, n50, n90, l50, l90, total_ns])

    log_progress("CSV files successfully written.", log_file, flowcell_id, input_fastq)
    print(f"CSV output and log file written to: {output_path}")

def main():
    file_prefix = f"${sample}_${flowcell}_${platform}"
    log_file = os.path.join('$output', f"{file_prefix}_processing_log.txt")
    process_fastq('$fastq_file', '$output', '$flowcell', '$platform', '$sample', log_file)

if __name__ == "__main__":
    main()

"""
}

process longread_plots {
    
    module 'Rlib'
    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=2,mem=8GB,storage=gdata/if89+gdata/xl04'

    input:
    val (qc_results)
    val (output)
    val (qc_files)

    script:
    """
    Rscript /g/data/xl04/ka6418/github/ausarg/nextflow/long_read_plots_jigsaw.R -t "$qc_results" -o "$output" -p || true
    Rscript /g/data/xl04/ka6418/github/ausarg/nextflow/long_read_plots_jigsaw.R -t "$qc_results" -o "$output" -u || true

    """


}

process shortread_qc {

    conda '/g/data/xl04/ka6418/miniconda/envs/genejigsaw'
    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=1,mem=8GB,storage=gdata/if89+gdata/xl04'

    input:
    tuple path (R1), path (R2),val (sample),val (flowcell),path (output),val (platform)

    output:
    path ("$output/*_ILLUMINA.csv")


    script:
    """

#!/usr/bin/env python
import csv
import gzip
from collections import defaultdict
import hashlib
import argparse
import os
import json
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

def phred_to_quality(phred_string):

    return [ord(char) - 33 for char in phred_string]

def get_file_size_gb(file_path):

    size_in_bytes = os.path.getsize(file_path)
    size_in_gb = size_in_bytes / (1024 ** 2)  # Corrected to return in GB
    return round(size_in_gb, 2)

def calculate_md5_compressed(fastq_file):
    hash_md5_compressed = hashlib.md5()
    with open(fastq_file, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5_compressed.update(chunk)
    return hash_md5_compressed.hexdigest()

def calculate_md5_uncompressed(fastq_file):
    hash_md5_uncompressed = hashlib.md5()
    with gzip.open(fastq_file, 'rt') as f:
        for line in f:
            hash_md5_uncompressed.update(line.encode())
    return hash_md5_uncompressed.hexdigest()

def parse_fastq(fastq_file):

    num_reads = 0
    total_bases = 0
    quality_scores_sum = defaultdict(int)
    quality_scores_count = defaultdict(int)
    total_bases_count = defaultdict(int)
    avg_qv_reads = [0] * 101
    total_gc_count = 0
    nucleotide_freq_data = defaultdict(lambda: defaultdict(int))
    overall_content_data = defaultdict(int)
    BASES_ORDER = ['A', 'C', 'G', 'T', 'N']
    bases = set()
    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            # Sequence lines are the second line in every set of 4 lines
            if i % 4 == 1:
                seq = line.strip()
                num_reads += 1
                total_bases += len(seq)
                total_gc_count += seq.count('G') + seq.count('C')
                for idx, base in enumerate(seq):
                    total_bases_count[idx] += 1
                    nucleotide_freq_data[idx][base] += 1
                    overall_content_data[base] += 1
                    bases.update(base)
            elif i % 4 == 3:
                qs = line.strip()
                scores = phred_to_quality(qs)
                avg_quality = sum(scores) / len(scores) if scores else 0
                bin_index = min(int(avg_quality), 100)
                avg_qv_reads[bin_index] += 1
                for idx, score in enumerate(scores):
                    quality_scores_sum[idx] += score
                    quality_scores_count[idx] += 1

    avg_quality_values = [round(quality_scores_sum[i] / quality_scores_count[i]) if quality_scores_count[i] != 0 else 0 for i in range(len(quality_scores_sum))]
    gc_content = round((total_gc_count / total_bases) * 100, 2) if total_bases else 0
    nucleotide_freq = [":".join(str(nucleotide_freq_data[idx].get(base, 0)) for base in BASES_ORDER) for idx in nucleotide_freq_data]
    overall_content = ":".join(str(overall_content_data.get(base, 0)) for base in BASES_ORDER)

    with ThreadPoolExecutor() as executor:
        future_md5_compressed = executor.submit(calculate_md5_compressed, fastq_file)
        future_md5_uncompressed = executor.submit(calculate_md5_uncompressed, fastq_file)
        md5_compressed = future_md5_compressed.result()
        md5_uncompressed = future_md5_uncompressed.result()

    return {
        'num_reads': num_reads,
        'total_bases': total_bases,
        'avg_read_length': total_bases / num_reads if num_reads > 0 else 0,
        'avg_quality_values': avg_quality_values,
        'gc_content': gc_content,
        'md5_compressed': md5_compressed,
        'md5_uncompressed': md5_uncompressed,
        'avg_qv_reads': avg_qv_reads,
        'file_size_gb': get_file_size_gb(fastq_file),
        'path': os.path.abspath(fastq_file),
        'nucleotide_freq': nucleotide_freq,
        'overall_content': overall_content
    }

def main():


    metrics = []
    fastq_files = ['$R1','$R2']

    num_cores_per_file = 2 // len(fastq_files)

    with ProcessPoolExecutor(max_workers=num_cores_per_file) as executor:
        for result in executor.map(parse_fastq, fastq_files):
            metrics.append(result)

    headers = [
        'File_path','Number_of_reads', 'Number_of_bases', 'Mean_read_length',
        'Mean_QV_at_read_position', 'Nucleotide_count_at_read_position',
        'Nucleotide_content', 'Mean_QV_per_read',
        'MD5_zipped', 'MD5_text', 'GC', 'File_size_in_MB'
    ]

    combined_metrics = []
    keys = [
        'path','num_reads', 'total_bases', 'avg_read_length', 'avg_quality_values',
        'nucleotide_freq', 'overall_content', 'avg_qv_reads', 'md5_compressed',
        'md5_uncompressed', 'gc_content', 'file_size_gb'
    ]

    for key in keys:
        if isinstance(metrics[0][key], list):
            combined_val = ",".join(str(val) for val in metrics[0][key])
            if len(metrics) > 1:
                combined_val += ";" + ",".join(str(val) for val in metrics[1][key])
        else:
            combined_val = f"{metrics[0][key]}"
            if len(metrics) > 1:
                combined_val += f";{metrics[1][key]}"
        combined_metrics.append(combined_val)

    metrics_dict = dict(zip(headers, combined_metrics))


    with open(os.path.join('$output', f"{'$sample'}_{'$flowcell'}_{'$platform'}.csv"), 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(headers)
        csvwriter.writerow(combined_metrics)



if __name__ == "__main__":
    main()

 """ 
}


process shortread_trimming {
    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '10h'
    clusterOptions = '-l ncpus=1,mem=8GB,storage=gdata/if89+gdata/xl04'

    input:
    tuple path (R1),path (R2),val (sample),val (flowcell),path (output),val (platform)

    output:
    tuple path ("$output/*_val_1.fq.gz"),path ("$output/*_val_2.fq.gz"),val (sample),val (flowcell),path (output),val (platform) 


    script:

    """
    source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
    conda activate genejigsaw
    trim_galore -o '${output}' --cores \${PBS_NCPUS} --paired '${R1}' '${R2}'

    """

}

process assembly {
    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=8,mem=64GB,storage=gdata/if89+gdata/xl04'

    input:
    tuple val (sample), path (pacbio)
    tuple val (sample), path (ont)
    tuple val (sample), path(R1), path(R2)
    val (output)

    output:
    val ("${output}/${sample}_HifiASM.fasta")



    script:
    """
    /g/data/xl04/ka6418/bassiana/hifiasm_bassiana/hifiasm/hifiasm -t \${PBS_NCPUS} -o "$output/$sample" $pacbio --ul $ont --h1 $R1 --h2 $R2 

    awk '/^S/{print ">"\$2;print \$3}' ${output}/${sample}*.p_ctg.gfa > ${output}/${sample}_HifiASM.fasta

    """
}

process kmers {

    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '3h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04,jobfs=400GB'

    input:
    tuple val (sample), path (files), val(klength), val(tech), val (output)

    output:
    val ("${output}/${sample}_${tech}_${klength}.histo")
    

    script:
    """
    joined_files=""

    for file in $files; do
        if [ -z "\$joined_files" ]; then
            joined_files="\$file"
        else
            joined_files="\${joined_files}:\${file}"
        fi
    done

    /g/data/xl04/ka6418/github/ausarg/nextflow/kmer_nf.sh -i \${joined_files} -s $sample -o $output -l $klength -t $tech 

    """


}

process kmer_plotting {
    
    module 'pythonlib'
    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=2,mem=8GB,storage=gdata/if89+gdata/xl04,jobfs=100GB'

    input:
    val (kmer_results)
    val (output)
    val (histo)

    output:
    val ("$output/kmer_plots.png")


    script:
    """
    #!/usr/bin/env python
    import os
    import matplotlib.pyplot as plt
    import seaborn as sns
    from collections import defaultdict

    # Load data from file
    def load_data(file_path):
        x_values, y_values = [], []
        with open(file_path, 'r') as f:
            for line in f:
                x, y = line.strip().split()
                x_values.append(int(x))
                y_values.append(int(y))
        return x_values, y_values

    # Extract info from filename
    def extract_info_from_filename(filename):
        base_name = os.path.basename(filename).replace('.histo', '')
        specie_name, sequencing_tech, kmer_size = base_name.split("_")[:3]
        kmer_size = int(kmer_size)  # Convert kmer size to integer
        return specie_name, sequencing_tech, kmer_size

    # Filter data based on x-axis limit
    def filter_data(x_data, y_data, limit=500):
        filtered_x = [x for x in x_data if x <= limit]
        filtered_y = y_data[:len(filtered_x)]
        return filtered_x, filtered_y
    def automatic_plot_combined(files,save_path=None):
        
        files = [os.path.join(directory_path, filename) for filename in os.listdir(directory_path) if filename.endswith('.histo')]
        file_info = [extract_info_from_filename(file) for file in files]

        # Group files by sequencing tech
        grouped_files = defaultdict(list)
        for file, (specie_name, tech, kmer_size) in zip(files, file_info):
            grouped_files[tech].append((file, specie_name, kmer_size))

        total_rows = len(grouped_files)
        max_subplots_per_row = max(len(tech_files) for tech_files in grouped_files.values())

        # Create a single figure with multiple subplots arranged in rows by sequencing tech
        fig, ax = plt.subplots(total_rows, max_subplots_per_row, figsize=(5 * max_subplots_per_row, 5 * total_rows))

        # To handle cases where there's only one sequencing tech
        if total_rows == 1:
            ax = [ax]

        row_idx = 0
        for tech, tech_files in grouped_files.items():
            tech_files = sorted(tech_files, key=lambda x: x[2])  # Sort by k-mer size

            for col_idx, (file, specie_name, kmer_size) in enumerate(tech_files):
                x_data, y_data = filter_data(*load_data(file), limit=200)  # Set default limit to 200
                ax[row_idx][col_idx].bar(x_data, y_data, width=1, align='center', color=sns.color_palette("flare")[col_idx])
                ax[row_idx][col_idx].set_title(f"{specie_name} ({tech}, k={kmer_size})")
                ax[row_idx][col_idx].set_xlabel("Kmer Coverage")
                ax[row_idx][col_idx].set_yscale('log')

            row_idx += 1

        # Common y-axis label
        fig.text(0.00, 0.5, 'Frequency', va='center', rotation='vertical')
        
        # Set the mega title
        mega_title = "Kmer Analysis"
        plt.suptitle(mega_title, fontsize=16)

        # Adjust layout and show plot
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust top spacing for mega title
        
        fig.set_facecolor("white")
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        plt.show()


    directory_path = "$kmer_results"

    save_path = "$output/kmer_plots.png"

    # For demonstration purposes using the previously loaded files
    automatic_plot_combined(directory_path,save_path)

    """


}


workflow {


    //Set up the directories in the given output folder
    setup_directory(params.topfolder)

    //Prepare long read files for data QC
    longread_fastqs_qc = channel
        .fromQuery('select title, flowcell, platform, filename from SRA where platform IN ("PACBIO_SMRT", "OXFORD_NANOPORE")', db: 'inputdb')
        .map { row ->
            def (title, flowcell, platform, filename) = row
            def output = "${params.topfolder}/rawdata/longread/qc" 
            return [filename, title, flowcell, output, platform]
        }.view()
    
    //Run long read QC
    longread_qcfiles = longread_qc(longread_fastqs_qc)
    
    //Error : the longread_qcfiles need to be flattneded and passed for qc, this automatic file path detection approach will not work, and keeps failing.

    //Plot the long read QC results
    //longread_plots("${params.topfolder}/rawdata/longread/qc","${params.topfolder}/rawdata/longread/qc",[longread_qcfiles[0]])

    //Prepare short read files for data QC
    illumina_fastqs_qc = channel
    .fromQuery('select title, flowcell, platform, filename from SRA where platform is "ILLUMINA"', db: 'inputdb')
    .map { row ->
        def (title, flowcell, platform, filename) = row
        def output = "${params.topfolder}/rawdata/shortread/qc"
        def (file1, file2) = filename.split(':')
        return [file1, file2, title, flowcell, output, platform]
    }.view()
    
    //Trim short reads and run QC
    trimmed_shortreads_qc = shortread_trimming(illumina_fastqs_qc)
    shortread_qc(trimmed_shortreads_qc)

    //Prepare channels for K-Mer analysis 
    //Analysis will run for as many kmer values as listed here
    def kmerValues = [21]

    pacbio_kmer = channel
        .fromQuery('select title, filename from SRA where platform is "PACBIO_SMRT"', db: 'inputdb')
        .map { row ->
            def (title, filename) = row
            def pacbio_file = file(filename)
            return [title, pacbio_file]
        }
        .groupTuple()
        .flatMap { title, filePairs ->
            def allFiles = filePairs.collect { it.toString().split(':') }.flatten()
            // Create a new entry for each kmerval
            kmerValues.collect { kmerval -> [title, allFiles, kmerval, "PacBio","${params.topfolder}/rawdata/kmers"] }
        }
        .view()
    
    ont_kmer = channel
        .fromQuery('select title, filename from SRA where platform is "OXFORD_NANOPORE"', db: 'inputdb')
        .map { row ->
            def (title, filename) = row
            def pacbio_file = file(filename)
            return [title, pacbio_file]
        }
        .groupTuple()
        .flatMap { title, filePairs ->
            def allFiles = filePairs.collect { it.toString().split(':') }.flatten()
            // Create a new entry for each kmerval
            kmerValues.collect { kmerval -> [title, allFiles, kmerval, "ONT", "${params.topfolder}/rawdata/kmers"] }
        }
        .view()

    
    illumina_kmer = channel
    .fromQuery('select title, filename from SRA where platform is "ILLUMINA"', db: 'inputdb')
    .map { row ->
        def (title, filename) = row
        def pacbio_file = file(filename)
        return [title, pacbio_file]
    }
    .groupTuple()
    .flatMap { title, filePairs ->
        def allFiles = filePairs.collect { it.toString().split(':') }.flatten()
        // Create a new entry for each kmerval
        kmerValues.collect { kmerval -> [title, allFiles, kmerval, "Illumina","${params.topfolder}/rawdata/kmers"] }
    }
    .view()

    //Run K_mer analysis
    kmer_channel = pacbio_kmer.mix(illumina_kmer,ont_kmer)
    kmer_histos = kmers(kmer_channel)
    
    //Plot the K-mer files
    //kmer_plotting("${params.topfolder}/rawdata/kmers","${params.topfolder}/rawdata/kmers",kmer_histos)

    //Prepare ONT files for HifiASM
    ont_asm = channel
        .fromQuery('select title, filename from SRA where platform is "OXFORD_NANOPORE"', db: 'inputdb')
        .map { row ->
            def (title, filename) = row
            def pacbio_file = file(filename)
            return [title, pacbio_file]
        }
        .groupTuple()
        .map { title, files -> 
            return [title, files.toList()]
        }
        .view()
    
    //Prepare HiFi files for HifiASM
    pacbio_asm = channel
        .fromQuery('select title, filename from SRA where platform is "PACBIO_SMRT"', db: 'inputdb')
        .map { row ->
            def (title, filename) = row
            def pacbio_file = file(filename)
            return [title, pacbio_file]
        }
        .groupTuple()
        .map { title, files -> 
            return [title, files.toList()]
        }
        .view()
    
    //Prepare HiC files for HifiASM
    hic_asm = channel
        .fromQuery('select title, filename from SRA where platform is "ILLUMINA"', db: 'inputdb')
        .map { row ->
            def (title, filename) = row
            def (file1, file2) = filename.split(':')
            return [title, file1, file2]
        }
        .groupTuple()
        .map { title, files1, files2 -> 
            return [title, files1.toList(), files2.toList()]
        }
        .view()

    //Currently ASM is deactivated since it fails on small test data (nothing to assemble)
    //and a corresponding test FASTA is provided for dev purpose
    //assembly(pacbio_asm,ont_asm,hic_asm,"${params.topfolder}/assembly")

    def fasta = "/g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta"

    //Prepare HiC files for scaffolding
    r1_scaff = hic_asm.map { title, files1, files2 ->
    return [files1]}.flatten()
    r2_scaff = hic_asm.map { title, files1, files2 ->
    return [files2]}.flatten()

    //Align HiC to assembly 
    hicalignment = arima_mapping(fasta,r1_scaff,r2_scaff,"${params.topfolder}/scaffolding/arima")

    //Scaffold with YAHS
    scaffolds = yahs(hicalignment,fasta,"${params.topfolder}/scaffolding/yahs")

    //Generate Juicer HiC map for JuiceBox visualisation
    hicmap = generate_hicmap(scaffolds,r1_scaff,r2_scaff,"${params.topfolder}/scaffolding/yahs_hicmap")


}
    







