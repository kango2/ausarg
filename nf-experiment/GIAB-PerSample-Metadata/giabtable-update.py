import os
import sqlite3

# Define the varying data for each sample in a dictionary
sample_data_dict = {
    "HG001": ("Homo sapiens (human)", "Blood", "B-Lymphocyte", "Female", "NA", "White", "Uta/Mormon", "Mother", "Healthy", "NA12878", "HG001", "PAIRED", "WGS PCR-free", "Genomics", "ILLUMINA", "ILLUMINA HiSeq 2500 Rapid SBS", "300x", "564.56 bp", "2 x 148bp"),
    "HG002": ("Homo sapiens (human)", "Blood", "B-Lymphocyte", "Male", "46YR", "White", "Ashkenazim Jewish", "Son", "Healthy", "NA24385", "HG002", "PAIRED", "WGS PCR-free", "Genomics", "ILLUMINA", "ILLUMINA HiSeq 2500 Rapid SBS", "300x", "564.56 bp", "2 x 148bp"),
    "HG003": ("Homo sapiens (human)", "Blood", "B-Lymphocyte", "Male", "90YR", "White", "Ashkenazim Jewish", "Father", "Healthy", "NA24149", "HG003", "PAIRED", "WGS PCR-free", "Genomics", "ILLUMINA", "ILLUMINA HiSeq 2500 Rapid SBS", "300x", "564.56 bp", "2 x 148bp"),
    "HG004": ("Homo sapiens (human)", "Blood", "B-Lymphocyte", "Male", "74YR", "White", "Ashkenazim Jewish", "Mother", "Healthy", "NA24143", "HG004", "PAIRED", "WGS PCR-free", "Genomics", "ILLUMINA", "ILLUMINA HiSeq 2500 Rapid SBS", "300x", "564.56 bp", "2 x 148bp"),
    "HG005": ("Homo sapiens (human)", "Blood", "B-Lymphocyte", "Male", "33YR", "Asian", "Chinese", "Son", "Healthy", "NA24631", "HG005", "PAIRED", "WGS PCR-free", "Genomics", "ILLUMINA", "ILLUMINA HiSeq 2500 Rapid SBS", "300x", "564.56 bp", "2 x 250 bp"),
    "HG006": ("Homo sapiens (human)", "Blood", "B-Lymphocyte", "Male", "64YR", "Asian", "Chinese", "Father", "Healthy", "NA24694", "HG006", "PAIRED", "WGS PCR-free", "Genomics", "ILLUMINA", "ILLUMINA HiSeq 2500 Rapid SBS", "100x", "564.56 bp", "2 x 148 bp"),
    "HG007": ("Homo sapiens (human)", "Blood", "B-Lymphocyte", "Male", "63YR", "Asian", "Chinese", "Mother", "Healthy", "NA24695", "HG007", "PAIRED", "WGS PCR-free", "Genomics", "ILLUMINA", "ILLUMINA HiSeq 2500 Rapid SBS", "100x", "564.56 bp", "2 x 148 bp")
}

def extract_info_from_fastq_name(fastq_path):
    fastq_name = os.path.basename(fastq_path)
    parts = fastq_name.split('_')
    date_of_sequencing = parts[1]
    flow_cell_number = parts[3]
    flow_cell_id = parts[4]
    return date_of_sequencing, flow_cell_number, flow_cell_id

def write_results_to_database(fastq_directory, db_path, sample_id):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    for file_name in os.listdir(fastq_directory):
        if file_name.endswith("_R1.fastq.gz"):
            fastq_file1 = os.path.join(fastq_directory, file_name)
            fastq_file2 = os.path.join(fastq_directory, file_name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            date_of_sequencing, flow_cell_number, flow_cell_id = extract_info_from_fastq_name(fastq_file1)

            # Use the sample_id to get the appropriate data
            sample_data = sample_data_dict[sample_id]

            cursor.execute(
            """
                INSERT OR REPLACE INTO giab_study (
                    fastq_file1, fastq_file2, organism, tissue_type, cell_type, sex, age, race, ethnicity, 
                    relationship, condition, coriell_id, nist_id, library_layout, 
                    library_strategy, library_source, instrument_platform, instrument_model, 
                    sequence_depth, mean_paired_end_distance, mean_read_length,
                    date_of_sequencing, flow_cell_number, flow_cell_id
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, 
            (
                fastq_file1, fastq_file2, *sample_data, date_of_sequencing, flow_cell_number, flow_cell_id
            )
            )

    conn.commit()
    conn.close()

def create_database(db_path):
    # Connect to the database (this will create it if it doesn't exist)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create the table structure
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS giab_study (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        fastq_file1 TEXT,
        fastq_file2 TEXT,
        organism TEXT,
        tissue_type TEXT,
        cell_type TEXT,
        sex TEXT,
        age TEXT,
        race TEXT,
        ethnicity TEXT,
        relationship TEXT,
        condition TEXT,
        coriell_id TEXT,
        nist_id TEXT,
        library_layout TEXT,
        library_strategy TEXT,
        library_source TEXT,
        instrument_platform TEXT,
        instrument_model TEXT,
        sequence_depth TEXT,
        mean_paired_end_distance TEXT,
        mean_read_length TEXT,
        date_of_sequencing TEXT,
        flow_cell_number TEXT,
        flow_cell_id TEXT
    )
    ''')

    # Commit the changes and close the connection
    conn.commit()
    conn.close()

# Specify the SQLite database path
db_path = "/g/data/xy86/giab/giab.db"

create_database(db_path)

# Loop over all sample IDs and call write_results_to_database for each
for sample_id in sample_data_dict.keys():
    fastq_directory = f"/g/data/xy86/giab/{sample_id}"
    write_results_to_database(fastq_directory, db_path, sample_id)
