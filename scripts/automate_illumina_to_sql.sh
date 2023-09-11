



for i in `sqlite3 /g/data/xl04/ka6418/ausarg/database/ausarg.db 'select file_path from Illumina_Metrics where status == "Pending"'`; do j=$(echo $i | sed 's/;/:/'); echo qsub -v filepair=$j,output=/g/data/xl04/bpadata/Bassiana_duperreyi/raw/evaluation/illumina_qc run_illumina_qc.sh; done