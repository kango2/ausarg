



for i in `sqlite3 /g/data/xl04/ka6418/ausarg/database/ausarg.db 'select file_path from Illumina_Metrics where status == "Pending"'`; do j=$(echo $i | sed 's/;/:/'); qsub -v filepair=$j,output=/g/data/xl04/bpadata/Pogona_vitticeps/raw/evaluation/illumina_qc /g/data/xl04/ka6418/ausarg/scripts/run_illumina_qc.sh; done