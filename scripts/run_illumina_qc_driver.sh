$ for i in `sqlite3 /g/data/xl04/ka6418/sql/ausarg.db 'select filename from SRA where library_strategy=="WGS" and platform=="ILLUMINA"'`; do j=$(echo $i | sed 's/;/:/'); qsub -v filepair=$j,output=/g/data/xl04/bpadata/Tiliqua_rugosa/raw/evaluation/illumina_qc /g/data/xl04/ka6418/ausarg/scripts/run_illumina_qc.sh; done