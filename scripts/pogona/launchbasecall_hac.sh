module load utils 

ontfiles=$(fetchPogona.sh reads POGVIT BLOW5)
logs=/g/data/xl04/genomeprojects/Pogona_vitticeps/logs
output=/g/data/xl04/bpadata/Pogona_vitticeps/raw/ONT/fastx_hac

IFS=':' 
for file in $ontfiles; do 
   qsub -o $logs -v MERGED_SLOW5=$file,BASECALL_OUT=$output /g/data/te53/ka6418/refgen/workflows/basecalling-bpa/bin/basecall.sh
done