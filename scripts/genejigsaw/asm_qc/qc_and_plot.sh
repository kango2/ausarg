

fasta=${ref}
illumina=${illumina}
pacbio=${pacbio}
ont=${ont}
outputdir=${outputdir}
gc=/g/data/xl04/ka6418/github/ausarg/scripts/gc_content.sh 
telomere=/g/data/xl04/ka6418/github/ausarg/scripts/find_telomeres.sh
centromere=/g/data/xl04/ka6418/github/ausarg/scripts/centromeres.sh


#Launch sequence table job 

#Launch GC content job 

qsub ${gc} #fasta #output

#Launch Telomere job

qsub ${telomere}

#Launch Centromere job 

#Launch Read Depth job  





