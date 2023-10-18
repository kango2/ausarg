## find_telomeres.sh

![Description of the Image](images/find_telomeres.png)

[Flowchart link](https://whimsical.com/detailed-parallelized-flowchart-for-find-telomeres-sh-2G623E2p7ubVTCZ8DGfJ4S@2Ux7TurymNGkJUXCfVvk)

This script is designed to identify telomeres in a given fasta sequence. It uses the Tandem Repeats Finder (TRF) to detect telomeric repeats and then processes the results to generate a CSV file with the identified telomeres. The gist is that the script will find all 6BP repeats in the FASTA and then finally filter out the possible repeats based on the possible variations of TTAGGG in the forward and reverse strand. For example, the combinations are generated like GTTAGG, GGTTAG and so on. More information about how the telomeres are filtered out is in ```clean_telomeres_csv.py ```

### Usage

To use this script, submit it with the required parameters:

``` 
qsub -l storage=gdata/if89+gdata/projectcode -o /path/to/stdouterr -P projectcode -v input=/path/to/fasta,output=/path/to/output/csv,permatch=90,copies=100 ./find_telomeres.sh
```

### Parameters 

- input: Path to the input fasta file.
- output: Path to the output directory where the CSV file will be saved.
- permatch: Percentage match for the telomeric repeats. The recommended setting is 90.
- copies: Minimum number of copies of the telomeric repeat. The recommended setting is 100.

### Dependencies

The script requires the following modules:

- **kentutils 0.0**
- **TRF (Tandem Repeats Finder) 4.09.1**
- **biopython 1.79**
- **parallel 20191022**
- Additionally, it uses a Python script ```trf2gff.py``` to convert TRF output to GFF3 format and another Python script ```clean_telomere_csv.py``` to process the results.

### Expected runtime 

When using 48 CPUS and 64 GB RAM - I observed runtime of anywhere from 2 hours to 5 hours for chromosome level skink assemblies. If you assembly is fragmented, the runtime decreases exponentially, since the script makes good use of GNU parallel. 

### Output

The script generates a CSV file with the following columns:

| Sequence_ID | Start | End | ID | period | copies | consensus_size | perc_match | perc_indels | align_score | entropy | cons_seq | repeat_seq | relative start | relative end |

Most of the output columns are self-explanatory, however I added the ```relative start``` and ```relative end``` columns for understanding where the telomeric sequences are located in each contig. They are just ``` telomeric start or end coordinate \ total length of contig ```. So, if the value is 0 - that means the telomere occurs at the start and similarly 1 indicates the end. In draft assemblies you might observes values like 0.3, 0.5 which indicates the telomeres lie in the 30%/50% positions. 



