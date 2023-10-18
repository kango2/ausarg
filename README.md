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

input: Path to the input fasta file.
output: Path to the output directory where the CSV file will be saved.
permatch: Percentage match for the telomeric repeats. The recommended setting is 90.
copies: Minimum number of copies of the telomeric repeat. The recommended setting is 100.

### Dependencies

The script requires the following modules:

- kentutils
- TRF (Tandem Repeats Finder)
- biopython
- parallel
- Additionally, it uses a Python script trf2gff.py to convert TRF output to GFF3 format and another Python script clean_telomere_csv.py to process the results.

### Output

The script generates a CSV file with the following columns:

| Column Name     |
|-----------------|
| Sequence_ID     |
| Start           |
| End             |
| ID              |
| period          |
| copies          |
| consensus_size  |
| perc_match      |
| perc_indels     |
| align_score     |
| entropy         |
| cons_seq        |
| repeat_seq      |
