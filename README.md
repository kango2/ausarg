# File organisation

To arrange and organise our files on NCI Gadi, we derived inspiration from the structure followed by VGP (https://vertebrategenomesproject.org). Here is a rough sturture followed by VGP -

  
![VGP structrure](https://i.ibb.co/cQpfG6n/Screenshot-2023-05-08-at-1-18-39-pm.png)

![VGP Nomenclature
](https://i.ibb.co/gjLVpbK/Screenshot-2023-05-08-at-1-27-11-pm.png)

We can observe the following key points :

 -  The automated assembly results and the curated assembly are stored in separate folders, resulting in an extensive representation of the process followed to reach the final results 
 - The automated assembly folder has a repository of intermediate files from all the tools and softwares used, along with a metadata file indicating the assembly process, software versions and some other details. 
 - The curated assembly folder has an evaluation folder with results from a variety of evaluation tools.
 - Lastly,  all the raw genomic data is stored along with the assembly folders. 

Taking inspiration from the above structure, we aspire to replicate it for our data. Since we are only getting started with the process, we have temporarily organised our data as follows:

![AusARG files](https://i.ibb.co/tZwhGJk/Screenshot-2023-05-08-at-1-40-12-pm.png%5B)

 - Our raw data is currently residing at "/g/data/xl04/bpadata" on NCI Gadi, and the processing is currently being conducted by Kirat Alreja (Bioinformatics Support Officer, ANU)  in his folder, "ka6418". We will eventually move this to a shared folder similar to "bpadata"
 - We are yet to download and organise all the raw data as well, so the whole process is going hand-in-hand with manual curation & evaluation. 

# Metadata Tables
We are working towards creating two comprehensive metadata table for our assembly file, which serve as a combined integrity and evaluation check. There are two tables, Sequence Table and Assembly Table. 

## Sequence Table
![Sequence Table](https://i.ibb.co/Lhj97V8/Screenshot-2023-05-08-at-2-12-57-pm.png)

 - **Assembly ID** : Name of the Assembly file, it acts as there foreign key and will eventually help us set up a connection between the Sequence Table and the Assembly Table 
 - **Sequence ID** : Name of every sequence of the Assembly file 
 - **Length** : Length of every individual sequence 
 - **Order In File** : The current order of the sequence inside the assembly file - will help us maintain a visualisation of how the sequences are rearranged
 - **MD5 Sumcheck** : MD5 for every sequence 

## Assembly Table 
![assembly table](https://i.ibb.co/ZJzLn91/Screenshot-2023-05-08-at-3-39-49-pm.png)
 
 

 - **Assembly ID** : Name of the Assembly File 
 - **Length** : Length of the entire assembly 
 - **Uncompressed MD5** : MD5 sumcheck of the .fasta / uncompressed file 
 - **Compressed MD5** : MD5 sumcheck of the .gz / compressed file, -   The reasoning behind having both Uncompressed MD5 and Compressed MD5 is that the sumcheck changes when you uncompress a .GZ file, and also when you recompress the same. Having both the metrics ensures that we can verify integrity in the case if anyone recompresses an existing .fasta file.
 - **Metrics** : Several columns with assembly evaluation metrics, discussed below 

### Metrics 

#### General Metrics
-   N50: The length of the shortest sequence that covers 50% of the total assembly, indicating the  **continuity**  of the assembly.
-   N90: The length of the shortest sequence that covers 90% of the total assembly, providing  **a more stringent measure of assembly continuity.**
-   L50: The number of sequences required to reach 50% of the assembly length, reflecting the  **fragmentation of the assembly**.
-   L90: The number of sequences required to reach 90% of the assembly length, providing  **a more stringent measure of assembly fragmentation.**
-   #N per 100kbp: The number of N bases (representing gaps) per 100 kilobase pairs, indicating the  **completeness of the assembly**  and the presence of gaps.
-   Largest scaffold: The length of the longest sequence in the assembly, reflecting the  **ability of the assembly to capture large genomic regions.**
-   Total number of scaffolds: The number of contigs/scaffolds in the assembly,  **reflecting the level of fragmentation of the assembly.**
-   GC%: The percentage of the genome that consists of the nucleotides G and C, which can affect gene expression and function, as well as provide information on the  **quality of the sequencing and assembly processes**.

#### Specialised Metrics 
- BUSCO Evaluation 
- Heterozygosity Score
- Error rate estimations 


