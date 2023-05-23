# Genome Assembly Pipeline 

A comprehensive overview of the pipeline, including the explanation of the steps involved - along with softwares used and their versions. 

Table of Contents
=================

- [Data Management](#data-management)
- [Data Prepatation](#data-preparation)
- [Assembly Generation](#assembly-generation)
- [Assembly Evaluation](#assembly-evaluation)

---
## Data Management 
### File Directories 
To arrange and organise our files on NCI Gadi, we derived inspiration from the structure followed by VGP (https://vertebrategenomesproject.org). Here is a rough sturture followed by VGP -

---

**VGP File Structure**

- species
  - specie_name
    - version
      - Automated Assembly
        - ...
        - Intemediates
          - hifiasm
          - ...
      - Curated Assembly
         - ...
        - Evaluation
          - busco
          - genomescope
          - merqury
          - pretest
          - quast
      - Genomic Data
        - raw data

---

**VGP Nomenclature**

- Specie Names : [amr..]SpeNme[1,2,3]
where 
  - [amr..] represents the specie type 
  - SpeNme are the initials of the specie 
  - [1,2,3] refers to the version number 

---

Following the above convention, an example name generated would be `bAmmCau1` 

- **Assembly Names** : [amr..]SpeNme[1,2,3] . org . asmType . date. fasta . gz 
where 
  - [amr..]SpeNme[1,2,3] is the specie name, as explained above 
  - org is the abbreviation of the organisation which worked on the assembly 
  - asmType is the type of assembly, such as "cur" for curated
  - date is the date on which the assembly was generated

Following the above convention, an example assembly name generated would be `bAmmCau1.pri.cur.20230231.fasta.gz`

**We can observe the following key points :**

 -  The automated assembly results and the curated assembly are stored in separate folders, resulting in an extensive representation of the process followed to reach the final results 
 - The automated assembly folder has a repository of intermediate files from all the tools and softwares used, along with a metadata file indicating the assembly process, software versions and some other details. 
 - The curated assembly folder has an evaluation folder with results from a variety of evaluation tools.
 - Lastly,  all the raw genomic data is stored along with the assembly folders. 

Taking inspiration from the above structure, we aspire to replicate it for our data. Since we are only getting started with the process, we have temporarily organised our data as follows:

---

**Processed Files**

- species
    - specie_name
        - version
            - assembly
            - metadata
            - evaluation

---

**Raw Files**

- bpadata
    - specie_name
        - raw
            - Pacbio
            - ONT
            - Illumina
        - clean
            - pacbio
            - ONT
            - Illumina


---
 - Our raw data is currently residing at `/g/data/xl04/bpadata` on NCI Gadi, and the processing is currently being conducted by Kirat Alreja (Bioinformatics Support Officer, ANU)  in his folder, `ka6418`. We will eventually move this to a shared folder similar to `bpadata`
 - We are yet to download and organise all the raw data as well, so the whole process is going hand-in-hand with manual curation & evaluation.

 ### MD5 Metadata Tables 

 We are working towards creating two comprehensive metadata table for our assembly file, which serve as a combined integrity and evaluation check. There are two tables, Sequence Table and Assembly Table. 

#### Sequence Table

| Assembly_ID (FK) | Sequence_ID (PK) | Length | Order In File | MD5 Sumcheck |
|----------|----------|----------|----------|----------|
|          |          |          |          |          |
|          |          |          |          |          |
|          |          |          |          |          |
|          |          |          |          |          |
|          |          |          |          |          |


 - **Assembly ID** : Name of the Assembly file, it acts as there foreign key and will eventually help us set up a connection between the Sequence Table and the Assembly Table 
 - **Sequence ID** : Name of every sequence of the Assembly file 
 - **Length** : Length of every individual sequence 
 - **Order In File** : The current order of the sequence inside the assembly file - will help us maintain a visualisation of how the sequences are rearranged
 - **MD5 Sumcheck** : MD5 for every sequence 

#### Assembly Table 

| Assembly_ID (PK) | Length | Uncompressed MD5 | Compressed MD5 | Metrics |
|----------|----------|----------|----------|----------|
|          |          |          |          |          |
|          |          |          |          |          |
|          |          |          |          |          |
|          |          |          |          |          |
|          |          |          |          |          |


 - **Assembly ID** : Name of the Assembly File 
 - **Length** : Length of the entire assembly 
 - **Uncompressed MD5** : MD5 sumcheck of the .fasta / uncompressed file 
 - **Compressed MD5** : MD5 sumcheck of the .gz / compressed file, -   The reasoning behind having both Uncompressed MD5 and Compressed MD5 is that the sumcheck changes when you uncompress a .GZ file, and also when you recompress the same. Having both the metrics ensures that we can verify integrity in the case if anyone recompresses an existing .fasta file.
 - **Metrics** : Several columns with assembly evaluation metrics, discussed below 

#### Metrics 

##### General Metrics
-   N50: The length of the shortest sequence that covers 50% of the total assembly, indicating the  **continuity**  of the assembly.
-   N90: The length of the shortest sequence that covers 90% of the total assembly, providing  **a more stringent measure of assembly continuity.**
-   L50: The number of sequences required to reach 50% of the assembly length, reflecting the  **fragmentation of the assembly**.
-   L90: The number of sequences required to reach 90% of the assembly length, providing  **a more stringent measure of assembly fragmentation.**
-   #N per 100kbp: The number of N bases (representing gaps) per 100 kilobase pairs, indicating the  **completeness of the assembly**  and the presence of gaps.
-   Largest scaffold: The length of the longest sequence in the assembly, reflecting the  **ability of the assembly to capture large genomic regions.**
-   Total number of scaffolds: The number of contigs/scaffolds in the assembly,  **reflecting the level of fragmentation of the assembly.**
-   GC%: The percentage of the genome that consists of the nucleotides G and C, which can affect gene expression and function, as well as provide information on the  **quality of the sequencing and assembly processes**.

##### Specialised Metrics 
-   **BUSCO score**  - The BUSCO score reflects the completeness and correctness of the genome assembly by comparing it to a set of conserved genes. A higher BUSCO score is preferred, indicating a more complete and accurate assembly.
-   **Error rate estimations**  - Error rate estimations can provide an estimate of the quality of the assembly in terms of base pair accuracy and gene structure.
-   **Heterozygosity rate**  - The heterozygosity rate can impact the quality of the assembly, especially for species with high levels of genetic diversity. High heterozygosity rates can lead to higher fragmentation and lower contiguity of the assembly.

## NCBI Metadata

>To Do, consult with Hardip 

---
# Data Prepatation
We are using DNA Sequencing data from five primary technologies, to benefit from the complementary strengths of each approach.
- **PacBio Hifi** - Long reads that can span repetitive and complex genomic regions. It is useful for resolving structural variations and producing high-quality consensus sequences. Used for accurate and contiguous genome assemblies. 
- **ONT Ultralong** - Long reads similar to PacBio, but offers additional sequence information to account for genomic gaps and resolving repetitive regions. 
- **Illumina** - Used for generating short reads at high throughput. We utilise it for generating large amounts of data to aid in error correction, genome polishing and validation of the final assembly.
- **HiC** - Used for scaffolding contigs and improving the overall accuracy and contiguity of the assembly. Provides informaiton about >to do
- **RNASeq** - Used for gene prediction, annotation and validation of the assembly and particularly useful for providing insights into the functional elements of the genome. 

Here is a list of tools we used to pre-process our sequencing data (adapter removal & quality filtering), along with explanations and their versions. 

| Tool            | Version     | Reads Type                              |
|-----------------|-------------|------------------------------------------|
| DeepConsensus     | TODO       | PacBioHifi |
| Guppy/Dorado   | TODO      | ONT       |
| Trimmomatic | TODO | Illumina |

