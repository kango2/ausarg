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

```
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
```
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
```
- species
    - specie_name
        - version
            - assembly
            - metadata
            - evaluation
```
---

**Raw Files**
```
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

```
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
 - **Metrics** : Several columns with assembly evaluation metrics, discussed in [Assembly Evaluation](#assembly-evaluation)


#### NCBI Metadata

>To Do, consult with Hardip 

---
## Data Preparation
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

## Assembly Generation 

Similar to the data preparation step, we use a number of tools & techniques to benefit from the complimentary strengths of each approach. 

1. **Verkko**: This is an iterative, graph-based pipeline for assembling complete, diploid genomes. It begins with a multiplex de Bruijn graph built from long, accurate reads, progressively simplifies this graph by integrating ultra-long reads and haplotype-specific markers. The result is a phased, diploid assembly of both haplotypes, with many chromosomes automatically assembled from telomere to telomere. In a test on the HG002 human genome, 20 of 46 diploid chromosomes were assembled without gaps at 99.9997% accuracy​. [[1]](https://www.nature.com/articles/s41587-023-01662-6)
2. **Hifiasm**: A de novo assembler that takes advantage of long high-fidelity sequence reads to represent the haplotype information in a phased assembly graph. It strives to preserve the contiguity of all haplotypes, which makes it particularly effective for studying sequence variations in a genome. [[2]](https://www.nature.com/articles/s41592-020-01056-5)
3. **Canu**: This tool specializes in assembling PacBio or Oxford Nanopore sequences and operates in three phases: correction, trimming, and assembly. It improves the accuracy of bases in reads, trims reads to high-quality sequences, and orders the reads into contigs. Canu can resume incomplete assemblies and will auto-detect computational resources and scale itself to fit, using all available resources. [[3]](https://canu.readthedocs.io/en/latest/quick-start.html)
4. **Flye**: This is a de novo assembler for single-molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies. It's designed for a wide range of datasets, from small bacterial projects to large mammalian-scale assemblies. It uses a repeat graph as the core data structure, which is built using approximate sequence matches and can tolerate the higher noise of SMS reads. However, it typically produces collapsed assemblies of diploid genomes, represented by a single mosaic haplotype. [[4]](https://github.com/fenderglass/Flye)


[1] : Telomere-to-Telomere consortium. (2023). "The Telomere-to-Telomere consortium recently assembled the first truly complete sequence of a human genome." Available at: https://www.nature.com/articles/s41587-023-01662-6 [Accessed 23 May 2023]​

[2] : Hifiasm Team. (2023). "Haplotype-resolved de novo assembly is the ultimate solution to the study of sequence variations in a genome." Available at: https://www.nature.com/articles/s41592-020-01056-5 [Accessed 23 May 2023]

[3] : Canu Documentation. (2023). "Canu specializes in assembling PacBio or Oxford Nanopore sequences." Available at: https://canu.readthedocs.io/en/latest/quick-start.html [Accessed 23 May 2023]

[4] : Flye Development Team. (2023). "Flye is a de novo assembler for single-molecule sequencing reads." Available at: https://github.com/fenderglass/Flye [Accessed 23 May 2023]​

Tool versions we used for development :

| Tool | Version  | 
|----------|----------|
|   Verkko       |          |
|   HifiASM       |          |
|   Canu       |          |
|   Flye       |          |



## Assembly Evaluation 

Our assembly evaluation workflow revolves around the 3C criterion of genome assembly evaluation : 

1. **Contiguity**: It evaluates the assembly in terms of number and size of contigs and scaffolds, the pieces found in an assembly. Metrics includes statistics related to maximum length, average length, combined total length, and contig N50 (length-weighted median of ordered contigs or scaffolds). However, contiguity metrics thereof need to be interpreted with caution due they do not contain information on assembly accuracy and completeness.

![Contiguity](https://i.imgur.com/rL5kL9n.png)

2. **Correctness**: it refers to how well those pieces accurately represent the genome sequenced and, in general is acceptable that it is essential to prioritize correctness rather than contiguity. However, correctness is difficult to evaluate if a preliminary reference genome is not available, which is a particular problem for de novo assembly. Mapping and comparison to reference or draft genome (or a consensus sequence) can be used to detect misassemblies, including mismatches, indels, and misjoins.

![Correctness](https://i.imgur.com/FIKfxxh.png)

3. **Completeness**: it assesses how much of the genome is represented by the pieces of the assembly. This implies the evaluation of ability to assembly not only all the genes, but also to solve all complicated regions, including repetitive sequences and, if it is expected, circularization of genome. The most important metric for this case is the “completeness score”, calculated by the examination of single-copy orthologs conserved genes. In addition, information of known sequences, unexpected variations in coverage, and remapping of reads allows to analyze the consistency of the genome and identification of potentially poorly assembled regions.

![Completeness](https://i.imgur.com/ksn7LS3.png)

TODO : reference - https://www.nature.com/articles/s41598-020-58319-6

We assess all these three criteries with a vareity of metrics and tools : 

| Tool Name                      | Evaluation Type   |
|-------------------------------|-------------------|
| N50                           | Contiguity        |
| N90                           | Contiguity        |
| L50                           | Contiguity        |
| L90                           | Contiguity        |
| #N per 100kbp                 | Completeness      |
| Largest scaffold              | Contiguity        |
| Total number of scaffolds     | Contiguity        |
| GC%                           | Completeness      |
| BUSCO score                   | Completeness      |
| Error rate estimations        | Correctness       |
| Heterozygosity rate           | Correctness       |

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








