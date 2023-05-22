#Genome Assembly Pipeline 

A comprehensive overview of the pipeline, including the explanation of the steps involved - along with softwares used and their versions. 

---
##Data Management 
###File Directories 
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

- Assembly Names : [amr..]SpeNme[1,2,3] . org . asmType . date. fasta . gz 
where 
  - [amr..]SpeNme[1,2,3] is the specie name, as explained above 
  - org is the abbreviation of the organisation which worked on the assembly 
  - asmType is the type of assembly, such as "cur" for curated
  - date is the date on which the assembly was generated

Following the above convention, an example assembly name generated would be `bAmmCau1.pri.cur.20230231.fasta.gz`

We can observe the following key points :

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