# New to Bioinformatics or Genome Assembly in general?

This starters guide will cover the basics and essential steps you need to know about. Look for the table of contents and jump to the topic you'd like to learn more about.

# Table of Contents

  

- [Sequence Alignment](#sequence-alignment)

  

# Sequence Alignment

  

## The Significance

  

By definition, Sequence Alignment is the process of comparing and detecting similarities between biological sequences. In the context of our task, we are using this process between the raw sequencing data and the genome assemblies. When used as a validation technique for the assembly, we align the raw reads back to the assembly - and get alignment patterns as the result. These patterns can be further analysed and examined to obtain errors and inconsistencies in the assembly. Some of these can include :

  

1.  **Incomplete coverage** : In this case, some regions of the assembly will have low read coverage or will not be covered by reads at all. Ideally, each base or region of the genome should be covered by multiple reads to ensure that the sequence is accurately represented.

2.  **Misassemblies** : During the assembly process, several contigs are assembled into longer scaffolds or chromosomes. Here, misassembly errors can occur at any step and and lead to incorrect contig connections and gaps in the sequence.

  

## The Process

  

Our raw reads come from mainly three DNA sequencing technologies - PacBio, ONT and Illumina. PacBio & ONT classify as long-read sequencing technologies, whereas Illumina produces short reads. The sequence alignment for both of them is essentially the same, but with minor variations in terms of how the file are handled. The result of running a sequence alignment between the raw read data and the reference (assembly FASTA, in our case) is a .SAM file - which means a **“Sequence Alignment Map”** file.

  

It stores the mapping of raw reads to the reference. Each line represents a single read that has been aligned to a reference genome. The file contains information about the read, its alignment position, and the quality of the alignment. However, a .SAM file by itself would be tricky to work with considering the large size and enourmous content associated with the file. For efficient processing and analysis of these files, we need to sort and index them.

  

-  **Sorting** : Sorting a SAM file involves reordering the reads in the file based on their alignment position. This means that reads that map to the same genomic region are grouped together. By sorting the SAM file in this way, it becomes much easier to identify reads that overlap specific regions of the genome.

-  **Indexing** : Indexing a SAM file involves creating an index file that maps the reference genome's position to the corresponding reads in the SAM file. Indexing allows faster access to specific regions of the reference genome.

  

Another caveat that comes into play is that we may have multiple .SAM file for the same specie. This might be a result of “Sample Multiplexing” or “Repetitive Sequencing”. In many sequencing experiments, multiple samples are sequenced together in a single run, and the resulting reads are split into different files based on their sample of origin. This process is called **“Sample Multiplexing”**. On the other hand, **“Repetitive Sequencing”** is done when the same sample is sequenced multiple times using different sequencing technologies or protocols, resulting in multiple SAM files. In both the cases, we need to merge these individual .SAM file to proceed with our analysis & evaluation

  

Lastly, when dealing with a large amount of data, we convert the resulting .SAM file into a .BAM file. The reason is to reduce the file size and improve processing speed. SAM files are text files that contain a lot of information about each read, including the read sequence, quality scores, alignment position, and other optional fields. While this level of detail is useful, it can also make the files quite large, especially for datasets with millions or billions of reads.

  

BAM files, on the other hand, are binary versions of SAM files. They contain the same information but are stored in a more compact format that makes them faster to read and write. BAM files use a combination of binary encoding and compression to reduce the file size without losing any information.

  

## Evaluation

  

One of the most significant evaluation measures using sequence alignment is **“read depth”**. Read depth is a measure of the number of times a given nucleotide in a DNA sequence has been read by a sequencing technology. In other words, it is the number of times a particular position in a genome has been sequenced or covered by reads. Read depth is important because it provides information about the accuracy and completeness of a sequence assembly and the reliability of variant detection.

  
Once we get the final .BAM file, we can proceed to dive into evaluation. The first step would be analyzing the “read depth”, which can be done using a .txt coverage file generated by `samtools`.

    samtools depth sample.bam > coverage.txt
The resulting .txt file outputs the read depth at each position in the genome, sorted by chromosome and position. The output file contains three columns: the chromosome name, the position in the chromosome, and the read depth at that position.

Here is how this file looks like: 

    haplotype1-0000001      1       36
    haplotype1-0000001      2       36
    haplotype1-0000001      3       36
    haplotype1-0000001      4       36
    haplotype1-0000001      5       36

