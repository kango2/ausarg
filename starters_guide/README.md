# New to Bioinformatics or Genome Assembly in general?
This starters guide will cover the basics and essential steps you need to know about. Look for the table of contents and jump to the topic you'd like to learn more about. 

# Table of Contents

- [Sequence Alignment](#sequence-alignment)


# Sequence Alignment 
## The Significance

By definition, Sequence Alignment is the process of comparing and detecting similarities between biological sequences [1]. In the context of our task, we are using this process between the raw sequencing data and the genome assemblies. When used as a validation technique for the assembly, we align the raw reads back to the assembly - and get alignment patterns as the result. These patterns can be further analysed and examined to obtain errors and inconsistencies in the assembly. Some of these can include :

1.  Incomplete coverage : In this case, some regions of the assembly will have low read coverage or will not be covered by reads at all. Ideally, each base or region of the genome should be covered by multiple reads to ensure that the sequence is accurately represented.
2.  Misassemblies : During the assembly process, several contigs are assembled into longer scaffolds or chromosomes. Here, misassembly errors can occur at any step and and lead to incorrect contig connections and gaps in the sequence.

## The Process

Our raw reads come from mainly three DNA sequencing technologies - PacBio, ONT and Illumina. PacBio & ONT classify as long-read sequencing technologies, whereas Illumina produces short reads. The sequence alignment for both of them is essentially the same, but with minor variations in terms of how the file are handled. The result of running a sequence alignment between the raw read data and the reference (assembly FASTA, in our case) is a .SAM file - which means a “Sequence Alignment Map” file. It stores the mapping of raw reads to the reference. Each line represents a single read that has been aligned to a reference genome. The file contains information about the read, its alignment position, and the quality of the alignment. However, a .SAM file by itself would be tricky to work with considering the large size and enourmous content associated with the file. For efficient processing and analysis of these files, we need to sort and index them.

-   Sorting : Sorting a SAM file involves reordering the reads in the file based on their alignment position. This means that reads that map to the same genomic region are grouped together. By sorting the SAM file in this way, it becomes much easier to identify reads that overlap specific regions of the genome.
-   Indexing : Indexing a SAM file involves creating an index file that maps the reference genome's position to the corresponding reads in the SAM file. Indexing allows faster access to specific regions of the reference genome.