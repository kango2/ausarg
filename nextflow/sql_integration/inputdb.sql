CREATE TABLE SRA (
    biosample_accession TEXT,
    library_ID TEXT,
    title TEXT,
    flowcell TEXT,
    library_strategy TEXT CHECK(library_strategy IN (
        'WGA', 'WGS', 'WXS', 'RNA-Seq', 'miRNA-Seq', 'WCS', 'CLONE', 'POOLCLONE', 'AMPLICON',
        'CLONEEND', 'FINISHING', 'ChIP-Seq', 'MNase-Seq', 'DNase-Hypersensitivity', 'Bisulfite-Seq',
        'Tn-Seq', 'EST', 'FL-cDNA', 'CTS', 'MRE-Seq', 'MeDIP-Seq', 'MBD-Seq', 'Synthetic-Long-Read',
        'ATAC-seq', 'ChIA-PET', 'FAIRE-seq', 'Hi-C', 'ncRNA-Seq', 'RAD-Seq', 'RIP-Seq', 'SELEX',
        'ssRNA-seq', 'Targeted-Capture', 'Tethered Chromatin Conformation Capture', 'OTHER'
    )),
    library_source TEXT CHECK(library_source IN (
        'GENOMIC', 'TRANSCRIPTOMIC', 'METAGENOMIC', 'METATRANSCRIPTOMIC', 'SYNTHETIC', 'VIRAL RNA',
        'GENOMIC SINGLE CELL', 'TRANSCRIPTOMIC SINGLE CELL', 'OTHER'
    )),
    library_selection TEXT CHECK(library_selection IN (
        'RANDOM', 'PCR', 'RT-PCR', 'HMPR', 'MF', 'CF-S', 'CF-M', 'CF-H', 'CF-T', 'MDA', 'MSLL',
        'cDNA', 'ChIP', 'MNase', 'DNAse', 'Hybrid Selection', 'Reduced Representation',
        'Restriction Digest', '5-methylcytidine antibody', 'MBD2 protein methyl-CpG binding domain',
        'CAGE', 'RACE', 'size fractionation', 'Padlock probes capture method', 'other', 'unspecified',
        'cDNA_oligo_dT', 'cDNA_randomPriming', 'Inverse rRNA', 'Oligo-dT', 'PolyA', 'repeat fractionation'
    )),
    library_layout TEXT,
    platform TEXT CHECK(platform IN (
        '_LS454', 'ABI_SOLID', 'BGISEQ', 'CAPILLARY', 'COMPLETE_GENOMICS', 'HELICOS', 'ILLUMINA',
        'ION_TORRENT', 'OXFORD_NANOPORE', 'PACBIO_SMRT'
    )),
    instrument_model TEXT CHECK(
        (platform = '_LS454' AND instrument_model IN (
            '454 GS', '454 GS 20', '454 GS FLX', '454 GS FLX+', '454 GS FLX Titanium', '454 GS Junior'
        ))
        OR (platform = 'ILLUMINA' AND instrument_model IN (
            'HiSeq X Five', 'HiSeq X Ten', 'Illumina Genome Analyzer', 'Illumina Genome Analyzer II',
            'Illumina Genome Analyzer IIx', 'Illumina HiScanSQ', 'Illumina HiSeq 1000',
            'Illumina HiSeq 1500', 'Illumina HiSeq 2000', 'Illumina HiSeq 2500', 'Illumina HiSeq 3000',
            'Illumina HiSeq 4000', 'Illumina iSeq 100', 'Illumina NovaSeq 6000', 'Illumina MiniSeq',
            'Illumina MiSeq', 'NextSeq 500', 'NextSeq 550', 'NextSeq 1000', 'NextSeq 2000',
            'Illumina HiSeq X'
        ))
        OR (platform = 'HELICOS' AND instrument_model = 'Helicos HeliScope')
        OR (platform = 'ABI_SOLID' AND instrument_model IN (
            'AB 5500 Genetic Analyzer', 'AB 5500xl Genetic Analyzer', 'AB 5500x-Wl Genetic Analyzer',
            'AB SOLiD 3 Plus System', 'AB SOLiD 4 System', 'AB SOLiD 4hq System', 'AB SOLiD PI System',
            'AB SOLiD System', 'AB SOLiD System 2.0', 'AB SOLiD System 3.0'
        ))
        OR (platform = 'COMPLETE_GENOMICS' AND instrument_model = 'Complete Genomics')
        OR (platform = 'ION_TORRENT' AND instrument_model IN (
            'Ion Torrent PGM', 'Ion Torrent Proton', 'Ion Torrent S5 XL', 'Ion Torrent S5',
            'Ion Torrent Genexus'
        ))
        OR (platform = 'PACBIO_SMRT' AND instrument_model IN (
            'PacBio RS', 'PacBio RS II', 'PacBio Sequel', 'PacBio Sequel II'
        ))
        OR (platform = 'CAPILLARY' AND instrument_model IN (
            'AB 310 Genetic Analyzer', 'AB 3130 Genetic Analyzer', 'AB 3130xL Genetic Analyzer',
            'AB 3500 Genetic Analyzer', 'AB 3500xL Genetic Analyzer', 'AB 3730 Genetic Analyzer',
            'AB 3730xL Genetic Analyzer'
        ))
        OR (platform = 'OXFORD_NANOPORE' AND instrument_model IN (
            'GridION', 'MinION', 'PromethION'
        ))
        OR (platform = 'BGISEQ' AND instrument_model IN (
            'BGISEQ-500', 'DNBSEQ-G400', 'DNBSEQ-T7', 'DNBSEQ-G50', 'MGISEQ-2000RS'
        ))
    ),
    design_description TEXT,
    filetype TEXT,
    filename TEXT,
    filename2 TEXT,
    filename3 TEXT,
    filename4 TEXT,
    assembly TEXT,
    fasta_file TEXT
);
