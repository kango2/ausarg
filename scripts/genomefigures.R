library(tidyverse)
library(scales)
library(ggrepel)
library(karyoploteR)
library(tidygraph)
library(Biostrings)
library(readxl)

projectdir <- "/g/data/xl04/genomeprojects/Pogona_vitticeps"

setwd(projectdir)
reffasta <- "./fasta/POGVIT.v2.1.fasta"
refseqtable <- "./analysis/seqtables/POGVIT.v2.1_seqtable.csv"
telomeres <- "./analysis/Telomeres/POGVIT.v2.1_Telomeres.csv"
ilmnrd <- "./analysis/depth/rearranged/POGVIT.v2.1.merged.illum.10000.depth.bed"
ontrd <- "./analysis/depth/rearranged/POGVIT.v2.1.merged.ont.10000.depth.bed"
hifird <- "./analysis/depth/rearranged/POGVIT.v2.1.merged.pb.10000.depth.bed"
gc <- "./analysis/gc/POGVIT.v2.1.fasta_GC.csv"
gaps <- "./analysis/gaps/POGVIT.v2.1_Nregions.csv"
hiczscores <- "././analysis/plots/hicinteractionszscores.tsv"
rdna <- "./analysis/rDNA/POGVIT.v2.1.rdnaregions.bed"
h1fasta <- "./fasta/yahs/POGVITdef.h1.corrected_YAHS.fasta"
h2fasta <- "./fasta/yahs/POGVITdef.h2.corrected_YAHS.fasta"
plotdir <- "./analysis/plots/"

###definitions
###tech = sequencing technology platforms, values = Illumina, PB, ONT, HiC
###

###################################################################################################
################Kmer analyses
kmerdir <- "./analysis/kmers"
setwd(kmerdir)
files <- list.files(".", pattern = "\\.histo$")
khisto <- list()
for (i in 1:length(files)) {
  f <- files[i]
  x <- str_split(str_split(f, "\\.", simplify = T)[1], "_", simplify = T)
  a <- read_delim(files[i], delim = " ", col_names = F)
  colnames(a) <- c("kcounts", "kfrequency")
  a$tech <- x[2]
  a$klen <- x[3]
  khisto[[i]] <- a
}
khisto <- bind_rows(khisto)

### get the value of the homozygous peak
peaks <- khisto %>% 
  filter(kcounts < 100) %>% 
  group_by(tech, klen) %>% 
  arrange(desc(kcounts)) %>% 
  mutate(d=kfrequency-lead(kfrequency)) %>% 
  group_by(tech, klen) %>% 
  arrange(desc(kcounts)) %>%
  filter(d>0) %>%
  filter(kcounts == max(kcounts))

### get the k-mer plot
kmerplot <- khisto %>% 
  mutate(kcounts = case_when(kcounts>200 ~ 200, TRUE ~ kcounts)) %>% 
  group_by(tech, klen, kcounts) %>% 
  summarise(kfrequency = sum(kfrequency)) %>%
  filter((kcounts>1 & kcounts <=100 & kfrequency < 5e7 & tech == "Illumina") | (kcounts>2 & kcounts <=100 & tech == "PB") | (kcounts>3 & kcounts <=100 & kfrequency < 5e7 & tech == "ONT")) %>%
  ggplot(aes(x=kcounts, y = kfrequency)) + 
  geom_area(stat="identity", fill = "lightblue") + 
  geom_vline(data = peaks, aes(xintercept = kcounts)) +
  geom_label(data = peaks, aes(x=kcounts, y = kfrequency * 1.5,  label = kcounts)) + 
  theme_bw() +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6, sep = "")) +
  facet_grid(tech~klen, scales = "free_y") +
  xlab("Kmer counts") + ylab("Frequency") +
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = 1))
kmerplot
### save the k-mer plots
setwd(projectdir)
ggsave(paste(plotdir, "kmerhisto.pdf", sep = ""), kmerplot, width = 8, height = 8)

### need to save peaks in excel or somewhere for genome size estimations
write_delim(select(peaks, Technology = tech, `K-mer Length` = klen, `K-mer Peak` = kcounts, Frequency = kfrequency) %>% arrange(`K-mer Length`, Technology), 
            paste(plotdir, "kmerpeaks.tsv", sep = ""), delim = "\t")
########################################################################################################
################Raw data QC metrics
setwd(projectdir)
rawqcdir <- "./analysis/raweval/"
setwd(rawqcdir)

pbfiles <- list.files("./PB", full.names = T, recursive = T, pattern = "*_quality_freq.csv")
ontfiles <- list.files("./ONT/", full.names = T, recursive = T, pattern = "*_quality_freq.csv")
qvlenfiles <- c(pbfiles, ontfiles)

qvlen <- NULL
qvlen <- list()
for (i in 1:length(qvlenfiles)) {
  x <- read_delim(qvlenfiles[i], delim = ",")
  qvlen[[i]] <- x  
}
qvlen <- bind_rows(qvlen)

qvlen <- qvlen %>% 
  separate(Sample, into = c("sample", "flowcell", "platform"), sep = "_", remove = FALSE) %>%
  dplyr::rename(readlength = Read_Length, readnumbers = Read_Numbers, qv = QV) %>%
  mutate(platform = case_when(platform == "PB" ~ "PacBio HiFi", TRUE ~ platform))

### contour plot for QV vs Read Lengths vs Read numbers
qvcountourreads <-  qvlen %>% 
  group_by(platform, readlength, qv) %>% 
  summarise(readnumbers = sum(readnumbers)) %>%
  group_by(platform) %>%
  mutate(readlength = readlength + 50, readfraction = readnumbers/sum(readnumbers)) %>% 
  ggplot(aes(x=readlength, y = qv, z = readfraction)) + 
  geom_contour_filled(alpha=1) +
  facet_wrap(~platform, scales = "free") +
  theme_bw() + 
  guides(fill = guide_legend(title = "Read Fraction")) +
  ylab("Average Read QV") + xlab("Read Length (bp)") +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(hjust=1,angle=45, size = 16), legend.position = "bottom")

qvcountourreads
## contour plot for QV vs read length vs base counts
qvcontourbases <- qvlen %>% 
  group_by(platform, readlength, qv) %>% 
  summarise(basenumbers = sum(readnumbers*(readlength+50))) %>%
  group_by(platform) %>%
  mutate(readlength = readlength + 50, basefraction = basenumbers/sum(basenumbers)) %>%
  ggplot(aes(x=readlength, y = qv, z = basefraction)) + 
  geom_contour_filled(alpha=1) +
  facet_wrap(~platform, scales = "free") +
  theme_bw() + 
  guides(fill = guide_legend(title = "Base Fraction")) +
  ylab("Average Read QV") + xlab("Read Length (bp)") +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(hjust=1,angle=45, size = 16), legend.position = "bottom")

qvcontourbases
##plot QV readnumbers
## for ONT, only few reads have >30QV. messes up the plot so we mutate >30 QV to 30 QV for ONT
qvlenplot <- qvlen %>% 
  group_by(platform, readlength, qv) %>% 
  summarise(readnumbers = sum(readnumbers)) %>% 
  mutate(readlength = readlength + 50) %>% 
  mutate(basecounts = readlength * readnumbers) %>% 
  mutate(qv = case_when(qv > 30 & platform == "ONT" ~ 30, TRUE ~ qv)) %>%
  ggplot(aes(x=qv,y=readnumbers)) + 
  geom_bar(stat = "identity") +
  xlab("Read average QV") + ylab("No. of reads") +
  facet_wrap(~platform, scales = "free")
qvlenplot
####plot QV
setwd(projectdir)
ggsave(paste(plotdir, "qvvsreadfraccontour.pdf", sep = ""), qvcountourreads, width = 8, height = 6)
ggsave(paste(plotdir, "qvvsbasefraccontour.pdf", sep = ""), qvcontourbases, width = 8, height = 6)
ggsave(paste(plotdir, "qvvsreadnumhistogram.pdf", sep = ""), qvlenplot, width = 8, height = 6)

#### need to change legend size and perhaps the ranges as well
#### perhaps to add some stats out of this in the future around read QV and number/proportion of bases at certain QV

### read in the stats.csv file to get the number of reads etc and generate supplementary and main report table excel
### metrics calculated for the long and short reads are different so have to process them separately
### we can restrict to only the total number of bases, number of reads and average read-length to report consistently across technology
### these metrics are used for calculating genome size and only the average read length and number of bases are required for the genome size calculations
### QV metrics can perhaps display median and N50 values of reads so we dont have to do it here again
### This code is very messy. Perhaps upstream reporting can be improved to make this easy

##tech <- c("PB", "ONT", "illumina")

setwd(projectdir)
rawqcdir <- "./analysis/raweval/"
setwd(rawqcdir)

pbfiles <- list.files("./PB", "*_stats.csv", full.names = T)
pbmetrics <- NULL
pbmetrics <- list()
for (i in 1:length(pbfiles)) {
  x <- read_delim(pbfiles[i], delim = ",")
  x$tech <- str_split_i(pbfiles[i], "/", 2)
  pbmetrics[[i]] <- x
}
pbmetrics <- bind_rows(pbmetrics)

ontfiles <- list.files("./ONT", "*_stats.csv", full.names = T)
ontmetrics <- NULL
ontmetrics <- list()
for (i in 1:length(ontfiles)) {
  x <- read_delim(ontfiles[i], delim = ",")
  x$tech <- str_split_i(ontfiles[i], "/", 2)
  ontmetrics[[i]] <- x
}
ontmetrics <- bind_rows(ontmetrics)

ilmnfiles <- list.files("./illumina", "*_QC.csv", full.names = T, recursive = T)
ilmnmetrics <- NULL
ilmnmetrics <- list()
for (i in 1:length(ilmnfiles)){
  x <- read_delim(ilmnfiles[i], delim = ",")
  x <- x %>% separate(Number_of_reads, c("r1nreads","r2nreads"), sep = ";") %>%
    separate(Number_of_bases, c("r1nbases","r2nbases"), sep = ";") %>% 
    select(r1nreads,r2nreads,r1nbases,r2nbases) %>% 
    mutate_all(as.numeric) %>%
    mutate(nreads = r1nreads + r2nreads, nbases = r1nbases + r2nbases) 
  x$fid <- i
  tech <- str_split_i(ilmnfiles[i], "/", 3)
  if (tech == "dnaseq") {
    tech <- "Illumina"
  }
  else if (tech == "hic") {
    tech <- "HiC"
  }
  x$tech <- tech
  ilmnmetrics[[i]] <- x
}

ilmnmetrics <- bind_rows(ilmnmetrics)

###collect read numbers, total bases per tech
ilmnmetrics <- ilmnmetrics %>% 
  select(tech, nreads, nbases) %>% 
  group_by(tech) %>%
  reframe(tech, nreads = sum(nreads), nbases = sum(nbases)) %>% 
  distinct() %>% 
  mutate(meanrl = nbases/nreads)

lrmetrics <- bind_rows(pbmetrics, ontmetrics)
lrmetrics <- lrmetrics %>% 
  group_by(tech) %>% 
  reframe(tech, nreads = sum(`Total_Reads`), nbases = sum(`Total_Bases`)) %>%
  distinct() %>% 
  mutate(meanrl = nbases/nreads)

### borrow info from peaks data to calculate genome size
selectk <- 17
readmetrics <- NULL
readmetrics <- bind_rows(lrmetrics, ilmnmetrics)
readmetrics <- left_join(readmetrics, 
                         peaks %>% 
                           filter(klen == selectk) %>% 
                           select(tech, mode = kcounts, k = klen) %>% 
                           mutate(k = as.numeric(k), mode = mode + 0.5))

readmetrics <- readmetrics %>% mutate(meanrd = mode * (meanrl / (meanrl - k + 1 )), genomesize = nbases / meanrd)

readmetrics <- readmetrics %>% 
  select("Sequencing Platform" = tech, 
         "Read count" = nreads, 
         "Base count" = nbases, 
         "Mean read length (bp)" = meanrl, 
         "Kmer mode" = mode,
         "K-length" = k,
         "Avg. read depth" = meanrd,
         "Genome size estimate" = genomesize)
##TODO: 
## round the numbers reasonably, i.e. number of bases can be in GB, mean values to be rounded to nearest integer etc
setwd(projectdir)
write_delim(readmetrics, paste(plotdir, "ReadMetrics.tsv", sep = ""), delim = "\t")

###########################################################################################
##TODO: Fasta file summary metrics

##read fasta index file to start the table
setwd(projectdir)
faidx <- paste(reffasta, ".fai", sep = "")
seqinfo <- read_delim(faidx, delim = "\t", col_names = F)
colnames(seqinfo) <- c("scaffoldid", "length", "offset", "width", "boffset")
seqinfo <- seqinfo %>% rownames_to_column(var = "seqorder")

##read depth file for obtaining those metrics
## read depth estimates for chromosome labels
##TODO: change the path later
rdfiles <- list.files("analysis/depth/rearranged/", full.names = T, recursive = T, pattern = "*.10000.depth.bed$")

mappeddepth <- NULL
mappeddepth <- list()
for (i in 1:length(rdfiles)) {
  x <- read_delim(rdfiles[i], col_names = F, delim = " ")
  x$platform <- str_extract(rdfiles[i], "\\S+\\.(\\S+)\\.10000\\.depth\\.bed$", group = 1)
  mappeddepth[[i]] <- x
}
mappeddepth <- bind_rows(mappeddepth)
colnames(mappeddepth) <- c("scaffoldid", "start", "end", "depth", "platform")

###attach mapped read-depth information to seqinfo
## this is used to identify sequences with half the median depth for sex chromosomes
seqinfo <- left_join(seqinfo, mappeddepth %>% group_by(scaffoldid, platform) %>% summarise(mediandepth=median(depth)) %>% pivot_wider(names_from = "platform", values_from = "mediandepth")) %>% View()

##mark sequences that are mostly ribosomal DNA. we can just add the number of bases that are rDNA bases
##TODO: make this into not so hardcoded path
rdna <- read_delim("./analysis/rDNA/POGVIT.v2.1.rdnaregions.bed", col_names = F, delim = "\t")
colnames(rdna) <- c("scaffoldid", "start", "end")
rdna <- rdna %>% mutate(rdnabases = end - start)

###join with the seqinfo
seqinfo <- left_join(seqinfo, rdna %>% group_by(scaffoldid) %>% summarise(rdnabases = sum(rdnabases)))
seqinfo %>% mutate(scftype = case_when(is.na(rdnabases) ~ "Other", !is.na(rdnabases) & rdnabases/length>0.9 ~ "rDNA", TRUE ~ "Other"))

##mark sequences that are mostly satellite repeats


rd <- left_join(rd, fai)

rd %>% 
  filter(scfid == "BASDUscf7.2") %>% 
  group_by(platform) %>% 
  summarise(d = median(depth))

###########################################################################################
y <- read_delim("/g/data/xl04/genomeprojects/Bassiana_duperreyi/AusARG_00001.1_BASDUv1.0_pseudohap_asmreport.txt", delim = "\t")
y %>% 
  filter(length > 1e6 & chrlabel != "putY") %>% 
  mutate(chrtype = case_when(order < 6 ~ "macro", order > 6 & order < 15 ~ "micro", TRUE ~ chrlabel)) %>%
  ggplot(aes(x=mediangc, y = length, label = scfid)) + 
  geom_point(size = 2) + 
  geom_label_repel(box.padding = 0.4, aes(fill = chrtype)) + 
  theme_bw() + 
  xlab("Median GC%") + 
  ylab("Scaffold Length (Mbp)") + 
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6, sep = "")) + 
  theme(text = element_text(size = 16))

#########################
##Y-chromosome labels
y %>% filter(! chrlabel %in% c("mtGenome","rDNA", "censat")) %>% 
  mutate(chrtype = case_when(order < 6 ~ "macro", order > 6 & order < 15 ~ "micro", TRUE ~ chrlabel)) %>%
  ggplot(aes(x=chrtype, y = illum)) + geom_boxplot()

#######
##repeat analyses
repeats <- read_delim("/g/data/xl04/genomeprojects/Bassiana_duperreyi/repeatmasker.BASDU-out.gff", delim = "\t", col_names = F, comment = "#")
colnames(repeats) <- c("scfid", "program", "feature", "start", "end", "score", "strand", "x", "info")

repeats <- repeats %>% mutate(info = str_remove(info, "Target \"Motif:"), info = str_remove(info, "\".*")) %>% 
  mutate(class = case_when(str_detect(info, "^\\(\\S+\\)") ~ "Tandem Repeat", 
                           str_detect(info, "-rich") ~ "Tandem Repeat",
                           str_detect(info, "polypurine") ~ "Tandem Repeat",
                           str_detect(info, "^rnd") ~ NA,
                           TRUE ~ NA))

rmlabel <- readDNAStringSet("/g/data/xl04/genomeprojects/Bassiana_duperreyi/repeatmodeler.BASDU-families.fa")
rmlabel <- data.frame(fn = names(rmlabel))
rmlabel <- rmlabel %>% mutate(info = str_extract(fn, "(rnd-\\d+_family-\\d+)#(\\S+)", group = 1), 
                   repeatname = str_extract(fn, "(rnd-\\d+_family-\\d+)#(\\S+)", group = 2),
                   class = str_split_i(repeatname, "/", i = 1), 
                   subclass = str_split_i(repeatname, "/", i = 2),
                   repeatid = str_extract(fn, "(rnd-\\d+_family-\\d+#\\S+)", group = 1))
rmlen <- read_delim("/g/data/xl04/genomeprojects/Bassiana_duperreyi/repeatmodeler.BASDU-families.fa.fai", delim = "\t", col_names = F)
colnames(rmlen) <- c("repeatid", "length", "offset", "width", "boffset")
rmlabel <- left_join(rmlabel, select(rmlen, repeatid, replength = length))

complexrepeats <- repeats %>% filter(is.na(class)) %>% select(-class)
complexrepeats <- left_join(complexrepeats, select(rmlabel, -fn))

simplerepeats <- repeats %>% filter(!is.na(class))
repeats <- bind_rows(complexrepeats, simplerepeats)
repeats <- repeats %>% mutate(class = case_when(class == "Simple_repeat" ~ "Simple Repeat", 
                                                class == "Tandem Repeat" ~ "Simple Repeat", TRUE ~ class))
repeats %>% mutate(l = end - start + 1) %>%
  group_by(class) %>%
  summarise(total = sum(l)) %>% 
  mutate(total = total * 100 / 1567894183) %>%
  ggplot(aes(x="", y=total, fill=class)) + 
  geom_bar(stat="identity", width=1) + 
  coord_polar("y", start=0) + 
  scale_fill_manual(values = c("#004949FF", "#009292FF", "#FF6DB6FF", "#FFB6DBFF", "#490092FF", "#006DDBFF", "#B66DFFFF", "#6DB6FFFF", "#B6DBFFFF", "#920000FF", "#924900FF", "#DB6D00FF"), 
                    name = "Repeat types") + 
  theme_void()

###check unknown overlaps with the exons
unknownrepeats <- filter(repeats, repeatname == "Unknown")
p <- unknownrepeats %>% mutate(width = end - start + 1, 
                          replengroup = case_when(replength < 2000 ~ floor(replength/100) * 100, 
                                                  replength >= 2000 ~ floor(replength/1000) * 1000)) %>% 
  group_by(replengroup) %>% 
  summarise(tw = sum(width)) %>% 
  mutate(rg = paste(replengroup, "-", replengroup + 99, sep = ""), 
         rg = case_when(replengroup == 2000 ~ "2000-2999", 
                        replengroup == 3000 ~ "3000-3999", 
                        replengroup == 4000 ~ "4000-4999", 
                        TRUE ~ rg)) %>% 
  arrange(replengroup) %>% 
  mutate(rg = factor(rg, levels=rg)) %>% 
  ggplot(aes(x=rg,y=tw)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), text = element_text(size=16)) +
  xlab("Repeat unit length (bp)") + 
  ylab("Summed length (Mbp)") + 
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6, sep = "", suffix = ""))


###
genes <- read_delim("/g/data/xl04/jc4878/Bassiana_publicationV2/genomeannotation/Augustus_annotation.gff3", delim = "\t", col_names = F, comment = "#")
colnames(genes) <- c("scfid", "program", "feature", "start", "end", "score", "strand", "x", "info")

exons <- genes %>% filter(feature %in% c("CDS", "five_prime_utr", "three_prime_utr"))
eranges <- GRanges(seqnames = pull(exons, scfid), ranges = IRanges(start = pull(exons, start), end = pull(exons, end)), strand = pull(exons, strand))
rranges <- GRanges(seqnames = pull(unknownrepeats, scfid), ranges = IRanges(start = pull(unknownrepeats, start), end = pull(unknownrepeats, end)), strand = pull(unknownrepeats, strand))
hits <- findOverlaps(rranges, eranges, minoverlap = 50)
overlaps <- pintersect(rranges[queryHits(hits)], eranges[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(rranges[queryHits(hits)])


###########
###Merqury redrawing
###
setwd(projectdir)
spectrafile <- list.files("analysis/merqury", "*.spectra-asm.hist", full.names = T)
mer <- read_delim(spectrafile, delim="\t")
merquryplot <- mer %>% 
  filter(kmer_multiplicity<150 & Count<2e7) %>% 
  mutate(Assembly = case_when(str_detect(Assembly, ".h1.") ~ "hap1-specific", 
                              str_detect(Assembly, ".h2.") ~ "hap2-specific", 
                              Assembly == "read-only" ~ "missing-in-asm", 
                              Assembly == "shared" ~ "found-in-both", 
                              TRUE ~ Assembly)) %>% 
  ggplot(aes(x=kmer_multiplicity, y = Count, fill = Assembly)) + 
  geom_area(alpha=0.4) + 
  theme_bw() +
  theme(text = element_text(size=16))

ggsave(paste(plotdir, "merqury-kmer-histogram.pdf", sep = ""), merquryplot, width = 8, height = 6)

#### Process SRF output to detect centromeric satellite repeats
a <- read_delim("analysis/SRF/POGVIT.v2.1.asm2srf.genome.bed", delim = "\t", col_names = F)
colnames(a) <- c("chr", "start", "end", "srid", "abun", "seqlen", "replen", "no")
a <- a %>% separate(srid, c("sampleid", "repid"), sep = "#", remove = F)
a <- a %>% mutate(width = end - start, rcopy = width / replen)
##get unique repeat IDs to iterate over it to find where they may and what part of the genome they cover
repids <- pull(a, repid) %>% unique()

repreduced <- NULL
repreduced <- list()
for (i in 1:length(repids)){
  gr <- makeGRangesFromDataFrame(a %>% filter(repid==repids[i]) %>% select(1:3))
  gr <- data.frame(reduce(gr, min.gapwidth = 10))
  gr$repid <- repids[i]
  repreduced[[i]] <- gr
}
repreduced <- bind_rows(repreduced)



a <- read_delim("analysis/SRF/POGVITdef.h1.corrected_YAHS.asm2srf.genome.bed", delim = "\t", col_names = F)
colnames(a) <- c("chr", "start", "end", "srid", "abun", "seqlen", "replen", "no")
a
a <- read_delim("analysis/SRF/POGVIT.v2.1.asm2srf.genome.bed", delim = "\t", col_names = F)
####process SRF output
rdna <- read_delim("/g/data/xl04/genomeprojects/Pogona_vitticeps/analysis/rDNA/POGVIT.v2.1.rdnaregions.bed", col_names = F, delim = "\t")
colnames(rdna) <- c("chr", "start", "end")
rdna <- makeGRangesFromDataFrame(rdna)

# read SRF bed file
srf <- read_delim("/g/data/xl04/genomeprojects/Pogona_vitticeps/analysis/SRF/POGVIT.v2.1.asm2srf.genome.bed", delim = "\t", col_names = F)
colnames(srf) <- c("chr", "start", "end", "srid", "abun", "seqlen", "replen", "no")
srf <- srf %>% separate(srid, c("asmid", "repid"), sep = "#")
repids <- pull(srf, repid) %>% unique()

##merge consecutive repeats of the same type into a window
repreduced <- NULL
repreduced <- list()
for (i in 1:length(repids)){
  gr <- makeGRangesFromDataFrame(srf %>% filter(repid==repids[i]) %>% select(1:3))
  gr <- data.frame(reduce(gr, min.gapwidth = 10))
  gr$repid <- repids[i]
  repreduced[[i]] <- gr
}
repreduced <- bind_rows(repreduced)
##join sequence length and repeat length
repreduced <- left_join(repreduced, srf %>% select(chr, seqlen) %>% distinct(), by = c("seqnames" = "chr"))
repreduced <- left_join(repreduced, select(srf, repid, replen) %>% distinct(), by = c("repid"))


repreduced %>% 
  mutate(rcopy = width/replen) %>% 
  filter(seqlen>10e6 & replen > 6 & rcopy > 25) %>% 
  ggplot(aes(x=start,y=rcopy, label = repid, shape = as.character(replen))) + 
  geom_point() + 
  geom_text_repel() +
  facet_wrap(~seqnames, scales = "free") +
  theme_bw()


#######srf processing from the asm 2 enlong.fa (paf)
srfpaf <- read_paf("analysis/SRF/POGVIT.v2.1.asm2srf.genome.paf", tibble = T)
srfpaf <- srfpaf %>% separate(tname, c("sampleid", "repid"), sep = "#", remove = F) %>% mutate(replen = as.numeric(str_extract(repid, "\\d+$")))
repids <- srfpaf %>% pull(tname) %>% unique()
##merge consecutive repeats of the same type into a window
repreduced <- NULL
repreduced <- list()
for (i in 1:length(repids)){
  repaln <- srfpaf %>% filter(tname==repids[i]) %>% select(chr = qname, start = qstart, end = qend, de)
  repaln <- makeGRangesFromDataFrame(repaln, keep.extra.columns = T)
  ##merge two consecutive alignments as one if there is <=10bp gap between them, this can be increased to 100bp or corresponding to the gap width of scaffolds
  ##merge can be based on length of repeats as well, i.e. if two consecutive repeats
  repregions <- reduce(repaln, min.gapwidth = 10)
  ###get overlapping alignments to collect metrics of how many alignments overlap the region and the mean de (divergence) score for all alignments
  olap <- findOverlapPairs(repregions, repaln)
  olap <- data.frame(olap) %>% group_by(first.seqnames, first.start, first.end, first.width) %>% summarise(alncounts = n(), meande = mean(second.de))
  olap <- ungroup(olap)
  colnames(olap) <- c("chr","start", "end", "width", "alncounts", "meande")
  olap$repid <- repids[i]
  olap$replen <- as.numeric(str_extract(repids[i], "\\d+$"))
  olap <- olap %>% mutate(rangeid = paste(chr, paste(start,end,sep = "-"), sep = ":"))
  repreduced[[i]] <- olap
}
repreduced <- bind_rows(repreduced)
print(paste(nrow(repreduced), " ranges before filtering", sep = ""))
##tried filtering the raw alignments. however, this can loose partial alignments on the edges of a motif
##so best to group them first and them remove them if they are too short
##remove alignments less than 100bp long since enlong extends repeat units to min 200bp
##remove alignments that are less than 10% of the repeat length to remove noisy alignments of long repeats
#& alen > replen * 0.1 & qend - qstart > 100
##get a copy estimate based on the width
repreduced <- repreduced %>% filter(width > 100 & width > replen * 0.1) %>% mutate(rcopy = round(width/replen, digits = 0))
print(paste(nrow(repreduced), " ranges after filtering", sep = ""))
##find overlapping alignments of different repeats
repranges <- makeGRangesFromDataFrame(repreduced, keep.extra.columns = T)
repolap <- data.frame(findOverlapPairs(repranges, repranges, minoverlap=10))
repolap <- repolap %>% group_by(first.rangeid, first.repid) %>% mutate(olapcount = n())
repolap <- repolap %>% ungroup()
##if there is only one overlapping repeat alignments, they are good to go as they will be self overlap
selection1 <- repolap %>%
  filter(olapcount == 1) %>%
  select(chr = first.X.seqnames,
         start = first.X.start,
         end = first.X.end,
         width = first.X.width,
         repid = first.repid,
         replen = first.replen,
         alncounts=first.alncounts,
         meande = first.meande,
         rcopy = first.rcopy,
         rangeid = first.rangeid,
         olapcount)
print(paste(pull(selection1,rangeid) %>% unique() %>% length(), " ranges with no other overlapping ranges", sep = ""))
### for overlapping repeats, discard the one with shorter overall width
### logic is that smaller alignments overall within a larger alignment are likely to be spurious hits
### also discard logic says that there is a bigger alignment compared to the discarded one. The biggest will never be discarded
discard1 <- repolap %>% filter(olapcount > 1) %>% filter(first.X.repid != second.X.repid & 
                                                           first.X.rangeid != second.X.rangeid & 
                                                           first.X.width > second.X.width &
                                                           first.X.start <= second.X.start &
                                                           first.X.end >= second.X.end
                                                         ) %>%
  pull(second.rangeid)
discard2 <- repolap %>% filter(olapcount > 1) %>% filter(first.X.repid != second.X.repid & 
                                                           first.X.rangeid != second.X.rangeid & 
                                                           first.X.width < second.X.width &
                                                           first.X.start >= second.X.start &
                                                           first.X.end <= second.X.end) %>%
  pull(first.rangeid)
discard <- c(discard1, discard2)
discard <- unique(discard)
print(paste(length(discard), " ranges discarded in favour of the larger ones overlapping them", sep = ""))
selection2 <- repolap %>% filter(olapcount > 1 & ! first.rangeid %in% discard) %>% select(chr = first.X.seqnames,
                                                                                          start = first.X.start,
                                                                                          end = first.X.end,
                                                                                          width = first.X.width,
                                                                                          repid = first.repid,
                                                                                          replen = first.replen,
                                                                                          alncounts=first.alncounts,
                                                                                          meande = first.meande,
                                                                                          rcopy = first.rcopy,
                                                                                          rangeid = first.rangeid,
                                                                                          olapcount)
# selection2 needs to be of unique ranges. redundancy removed by selecting alignments with maximum overlap count, minimum row number and smallest meande
selection2 <- selection2 %>% group_by(rangeid) %>% filter(olapcount == max(olapcount) & meande == min(meande)) %>% filter(row_number()==1)
print(paste(nrow(selection2), " ranges selected in favour of the shorter/less optimal ones overlapping them", sep = ""))
selection <- bind_rows(selection1, selection2) %>% arrange(chr, start)
###
chrlen <- srfpaf %>% select(chr=qname, chrlen = qlen) %>% distinct()
###
selection <- left_join(selection, chrlen)
###
##extra filters can be applied here to focus on truely repetetive regions
##srf can capture repeats that are scattered and dont fold on themselves
write_delim(selection, "analysis/SRF/POGVIT.v2.1.asm2srf.genome.rscript.txt", delim = "\t")


#####merge the class together
srfpaf <- read_paf("analysis/SRF/POGVIT.v2.1.asm2srf.genome.paf", tibble = T)
srfpaf <- srfpaf %>% separate(tname, c("sampleid", "repid"), sep = "#", remove = F) %>% mutate(replen = as.numeric(str_extract(repid, "\\d+$")))

classpaf <- read_paf("analysis/SRF/POGVIT.v2.1.enlong1K.paf", tibble = T)
reppairs <- classpaf %>% select(from = qname, to = tname) %>% distinct()
repgraph <- tbl_graph(nodes = NULL, edges = reppairs, directed = FALSE)
repclass <- repgraph %>% mutate(component = group_components()) %>% as_tibble()
repclass <- repclass %>% mutate(repclass = paste("srfclass", component, sep = "-"))

srfpaf <- left_join(srfpaf, repclass, by = c("tname" = "name"))
repclasses <- srfpaf %>% pull(repclass) %>% unique()

##merge consecutive repeats of the same type into a window
repreduced <- NULL
repreduced <- list()
for (i in 1:length(repclasses)){
  repaln <- srfpaf %>% filter(repclass==repclasses[i]) %>% select(chr = qname, start = qstart, end = qend, de)
  repaln <- makeGRangesFromDataFrame(repaln, keep.extra.columns = T)
  ##merge two consecutive alignments as one if there is <=10bp gap between them, this can be increased to 100bp or corresponding to the gap width of scaffolds
  ##merge can be based on length of repeats as well, i.e. if two consecutive repeats
  repregions <- reduce(repaln, min.gapwidth = 10)
  ###get overlapping alignments to collect metrics of how many alignments overlap the region and the mean de (divergence) score for all alignments
  olap <- findOverlapPairs(repregions, repaln)
  olap <- data.frame(olap) %>% group_by(first.seqnames, first.start, first.end, first.width) %>% summarise(alncounts = n(), meande = mean(second.de))
  olap <- ungroup(olap)
  colnames(olap) <- c("chr","start", "end", "width", "alncounts", "meande")
  olap$repid <- repclasses[i]
  tmprid <- srfpaf %>% filter(repclass==repclasses[i]) %>% pull(tname) %>% unique() %>% sort()
  olap$replen <- as.numeric(str_extract(tmprid[1], "\\d+$"))
  olap <- olap %>% mutate(rangeid = paste(chr, paste(start,end,sep = "-"), sep = ":"))
  repreduced[[i]] <- olap
}
repreduced <- bind_rows(repreduced)
print(paste(nrow(repreduced), " ranges before filtering", sep = ""))

##find overlapping alignments of different repeats
repranges <- makeGRangesFromDataFrame(repreduced, keep.extra.columns = T)
repolap <- data.frame(findOverlapPairs(repranges, repranges, minoverlap=10))
repolap <- repolap %>% group_by(first.rangeid, first.repid) %>% mutate(olapcount = n())
repolap <- repolap %>% ungroup()
##if there is only one overlapping repeat alignments, they are good to go as they will be self overlap
selection1 <- repolap %>%
  filter(olapcount == 1) %>%
  select(chr = first.X.seqnames,
         start = first.X.start,
         end = first.X.end,
         width = first.X.width,
         repid = first.repid,
         replen = first.replen,
         alncounts=first.alncounts,
         meande = first.meande,
         rcopy = first.rcopy,
         rangeid = first.rangeid,
         olapcount)
print(paste(pull(selection1,rangeid) %>% unique() %>% length(), " ranges with no other overlapping ranges", sep = ""))
### for overlapping repeats, discard the one with shorter overall width
### logic is that smaller alignments overall within a larger alignment are likely to be spurious hits
### also discard logic says that there is a bigger alignment compared to the discarded one. The biggest will never be discarded
discard1 <- repolap %>% filter(olapcount > 1) %>% filter(first.X.repid != second.X.repid & 
                                                           first.X.rangeid != second.X.rangeid & 
                                                           first.X.width > second.X.width &
                                                           first.X.start <= second.X.start &
                                                           first.X.end >= second.X.end
) %>%
  pull(second.rangeid)
discard2 <- repolap %>% filter(olapcount > 1) %>% filter(first.X.repid != second.X.repid & 
                                                           first.X.rangeid != second.X.rangeid & 
                                                           first.X.width < second.X.width &
                                                           first.X.start >= second.X.start &
                                                           first.X.end <= second.X.end) %>%
  pull(first.rangeid)
discard <- c(discard1, discard2)
discard <- unique(discard)
print(paste(length(discard), " ranges discarded in favour of the larger ones overlapping them", sep = ""))
selection2 <- repolap %>% filter(olapcount > 1 & ! first.rangeid %in% discard) %>% select(chr = first.X.seqnames,
                                                                                          start = first.X.start,
                                                                                          end = first.X.end,
                                                                                          width = first.X.width,
                                                                                          repid = first.repid,
                                                                                          replen = first.replen,
                                                                                          alncounts=first.alncounts,
                                                                                          meande = first.meande,
                                                                                          rcopy = first.rcopy,
                                                                                          rangeid = first.rangeid,
                                                                                          olapcount)
# selection2 needs to be of unique ranges. redundancy removed by selecting alignments with maximum overlap count, minimum row number and smallest meande
selection2 <- selection2 %>% group_by(rangeid) %>% filter(olapcount == max(olapcount) & meande == min(meande)) %>% filter(row_number()==1)
print(paste(nrow(selection2), " ranges selected in favour of the shorter/less optimal ones overlapping them", sep = ""))
selection <- bind_rows(selection1, selection2) %>% arrange(chr, start)
###
chrlen <- srfpaf %>% select(chr=qname, chrlen = qlen) %>% distinct()
###
selection <- left_join(selection, chrlen)
###
##extra filters can be applied here to focus on truely repetetive regions
##srf can capture repeats that are scattered and dont fold on themselves
write_delim(selection, "analysis/SRF/POGVIT.v2.1.asm2srf.genome.rscript.txt", delim = "\t")




genome <- toGRanges(data.frame(mutate(selection, start = 0) %>% filter(chrlen>100e6) %>% select(chr, start, end = chrlen) %>% distinct() %>% arrange(desc(end))))
dev.off()
kp <- plotKaryotype(genome, chromosomes = "all")
satrep <- makeGRangesFromDataFrame(data.frame(selection %>% 
                                                filter(chr %in% seqnames(genome)) %>%
                                                filter(replen > 6 & replen < 1000 & (rcopy > 25 | width >10000)) %>% 
                                                select(chr, start, end, repid, y = width, rcopy)), 
                                   keep.extra.columns = T)
kpPlotMarkers(kp, satrep, labels = satrep$repid)
#kpPlotRegions(kp, satrep)
kpPoints(kp, data = satrep)

dev.off()

selection %>% filter(rcopy>10 & replen %in% c(151,115,230,300,368))
selection %>% filter(rcopy>10 & replen %in% c(151,115,230,300,368)) %>% View()
selection %>% filter(rcopy>10 & replen %in% c(151,115,230,300,368)) %>% pull(chr) %>% unique() %>% length()
selection %>% filter(replen > 6 & replen < 1000 & (rcopy > 25 | width >25000)) %>% pull(chr) %>% unique() %>% length()
selection %>% filter(replen > 6 & replen < 1000 & (rcopy > 25 | width >25000)) %>% View()
write_delim(selection %>% filter(replen > 6 & replen < 1000 & (rcopy > 25 | width >25000)) %>% select(1,4,2,3)
)
selection %>% filter(replen > 6 & replen < 1000 & (rcopy > 25 | width >25000)) %>% select(1,4,2,3)
write_delim(selection %>% filter(replen > 6 & replen < 1000 & (rcopy > 25 | width >25000)) %>% select(1,4,2,3), "/g/data/xl04/genomeprojects/Pogona_vitticeps/analysis/SRF/putative.centro.csv", delim = ",")
