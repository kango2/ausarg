library(tidyverse)
library(scales)
library(ggrepel)

projectdir <- "/g/data/xl04/genomeprojects/Pogona_vitticeps"

setwd(projectdir)

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

### save the k-mer plots
setwd(projectdir)
ggsave(paste(plotdir, "kmerhisto.pdf", sep = ""), kmerplot, width = 8, height = 8)

### need to save peaks in excel or somewhere for genome size estimations

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

setwd(projectdir)
write_delim(readmetrics, paste(plotdir, "ReadMetrics.tsv", sep = ""), delim = "\t")


###########################################################################################
## read depth estimates for chromosome labels
setwd(projectdir)
faifiles <- list.files("./fasta/yahs", "*.fai", full.names = T)

for (i in 1:length(faifiles)){
  fai <- read_delim(faifiles[i], delim = "\t", col_names = F)
  colnames(fai) <- c("scfid", "length", "offset", "width", "boffset")
  asmid <- str_split_i(faifiles[i], "/", 4)
  asmid <- 
}
fai <- read_delim("/g/data/xl04/ka6418/bassiana/publication-v2/fasta/BASDU_SUP_DEEP_hifiasm_yahs_curv1.fasta.fai", delim = "\t", col_names = F)
colnames(fai) <- c("scfid", "length", "offset", "width", "boffset")
rdfiles <- list.files("/g/data/xl04/ka6418/bassiana/publication-v2/eval/depth/curv1/", full.names = T, recursive = T, pattern = "*.10000.depth.bed")

rd <- NULL
rd <- list()
for (i in 1:length(rdfiles)) {
  x <- read_delim(rdfiles[i], col_names = F, delim = " ")
  x$platform <- str_match(rdfiles[i], "BASDU_curv1.merged.(\\S+).10000")[,2]
  rd[[i]] <- x
}
rd <- bind_rows(rd)
colnames(rd) <- c("scfid", "start", "end", "depth", "platform")

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
rm <- read_delim("/g/data/xl04/genomeprojects/Bassiana_duperreyi/repeatmasker.BASDU-out.gff", delim = "\t", col_names = F, comment = "#")
colnames(rm) <- c("scfid", "program", "feature", "start", "end", "score", "strand", "x", "info")

rm <- rm %>% mutate(info = str_remove(info, "Target \"Motif:"), info = str_remove(info, "\".*")) %>% 
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
                   subclass = str_split_i(repeatname, "/", i = 2))

complexrepeats <- rm %>% filter(is.na(class)) %>% select(-class)
simplerepeats <- rm %>% filter(!is.na(class))
complexrepeats <- left_join(complexrepeats, select(rmlabel, -1))
rm <- bind_rows(complexrepeats, simplerepeats)
rm <- rm %>% mutate(class = case_when(class == "Simple_repeat" ~ "Simple Repeat", class == "Tandem Repeat" ~ "Simple Repeat", TRUE ~ class))
rm %>% mutate(l = end - start + 1) %>%
  group_by(class) %>%
  summarise(total = sum(l)) %>% 
  mutate(total = total * 100 / 1567894183) %>%
  ggplot(aes(x="", y=total, fill=class)) + 
  geom_bar(stat="identity", width=1) + 
  coord_polar("y", start=0) + 
  scale_fill_manual(values = c("#004949FF", "#009292FF", "#FF6DB6FF", "#FFB6DBFF", "#490092FF", "#006DDBFF", "#B66DFFFF", "#6DB6FFFF", "#B6DBFFFF", "#920000FF", "#924900FF", "#DB6D00FF"), 
                    name = "Repeat types") + 
  theme_void()

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

a <- read_delim("analysis/SRF/POGVITdef.h1.corrected_YAHS.asm2srf.genome.bed", delim = "\t", col_names = F)
colnames(a) <- c("chr", "start", "end", "srid", "abun", "seqlen", "replen", "no")
a <- a %>% separate(srid, c("asmid", "repid"), sep = "#")
##get unique repeat IDs to iterate over it to find where they may and what part of the genome they cover
repids <- pull(a, repid) %>% unique()




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


