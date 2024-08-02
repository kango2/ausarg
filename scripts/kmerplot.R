library(tidyverse)
library(scales)
library(ggrepel)
setwd("/g/data/xl04/ka6418/bassiana/publication-v2/eval/kmers/")
files <- list.files("./","*.histo")

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

#1061100027
##get the value of the homozygous peak
##logic, multiple the kmer frequency with counts to get density of data, filter k-counts <= 3 as errors
# peaks <- khisto %>% 
#   mutate(m=kcounts*kfrequency) %>% 
#   filter(kcounts>5 & kcounts < 100) %>% 
#   group_by(tech,klen) %>% 
#   filter(m == max(m))  
#   filter(kcounts>10)

peaks <- khisto %>% 
  filter(kcounts < 100) %>% 
  group_by(tech, klen) %>% 
  arrange(desc(kcounts)) %>% 
  mutate(d=kfrequency-lead(kfrequency)) %>% 
  group_by(tech, klen) %>% 
  arrange(desc(kcounts)) %>%
  filter(d>0) %>%
  filter(kcounts == max(kcounts))


khisto %>% 
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

####
########################################################################################################
##------------------------
pbfiles <- list.files("/g/data/xl04/ka6418/bassiana/publication-v2/eval/rawdata-qc/longread/pb/", full.names = T, recursive = T, pattern = "*_quality_freq.csv")
ontfiles <- list.files("/g/data/xl04/ka6418/bassiana/publication-v2/eval/rawdata-qc/longread/ont/", full.names = T, recursive = T, pattern = "*_quality_freq.csv")
qvlenfiles <- c(pbfiles, ontfiles)

qvlen <- NULL
qvlen <- list()

for (i in 1:length(qvlenfiles)) {
  x <- read_delim(qvlenfiles[i], delim = ",")
  qvlen[[i]] <- x  
}
qvlen <- bind_rows(qvlen)

## contour plot for QV vs Read Lengths vs Read numbers
qvlen %>% 
  group_by(platform, readlength, qv) %>% 
  summarise(readnumbers = sum(readnumbers)) %>%
  group_by(platform) %>%
  mutate(readlength = readlength + 50, readfraction = readnumbers/sum(readnumbers)) %>% 
  ggplot(aes(x=readlength, y = qv, z = readfraction)) + 
  geom_contour_filled(alpha=1) +
  facet_wrap(~platform, scales = "free") +
  theme_bw() + 
  scale_fill_brewer(palette = "heat") +
  guides(fill = guide_legend(title = "Read Fraction")) +
  ylab("Average Read QV") + xlab("Read Length (bp)")


## contour plot for QV vs read length vs base counts
qvlen %>% 
  group_by(platform, readlength, qv) %>% 
  mutate(platform = case_when(platform == "pacbio" ~ "PacBio HiFi", TRUE ~ platform)) %>%
  summarise(basenumbers = sum(readnumbers*(readlength+50))) %>%
  group_by(platform) %>%
  mutate(readlength = readlength + 50, basefraction = basenumbers/sum(basenumbers)) %>%
  ggplot(aes(x=readlength, y = qv, z = basefraction)) + 
  geom_contour_filled(alpha=1) +
  facet_wrap(~platform, scales = "free") +
  theme_bw() + 
  guides(fill = guide_legend(title = "Base Fraction")) +
  ylab("Average Read QV") + xlab("Read Length (bp)") +
  theme(text = element_text(size=16), axis.text.x = element_text(hjust=1,angle=45), legend.position = "bottom")

##plot QV readnumbers
qvlen %>% 
  group_by(platform, readlength, qv) %>% 
  summarise(readnumbers = sum(readnumbers)) %>% 
  mutate(readlength = readlength + 50) %>% 
  mutate(basecounts = readlength * readnumbers) %>% 
  mutate(qv = case_when(qv > 30 & platform == "ONT" ~ 30, TRUE ~ qv)) %>%
  ggplot(aes(x=qv,y=readnumbers)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~platform, scales = "free")

###########################################################################################
## read depth estimates for chromosome labels
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
mer <- read_delim("/g/data/xl04/ka6418/bassiana/publication-v2/eval/merqury/curv1/illumina/trimmed/17/haplotype/BASDU_H1H2.spectra-asm.hist", delim="\t")
head(mer)