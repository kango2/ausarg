.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))
library(tidyverse)
library(karyoploteR)

setwd("/g/data/xl04/ka6418/chromosome_graph/")

minseqlen <- 5e5

mygenome <- read_delim("rTilRug_HiC_pctg_seqtable.csv", delim = ",") 
mygenomegr <- makeGRangesFromDataFrame(data.frame(mutate(mygenome, start = 1) %>% 
                                                    filter(Length >= minseqlen) %>%
                                                    arrange(desc(Length)) %>%
                                                    select(chr = `Sequence ID`, start, end = Length)
                                                  ))
telomeresdf <- read_delim("rTilRug_HiC_pctg_Telomeres.csv", delim = ",") 
padtlen <- 1e5
telomeresgr <- toGRanges(data.frame(select(telomeresdf, chr=Sequence_ID, start = Start, end = End))) 

cendf <- read_delim("TRASH_output/Summary.of.repetitive.regions.rTilRug_HiC_pctg.fasta.csv", delim = ",") 
padtlen <- 1e5
cengr <- toGRanges(data.frame(select(cendf, chr=name, start = start, end = end))) 

ilmnrd <- read_delim("rTilRug_HiC_pctg_Illumina_sorted.bam.binned.fixed.depth.csv", delim = ",")
ilmnrdgr <- toGRanges(data.frame(mutate(ilmnrd, AverageDepth = case_when(AverageDepth > 200 ~ 200, TRUE ~ AverageDepth)) %>%
                               select(chr = Chromosome, 
                                      start = Start, 
                                      end = End, 
                                      y = AverageDepth)
                             ))

ontrd <- read_delim("rTilRug_HiC_pctg_ONT_sorted.bam.binned.fixed.depth.csv", delim = ",")
ontgr <- toGRanges(data.frame(mutate(ilmnrd, AverageDepth = case_when(AverageDepth > 200 ~ 200, TRUE ~ AverageDepth)) %>%
                                   select(chr = Chromosome, 
                                          start = Start, 
                                          end = End, 
                                          y = AverageDepth)
))

hifird <- read_delim("rTilRug_HiC_pctg_PacBio_sorted.bam.binned.fixed.depth.csv", delim = ",")
hifigr <- toGRanges(data.frame(mutate(ilmnrd, AverageDepth = case_when(AverageDepth > 200 ~ 200, TRUE ~ AverageDepth)) %>%
                                select(chr = Chromosome, 
                                       start = Start, 
                                       end = End, 
                                       y = AverageDepth)
))

gc <- read_delim("rTilRug_HiC_pctg.fasta_GC.csv", delim = ",") 
gcgr <- toGRanges(data.frame(gc %>% mutate(y = `GC Count`*100/10000) %>%
                     select(chr = `Sequence Header`, 
                            start = `Position Start`, 
                            end = `Position End`, 
                            y)))

dev.off()
pp <- getDefaultPlotParams(plot.type=1)
pp$data1inmargin <- 20
pp$data1outmargin <- 120


kp <- plotKaryotype(genome = mygenomegr,plot.params = pp, chromosomes="ptg000010l")
#kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=0.8, numticks = 5, col="#666666", cex=0.5)
kpAddBaseNumbers(kp)
##plot telomeres
kpPlotRegions(kp, data=cengr, col="#FF000080", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=telomeresgr, col="orange", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)
#plotting read depths
kpLines(kp, data = ilmnrdgr,
        r0=0, r1=0.25, ymin = 0, ymax = 200,
        col = "#6495ED", lwd = 0.3) 
kpLines(kp, data = ontgr,
        r0=0.25, r1=0.4, ymin = 0, ymax = 200,
        col = "#FA8072", lwd = 0.3)
kpLines(kp, data = hifigr,
        r0=0.45, r1=0.65, ymin = 0, ymax = 200,
        col = "#3CB371", lwd = 0.3)
kpLines(kp, data = gcgr,
        r0=0.65, r1=0.90, ymin = 40, ymax = 60,
        col = "#6A4DA8", lwd = 0.6)




                                                                                                                                                     y = k$`GC Count`/100, ymin = 0, ymax = 200, r0=0.55, r1=1, col = "blue")

