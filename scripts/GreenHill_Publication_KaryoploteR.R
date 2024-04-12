myLibPath <- "/g/data/xl04/ka6418/temp/Rlibraries"
.libPaths(c(myLibPath, .libPaths()))
.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))
.libPaths(c("/g/data/te53/software/Rpackages/are-4.2.2g/", .libPaths()))
library(karyoploteR)
library(tidyverse)
library(optparse)


mygenome <- read_delim("/g/data/xl04/ka6418/bassiana/publication/eval/curation_greenhill/fixing_curated/eval/seqtable/BASDU_HifiASM_GreenHill_SUP_CurV1.2_seqtable.csv", delim = ",") 

minseqlen <- 4e6

#mygenome$`Sequence ID` <- sub("_.*", "", mygenome$`Sequence ID`)
mygenomegr <- makeGRangesFromDataFrame(data.frame(mutate(mygenome, start = 1) %>% 
                                                    filter(Length >= minseqlen) %>%
                                                    arrange(desc(Length)) %>%
                                                    select(chr = `Sequence ID`, start, end = Length)
))

telomeresdf <- read_delim("/g/data/xl04/ka6418/bassiana/publication/eval/curation_greenhill/fixing_curated/eval/CurV1.2.fixed.sorted.renamed_Telomeres.csv", delim = ",") 

#telomeresdf$Sequence_ID <- sub("_.*", "", telomeresdf$Sequence_ID)
padtlen <- 1e5
telomeresgr <- toGRanges(data.frame(select(telomeresdf, chr=Sequence_ID, start = Start, end = End))) 

cendf <- read_delim("/g/data/xl04/ka6418/bassiana/publication/eval/curation_greenhill/fixing_curated/eval/TRASH/Summary.of.repetitive.regions.BASDU_HifiASM_GreenHill_SUP_CurV1.2.fasta.csv", delim = ",") 
#cendf$name <- sub("_.*", "", cendf$name)
padtlen <- 1e5
filtered_cendf <- filter(cendf, !is.na(class), width > 1e5)
cengr <- toGRanges(data.frame(filter(filtered_cendf, width > 1e5) %>% select(chr=name, start = start, end = end))) 


# Split data based on the class values
list_of_dataframes <- split(filtered_cendf, filtered_cendf$class)

# Convert each dataframe in the list to GRanges and store back in the list
list_of_granges <- lapply(list_of_dataframes, function(df) {
  toGRanges(data.frame(select(df, chr = name, start = start, end = end)))
})

# If you want to store them in separate variables
names(list_of_granges) <- paste0("cengr_", unique(filtered_cendf$class))
list2env(list_of_granges, envir = .GlobalEnv)


column_names <- c("Chromosome", "Start", "End", "AverageDepth")
ilmnrd <- read_delim("/g/data/xl04/ka6418/bassiana/publication/eval/depth/BASDU_HifiASM_GH_SUP_CurV1.2_fixed/BASDU_HifiASM_GH_SUP_CurV1.2.merged.illum.depth.bed", delim = " ",col_names = column_names)
#ilmnrd$Chromosome <- sub("_.*", "", ilmnrd$Chromosome)
ilmnrdgr <- toGRanges(data.frame(mutate(ilmnrd, AverageDepth = case_when(AverageDepth > 100 ~ 100, TRUE ~ AverageDepth)) %>%
                                   select(chr = Chromosome, 
                                          start = Start, 
                                          end = End, 
                                          y = AverageDepth)
))

column_names <- c("Chromosome", "Start", "End", "AverageDepth")
ontrd <- read_delim("/g/data/xl04/ka6418/bassiana/publication/eval/depth/BASDU_HifiASM_GH_SUP_CurV1.2_fixed/BASDU_HifiASM_GH_SUP_CurV1.2.merged.ont.depth.bed", delim = " ",col_names = column_names)
#ontrd$Chromosome <- sub("_.*", "", ontrd$Chromosome)
ontgr <- toGRanges(data.frame(mutate(ontrd, AverageDepth = case_when(AverageDepth > 100 ~ 100, TRUE ~ AverageDepth)) %>%
                                select(chr = Chromosome, 
                                       start = Start, 
                                       end = End, 
                                       y = AverageDepth)
))

hifird <- read_delim("/g/data/xl04/ka6418/bassiana/publication/eval/depth/BASDU_HifiASM_GH_SUP_CurV1.2_fixed/BASDU_HifiASM_GH_SUP_CurV1.2.merged.pb.depth.bed", delim = " ",col_names = column_names)
#hifird$Chromosome <- sub("_.*", "", hifird$Chromosome)
hifigr <- toGRanges(data.frame(mutate(hifird, AverageDepth = case_when(AverageDepth > 60 ~ 60, TRUE ~ AverageDepth)) %>%
                                 select(chr = Chromosome, 
                                        start = Start, 
                                        end = End, 
                                        y = AverageDepth)
))

gc <- read_delim("/g/data/xl04/ka6418/bassiana/publication/eval/curation_greenhill/fixing_curated/eval/GC/BASDU_HifiASM_GreenHill_SUP_CurV1.2.fasta_GC.csv", delim = ",") 
#gc$`Sequence Header` <- sub("_.*", "", gc$`Sequence Header`)
gcgr <- toGRanges(data.frame(gc %>% mutate(y = `GC Count`*100/10000, y = case_when(y > 60 ~ 60, TRUE ~ y), y = case_when(y < 30 ~ 30, TRUE ~ y))  %>%
                               select(chr = `Sequence Header`, 
                                      start = `Position Start`, 
                                      end = `Position End`, 
                                      y)))

column_names <- c("chromosome", "start", "end", "gaplen")
Nregions <- read_delim("/g/data/xl04/ka6418/bassiana/publication/eval/Ns/BASDU_HifiASM_GreenHill_SUP_CurV1.2_Nregions.csv",delim = ",", col_names = column_names)
Nregions_gr <- toGRanges(data.frame(Nregions) %>% select(chr = chromosome, start = start, end = end) )


#dev.off()
pp <- getDefaultPlotParams(plot.type=1)
pp$data1inmargin <- 20
pp$data1outmargin <- 120
pp$leftmargin <- 0.15


kp <- plotKaryotype(genome = mygenomegr,plot.params = pp)
#kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=0.8, numticks = 5, col="#666666", cex=0.5)
kpAddBaseNumbers(kp)
##plot telomeres
kpPlotRegions(kp, data=cengr_CEN187, col="#FF000080", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=cengr_CEN199, col="#80008080", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)
kpPlotRegions(kp, data=telomeresgr, col="orange", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE)

kpPlotRegions(kp, data=Nregions_gr, col="black", r0=-0.1, r1=-0.35,clipping = TRUE,avoid.overlapping=FALSE,lwd = 0.1)

#plotting read depths


kpDataBackground(kp, data.panel = 1,r0=0.0, r1=0.15)
kpAxis(kp,r0=0.0, r1=0.15,cex=0.2,numticks = 2,ymin=0,ymax=100,data.panel = 1,mai=c(0,0,0,1))
kpAddLabels(kp, labels="Illumina", r0=0.0, r1=0.15, data.panel = 1,cex=0.2,label.margin=0.003,side="right")

kpDataBackground(kp, data.panel = 1,r0=0.25, r1=0.40)
kpAxis(kp, r0=0.25, r1=0.40,cex=0.2,numticks = 2,ymin=0,ymax=100,data.panel = 1)
kpAddLabels(kp, labels="ONT", r0=0.25, r1=0.40, data.panel = 1,cex=0.2,label.margin=0.003,side="right")

kpDataBackground(kp, data.panel = 1,r0=0.50, r1=0.65)
kpAxis(kp, r0=0.50, r1=0.65,cex=0.2,numticks = 2,ymin=0,ymax=60,data.panel = 1)
kpAddLabels(kp, labels="HiFi", r0=0.50, r1=0.65, data.panel = 1,cex=0.2,label.margin=0.003,side="right")

kpDataBackground(kp, data.panel = 1,r0=0.75, r1=0.90)
kpAxis(kp, r0=0.75, r1=0.90,cex=0.2,numticks = 2,ymin=30,ymax=60,data.panel = 1)
kpAddLabels(kp, labels="GC", r0=0.75, r1=0.90, data.panel = 1,cex=0.2,label.margin=0.003,side="right")



kpLines(kp, data = ilmnrdgr,
        r0=0.0, r1=0.15, ymin = 0, ymax = 100,
        col = "#6495ED", lwd = 0.1,data.panel = 1) 
kpLines(kp, data = ontgr,
        r0=0.25, r1=0.40, ymin = 0, ymax = 100,
        col = "#FA8072", lwd = 0.1,data.panel = 1)
kpLines(kp, data = hifigr,
        r0=0.50, r1=0.65, ymin = 0, ymax = 60,
        col = "#3CB371", lwd = 0.1,data.panel = 1)
kpLines(kp, data = gcgr,
        r0=0.75, r1=0.90, ymin = 30, ymax = 60,
        col = "#6A4DA8", lwd = 0.1,data.panel = 1)




#y = k$`GC Count`/100, ymin = 0, ymax = 200, r0=0.55, r1=1, col = "blue")

